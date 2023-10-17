//
//  TriMesh.swift
//  OptCuts
//
//  Created by Reza on 6/20/23.
//

import Foundation
import Matrix
import GeometryProcessing
import Accelerate
import MetalKit

// Duplicate the vertices and edges of a mesh to separate its triangles,
// adjacent triangles in the original mesh will have a cohesive edge structure
// to indicate the connectivity
class TriMesh {
    // MARK: - Properties
    
    // Owned data
     var V_rest: Mat<Double> = .init() // duplicated rest vertex coordinates in 3D
     var V: Mat<Double> = .init() // duplicated vertex coordinates, the dimension depends on the search space
     var F: Mat<Int> = .init() // reordered triangle draw list (0, 1, 2, ...) indices based on V
     var cohE: Mat<Int> = .init() // cohesive edge pairs with the 4 end vertex indices based on V
     var initSeams: Mat<Int> = .init() // initial cohesive edge pairs actually
    
     weak var scaffold: Scaffold? = nil
     var areaThres_AM: Double = .zero // for preventing degeneracy of air mesh triangles
    
    // owned features
     var boundaryEdge: Vec<Int> = .init() // 1: boundary edge, 0: interior edge
     var edgeLen: Vec<Double> = .init() // cohesive edge rest length, used as weights
     var LaplacianMtr: SparseMatrix<Double> = .init() // 2 * V.rows() wide
     var triArea: Vec<Double> = .init() // triangle rest area
     var triNormal: Mat<Double> = .init()
     var surfaceArea: Double = .zero
     var triAreaSq: Vec<Double> = .init() // triangle rest squared area
     var e0dote1: Vec<Double> = .init() // triangle rest edge dot product
     var e0SqLen: Vec<Double> = .init()
     var e1SqLen: Vec<Double> = .init() // triangle edge rest squared length
     var e0SqLen_div_dbAreaSq: Vec<Double> = .init()
     var e1SqLen_div_dbAreaSq: Vec<Double> = .init()
     var e0dote1_div_dbAreaSq: Vec<Double> = .init()
     var avgEdgeLen: Double = .zero
     var virtualRadius: Double = .zero
     var validSplit: [Set<Pair<Int, Int>>] = []
     var fixedVert: Set<Int> = .init() // for linear solve
     var bbox: Mat<Double> = .init(2, 3) // 2x3 matrix
     var vertWeight: Vec<Double> = .init() // for regional seam placement
    
    // indices for fast access
     var edge2Tri: [Pair<Int, Int>: Int] = [:]
     var vNeighbor: [Set<Int>] = []
     var cohEIndex: [Pair<Int, Int>: Int] = [:]
    
     var fracTail: Set<Int> = .init()
     var curFracTail: Int = .zero
     var curInteriorFracTails: Pair<Int, Int> = .init(0, 0)
     var initSeamLen: Double = .zero
    
    // query split ivars
    var EwDecs: [Double] = []
    var paths_p: [[Int]] = [[]]
    var newVertPoses_p: [Matd] = []
    var energyChanges_p: [Pair<Double, Double>] = []
    
    // MARK: - Initialization
     init() {}
    
    // initialize from a triangle mesh, V will be constructed from UV_mesh in 2D,
    // V_mesh will be used to initialize restShape
     init(_ V_mesh: Mat<Double>,
         _ F_mesh: Mat<Int>,
         _ UV_mesh: Mat<Double>,
         _ FUV_mesh: Mat<Int> = Mati(),
         _ separateTri: Bool = true,
         _ p_initSeamLen: Double = 0.0,
         _ p_areaThres_AM: Double = 0.0) {
        
        initSeamLen = p_initSeamLen
        areaThres_AM = p_areaThres_AM
        
        var multiComp: Bool = false // TODO: detect whether the mesh is multi-component
        if (separateTri) {
            // duplicate vertices and edges, use new face vertex indices,
            // construct cohesive edge pairs,
            // compute triangle matrix to save rest shapes
            V_rest = Matd(F_mesh.rows * F_mesh.cols, 3)
            V = Matd(F_mesh.rows * F_mesh.cols, 2)
            F = Mati(F_mesh.rows, F_mesh.cols)
            var edge2DupInd: [Pair<Int, Int>: Vec3i] = [:]
            var cohEAmt: Int = 0
            for triI in 0..<F_mesh.rows {
                let vDupIndStart = triI * 3
                
                if (UV_mesh.rows == V_mesh.rows) {
                    // bijective map without seams, usually Tutte
                    V.row(vDupIndStart) <<== UV_mesh.row(F_mesh.row(triI)[0])
                    V.row(vDupIndStart + 1) <<== UV_mesh.row(F_mesh.row(triI)[1])
                    V.row(vDupIndStart + 2) <<== UV_mesh.row(F_mesh.row(triI)[2])
                }
                
                V_rest.row(vDupIndStart) <<== V_mesh.row(F_mesh.row(triI)[0])
                V_rest.row(vDupIndStart + 1) <<== V_mesh.row(F_mesh.row(triI)[1])
                V_rest.row(vDupIndStart + 2) <<== V_mesh.row(F_mesh.row(triI)[2])
                
                F[triI, 0] = vDupIndStart
                F[triI, 1] = vDupIndStart + 1
                F[triI, 2] = vDupIndStart + 2
                
                for vI in 0..<3 {
                    let vsI = F_mesh.row(triI)[vI]
                    let veI = F_mesh.row(triI)[(vI + 1) % 3]
                    
                    if let cohE = edge2DupInd[Pair(veI, vsI)] {
                        edge2DupInd[Pair(vsI, veI)] = Vec3i(cohE[0], F[triI, vI], F[triI, (vI + 1) % 3])
                    } else {
                        cohEAmt += 1
                        edge2DupInd[Pair(vsI, veI)] = Vec3i(cohEAmt, F[triI, vI], F[triI, (vI + 1) % 3])
                    }
                }
            }
            
            cohE = Mati(cohEAmt, 4)
            cohE.setConstant(-1)
            
            for (_, vec) in edge2DupInd {
                if vec[0] > 0 {
                    cohE.row(vec[0] - 1).valuesPtr.pointer[0] = vec[1]
                    cohE.row(vec[0] - 1).valuesPtr.pointer[1] = vec[2]
                } else {
                    cohE.row(-vec[0] - 1).valuesPtr.pointer[2] = vec[2]
                    cohE.row(-vec[0] - 1).valuesPtr.pointer[3] = vec[1]
                }
            }
            
            if (UV_mesh.rows == 0) {
                // no input UV
                initRigidUV()
            } else if (UV_mesh.rows != V_mesh.rows) {
                // input UV with seams
                assert(false, "TODO: separate each triangle in UV space according to FUV!")
            }
        } else {
            // deal with mesh
            if (UV_mesh.rows == V_mesh.rows) {
                // same vertex and UV index
                V_rest = .init(V_mesh, V_mesh.rows, V_mesh.cols)
                V = .init(UV_mesh, UV_mesh.rows, UV_mesh.cols)
                F = .init(F_mesh, F_mesh.rows, F_mesh.cols)
            } else if (UV_mesh.rows != 0) {
                // different vertex and uv index, split 3D surface according to UV and merge back while saving into files
                assert(F_mesh.rows == FUV_mesh.rows)
                // UV map contains seams
                // Split triangles along the seams on the surface (construct cohesive edges there)
                // to construct a bijective map
                var HE_UV: Set<Pair<Int, Int>> = .init()
                var HE: [Pair<Int, Int>: Pair<Int, Int>] = [:]
                for triI in 0..<FUV_mesh.rows {
                    let triVInd_UV: RVec3i = FUV_mesh.row(triI)
                    HE_UV.insert(Pair(triVInd_UV[0], triVInd_UV[1]))
                    HE_UV.insert(Pair(triVInd_UV[1], triVInd_UV[2]))
                    HE_UV.insert(Pair(triVInd_UV[2], triVInd_UV[0]))
                    let triVInd: RVec3i = F_mesh.row(triI)
                    HE[Pair(triVInd[0], triVInd[1])] = Pair(triI, 0)
                    HE[Pair(triVInd[1], triVInd[2])] = Pair(triI, 1)
                    HE[Pair(triVInd[2], triVInd[0])] = Pair(triI, 2)
                }
                var cohEdges: [[Int]] = []
                for triI in 0..<FUV_mesh.rows {
                    let triVInd_UV: RVec3i = FUV_mesh.row(triI)
                    let triVInd: RVec3i = F_mesh.row(triI)
                    for eI in 0..<3 {
                        let vI = eI
                        let vI_post: Int = (eI + 1) % 3
                        if !HE_UV.contains(Pair(triVInd_UV[vI_post], triVInd_UV[vI])) {
                            // boundary edge in UV space
                            if let pair = HE[Pair(triVInd[vI_post], triVInd[vI])] {
                                // non - boundary edge on the surface
                                // construct cohesive edge pair
                                cohEdges.append([])
                                cohEdges[cohEdges.endIndex - 1].append(triVInd_UV[vI])
                                cohEdges[cohEdges.endIndex - 1].append(triVInd_UV[vI_post])
                                cohEdges[cohEdges.endIndex - 1].append(FUV_mesh[pair.first, (pair.second + 1) % 3])
                                cohEdges[cohEdges.endIndex - 1].append(FUV_mesh[pair.first, pair.second])
                                HE.removeValue(forKey: Pair(triVInd[vI], triVInd[vI_post])) // prevent from inserting again
                            }
                        }
                    }
                }
                let makeCoh: Bool = true
                if (makeCoh) {
                    list_to_matrix(cohEdges, &cohE)
                }
                
                V_rest.resize(UV_mesh.rows, 3)
                V = .init(UV_mesh, UV_mesh.rows, UV_mesh.cols)
                F = .init(FUV_mesh, FUV_mesh.rows, FUV_mesh.cols)
                var updated: [Bool] = .init(repeating: false, count: UV_mesh.rows)
                for triI in 0..<F_mesh.rows {
                    let triVInd: RVec3i = F_mesh.row(triI)
                    let triVInd_UV: RVec3i = FUV_mesh.row(triI)
                    for vI in 0..<3 {
                        if (!updated[triVInd_UV[vI]]) {
                            V_rest.row(triVInd_UV[vI]) <<== V_mesh.row(triVInd[vI])
                            updated[triVInd_UV[vI]] = true
                        }
                    }
                }
                
                if (!makeCoh) {
                    for cohI in cohEdges {
                        initSeamLen += (V_rest.row(cohI[0]) - V_rest.row(cohI[1])).norm()
                    }
                } else {
                    initSeams = cohE
                }
            } else {
                assert(V_mesh.rows > 0)
                assert(F_mesh.rows > 0)
                V_rest = .init(V_mesh, V_mesh.rows, V_mesh.cols)
                F = .init(F_mesh, F_mesh.rows, F_mesh.cols)
                V = .Zero(V_rest.rows, 2)
                print("No UV provided, initialized to all 0")
            }
        }
        computeFeatures(false, true)
        vertWeight = .Ones(V.rows)
    }
    
     convenience init(primitive: Primitive, size: Double = 1.0, spacing: Double = 0.1, separateTri: Bool = true) {
        var V_rest = Mat<Double>()
        var V = Mat<Double>()
        var F = Mat<Int>()
        
        switch primitive {
        case .P_SQUARE:
            assert(size >= spacing)
            let gridSize: Int = Int(size / spacing) + 1
            var spacing = size / Double(gridSize - 1)
            V_rest.resize(gridSize * gridSize, 3)
            V.resize(gridSize * gridSize, 2)
            for rowI in 0..<gridSize {
                for colI in 0..<gridSize {
                    V_rest.row(rowI * gridSize + colI) <<== Vec3d(spacing * Double(colI), spacing * Double(rowI), 0.0)
                    V.row(rowI * gridSize + colI) <<== (spacing * Vec2d(colI, rowI))
                }
            }
            
            F.resize((gridSize - 1) * (gridSize - 1) * 2, 3)
            for rowI in 0..<(gridSize - 1) {
                for colI in 0..<(gridSize - 1) {
                    let squareI: Int = rowI * (gridSize - 1) + colI
                    F.row(squareI * 2) <<== Vec3i(
                        rowI * gridSize + colI, (rowI + 1) * gridSize + colI + 1, (rowI + 1) * gridSize + colI)
                    F.row(squareI * 2 + 1) <<== Vec3i(
                        rowI * gridSize + colI, rowI * gridSize + colI + 1, (rowI + 1) * gridSize + colI + 1)
                }
            }
            break
            
        case .P_CYLINDER:
            withUnsafeMutablePointer(to: &V) { vPtr in
                initCylinder(0.5, 0.5, 1.0, 1.0, 1.0, 20, 20, &V_rest, &F, vPtr)
            }
            break
        }
        
        if (separateTri) {
            self.init(V_rest, F, V)
        } else {
            self.init()
            computeFeatures(false, true)
        }
        initSeamLen = 0.0
    }
    
    // MARK: - Methods
     func saveAsMesh(F0: Mati, scaleUV: Bool) -> Data {
        assert(F0.rows == F.rows)
        assert(F0.cols == 3)
        
        var V_mesh: Matd = .init()
        for fI in 0..<F0.rows {
            for localVI in 0..<3 {
                let vI = F0[fI, localVI]
                if vI >= V_mesh.rows {
                    V_mesh.conservativeResize(vI + 1, 3)
                }
                V_mesh.row(vI) <<== V_rest.row(F[fI, localVI])
            }
        }
        
        var UV_mesh: Matd = .init(V, V.rows, V.cols)
        
        if scaleUV {
            let u: Vec<Double> = UV_mesh.col(0)
            let v: Vec<Double> = UV_mesh.col(1)
            let uMin = u.minCoeff()
            let vMin = v.minCoeff()
            let uScale = u.maxCoeff() - uMin
            let vScale = v.maxCoeff() - vMin
            let scale = max(uScale, vScale)
            for uvI in 0..<UV_mesh.rows {
                UV_mesh[uvI, 0] = (UV_mesh[uvI, 0] - uMin) / scale
                UV_mesh[uvI, 1] = (UV_mesh[uvI, 1] - vMin) / scale
            }
        }
        
        return save(V: V_mesh, F: F0, UV: UV_mesh, FUV: F)
    }
    
     func save(V: Matd, F: Mati, UV: Matd, FUV: Mati) -> Data {
        let formatter = NumberFormatter()
        formatter.numberStyle = .scientific
        formatter.maximumFractionDigits = 6
        
        var out: String = .init()
        for vI in 0..<V.rows {
            let v: RVec3d = V.row(vI)
            out += "v "
            out += formatter.string(for: v[0])!
            out += " "
            out += formatter.string(for: v[1])!
            out += " "
            out += formatter.string(for: v[2])!
            out += "\n"
        }
        
        for vI in 0..<UV.rows {
            let uv: RVec2d = UV.row(vI)
            out += "vt "
            out += formatter.string(for: uv[0])!
            out += " "
            out += formatter.string(for: uv[1])!
            out += "\n"
        }
        
        formatter.numberStyle = .none
        for triI in 0..<F.rows {
            let tri: RVec3i = F.row(triI)
            let tri_UV: RVec3i = FUV.row(triI)
            out += "f "
            out += formatter.string(for: tri[0] + 1)!
            out += "/"
            out += formatter.string(for: tri_UV[0] + 1)!
            out += " "
            out += formatter.string(for: tri[1] + 1)!
            out += "/"
            out += formatter.string(for: tri_UV[1] + 1)!
            out += " "
            out += formatter.string(for: tri[2] + 1)!
            out += "/"
            out += formatter.string(for: tri_UV[2] + 1)!
            out += "\n"
        }
        
        let data = out.data(using: .utf8)!
        return data
    }
    
     func clone() -> TriMesh {
        var mesh = TriMesh()
        mesh.V_rest = .init(V_rest, V_rest.rows, V_rest.cols)
        mesh.V = .init(V, V.rows, V.cols)
        mesh.F = .init(F, F.rows, F.cols)
        mesh.cohE = .init(cohE, cohE.rows, cohE.cols)
        mesh.initSeams = .init(initSeams, initSeams.rows, initSeams.cols)
        mesh.scaffold = scaffold
        mesh.areaThres_AM = areaThres_AM
        
        mesh.boundaryEdge = .init(boundaryEdge, boundaryEdge.rows, boundaryEdge.cols)
        mesh.edgeLen = .init(edgeLen, edgeLen.rows, edgeLen.cols)
        mesh.LaplacianMtr = LaplacianMtr
        mesh.triArea = .init(triArea, triArea.rows, triArea.cols)
        mesh.triNormal = .init(triNormal, triNormal.rows, triNormal.cols)
        mesh.surfaceArea = surfaceArea
        mesh.triAreaSq = .init(triAreaSq, triAreaSq.rows, triAreaSq.cols)
        mesh.e0dote1 = .init(e0dote1, e0dote1.rows, e0dote1.cols)
        mesh.e0SqLen = .init(e0SqLen, e0SqLen.rows, e0SqLen.cols)
        mesh.e1SqLen = .init(e1SqLen, e1SqLen.rows, e1SqLen.cols)
        mesh.e0SqLen_div_dbAreaSq = .init(e0SqLen_div_dbAreaSq, e0SqLen_div_dbAreaSq.rows, e0SqLen_div_dbAreaSq.cols)
        mesh.e1SqLen_div_dbAreaSq = .init(e1SqLen_div_dbAreaSq, e1SqLen_div_dbAreaSq.rows, e1SqLen_div_dbAreaSq.cols)
        mesh.e0dote1_div_dbAreaSq = .init(e0dote1_div_dbAreaSq, e0dote1_div_dbAreaSq.rows, e0dote1_div_dbAreaSq.cols)
        mesh.avgEdgeLen = avgEdgeLen
        mesh.virtualRadius = virtualRadius
        mesh.validSplit = validSplit
        mesh.fixedVert = fixedVert
        mesh.bbox = .init(bbox, bbox.rows, bbox.cols)
        mesh.vertWeight = .init(vertWeight, vertWeight.rows, vertWeight.cols)
        
        mesh.edge2Tri = edge2Tri
        mesh.vNeighbor = vNeighbor
        mesh.cohEIndex = cohEIndex
        mesh.fracTail = fracTail
        mesh.curFracTail = curFracTail
        mesh.curInteriorFracTails = curInteriorFracTails
        mesh.initSeamLen = initSeamLen
        
        return mesh
    }
    
     func computeFeatures(_ multiComp: Bool = false, _ resetFixedV: Bool = false) {
        if (resetFixedV) {
            fixedVert.removeAll()
            fixedVert.insert(0)
        }
        
        boundaryEdge.resize(cohE.rows)
        edgeLen.resize(cohE.rows)
        for cohI in 0..<cohE.rows {
            if (cohE.row(cohI).minCoeff() >= 0) {
                boundaryEdge[cohI] = 0
            } else {
                boundaryEdge[cohI] = 1
            }
            edgeLen[cohI] = (V_rest.row(cohE[cohI, 0]) - V_rest.row(cohE[cohI, 1])).norm()
        }
        
        triNormal.resize(F.rows, 3)
        triArea.resize(F.rows)
        surfaceArea = 0.0
        triAreaSq.resize(F.rows)
        e0SqLen.resize(F.rows)
        e1SqLen.resize(F.rows)
        e0dote1.resize(F.rows)
        e0SqLen_div_dbAreaSq.resize(F.rows)
        e1SqLen_div_dbAreaSq.resize(F.rows)
        e0dote1_div_dbAreaSq.resize(F.rows)
        var vertNormals: [RVec3d] = .init(repeating: .init(), count: V_rest.rows)
        var isMeshInvalid: Bool = false
        for triI in 0..<F.rows {
            let triVInd: Vec3i = F.row(triI)
            
            let P1: Vec3d = V_rest.row(triVInd[0])
            let P2: Vec3d = V_rest.row(triVInd[1])
            let P3: Vec3d = V_rest.row(triVInd[2])
            
            let P2m1: Vec3d = P2 - P1
            let P3m1: Vec3d = P3 - P1
            let normalVec: RVec3d = P2m1.cross(P3m1)
            
            triNormal.row(triI) <<== normalVec.normalized()
            triArea[triI] = 0.5 * normalVec.norm()
            if (triArea[triI] == 0.0) {
                print("area of triangle No. \(triI) (0-indx) is 0")
                print("vertex indices (0-index) are \(triVInd.transpose())")
                isMeshInvalid = true
            }
            if (triArea[triI] < areaThres_AM) {
                // air mesh triangle degeneracy prevention
                triArea[triI] = areaThres_AM
                surfaceArea += areaThres_AM
                triAreaSq[triI] = areaThres_AM * areaThres_AM
                e1SqLen[triI] = 4.0 / sqrt(3.0) * areaThres_AM // TODO: Possible bug (operator order)
                e0SqLen[triI] = e1SqLen[triI]
                e0dote1[triI] = e0SqLen[triI] / 2.0
                
                e1SqLen_div_dbAreaSq[triI] = 2.0 / sqrt(3.0) / areaThres_AM
                e0SqLen_div_dbAreaSq[triI] = e1SqLen_div_dbAreaSq[triI]
                e0dote1_div_dbAreaSq[triI] = e0SqLen_div_dbAreaSq[triI] / 2.0
            } else {
                surfaceArea += triArea[triI]
                triAreaSq[triI] = triArea[triI] * triArea[triI]
                e0SqLen[triI] = P2m1.squaredNorm()
                e1SqLen[triI] = P3m1.squaredNorm()
                e0dote1[triI] = P2m1.dot(P3m1)
                
                e0SqLen_div_dbAreaSq[triI] = e0SqLen[triI] / 2.0 / triAreaSq[triI]
                e1SqLen_div_dbAreaSq[triI] = e1SqLen[triI] / 2.0 / triAreaSq[triI]
                e0dote1_div_dbAreaSq[triI] = e0dote1[triI] / 2.0 / triAreaSq[triI]
            }
            vertNormals[triVInd[0]] += normalVec
            vertNormals[triVInd[1]] += normalVec
            vertNormals[triVInd[2]] += normalVec
        }
        if (isMeshInvalid) {
            fatalError("please clean up the mesh and retry.")
        }
        avgEdgeLen = avg_edge_length(V_rest, F)
        virtualRadius = sqrt(surfaceArea / .pi)
        for i in 0..<vertNormals.count {
            vertNormals[i].normalize()
        }
        
        //computeLabpaciamMtr()
        
        bbox.block(0, 0, 1, 3) <<== V_rest.row(0)
        bbox.block(1, 0, 1, 3) <<== V_rest.row(0)
        for vI in 1..<V_rest.rows {
            let v: RVec3d = V_rest.row(vI)
            for dimI in 0..<3 {
                if (v[dimI] < bbox[0, dimI]) {
                    bbox[0, dimI] = v[dimI]
                }
                if (v[dimI] > bbox[1, dimI]) {
                    bbox[1, dimI] = v[dimI]
                }
            }
        }
        
        edge2Tri.removeAll()
        vNeighbor.removeAll()
        vNeighbor = .init(repeating: Set<Int>(), count: V_rest.rows)
        for triI in 0..<F.rows {
            let triVInd: RVec3i = F.row(triI)
            for vI in 0..<3 {
                let vI_post: Int = (vI + 1) % 3
                edge2Tri[Pair(triVInd[vI], triVInd[vI_post])] = triI
                vNeighbor[triVInd[vI]].insert(triVInd[vI_post])
                vNeighbor[triVInd[vI_post]].insert(triVInd[vI])
            }
        }
        cohEIndex.removeAll()
        for cohI in 0..<cohE.rows {
            let cohEI: RVec4i = cohE.row(cohI)
            if (cohEI.minCoeff() >= 0) {
                cohEIndex[Pair(cohEI[0], cohEI[1])] = cohI
                cohEIndex[Pair(cohEI[3], cohEI[2])] = -cohI - 1
            }
        }
                
        validSplit = .init(repeating: Set(), count: V_rest.rows)
        for vI in 0..<V_rest.rows {
            validSplit[vI].removeAll()
            if (isBoundaryVert(vI)) {
                continue
            }
            
            let nbVs: [Int] = .init(vNeighbor[vI])
            //let projectedEdge: UnsafeMutablePointer<RVec3d> = .allocate(capacity: nbVs.count)
            var projectedEdge: [RVec3d?] = .init(repeating: nil, count: nbVs.count)
            /*
            defer {
                projectedEdge.deallocate()
            }*/
            for nbI in 0..<nbVs.count {
                let edge: RVec3d = V_rest.row(nbVs[nbI]) - V_rest.row(vI)
                projectedEdge[nbI] = (edge - edge.dot(vertNormals[vI]) * vertNormals[vI]).normalized()
            }
            
            for nbI in 0..<(nbVs.count - 1) {
                for nbJ in (nbI + 1)..<nbVs.count {
                    if projectedEdge[nbI]!.dot(projectedEdge[nbJ]!) <= 0.0 {
                        validSplit[vI].insert(Pair(nbVs[nbI], nbVs[nbJ]))
                        validSplit[vI].insert(Pair(nbVs[nbJ], nbVs[nbI]))
                    }
                }
            }
        }
        // init fracture tail record
        fracTail.removeAll()
        for cohI in 0..<cohE.rows {
            if (cohE[cohI, 0] == cohE[cohI, 2]) {
                fracTail.insert(cohE[cohI, 0])
            } else if (cohE[cohI, 1] == cohE[cohI, 3]) {
                fracTail.insert(cohE[cohI, 1])
            }
        }
        
        // tails of initial seams doesn't count as fracture tails for propogation
        for initSeamI in 0..<initSeams.rows {
            if (initSeams[initSeamI, 0] == initSeams[initSeamI, 2]) {
                fracTail.remove(initSeams[initSeamI, 0])
            } else if (initSeams[initSeamI, 1] == initSeams[initSeamI, 3]) {
                fracTail.remove(initSeams[initSeamI, 1])
            }
        }
    }
    
     func updateFeatures() {
        let nCE: Int = boundaryEdge.count
        boundaryEdge.conservativeResize(cohE.rows)
        edgeLen.conservativeResize(cohE.rows)
        for cohI in nCE..<cohE.rows {
            if (cohE.row(cohI).minCoeff() >= 0) {
                boundaryEdge[cohI] = 0
            } else {
                boundaryEdge[cohI] = 1
            }
            edgeLen[cohI] = (V_rest.row(cohE[cohI, 0]) - V_rest.row(cohE[cohI, 1])).norm()
        }
        computeLabpaciamMtr()
    }
    
     func resetFixedVert(_ p_fixedVert: Set<Int>) {
        for vI in p_fixedVert {
            assert(vI < V.rows)
        }
        
        fixedVert = p_fixedVert
        computeLabpaciamMtr()
    }
    
     func buildCohEfromRecord(_ cohERecord: Mat<Int>) {
        assert(cohERecord.cols == 4)
        
        cohE.resize(cohERecord.rows, 4)
        for cohEI in 0..<cohERecord.rows {
            let triI1 = cohERecord[cohEI, 0]
            let vI11 = F[triI1, cohERecord[cohEI, 1]]
            let vI12 = F[triI1, (cohERecord[cohEI, 1] + 1) % 3]
            let triI2 = cohERecord[cohEI, 2]
            let vI21 = F[triI2, cohERecord[cohEI, 3]]
            let vI22 = F[triI2, (cohERecord[cohEI, 3] + 1) % 3]
            cohE.row(cohEI) <<== [vI11, vI12, vI22, vI21]
        }
        
        computeFeatures()
    }
    
     func querySplit(_ lambda_t: Double,
                           _ propogate: Bool,
                           _ splitInterior: Bool,
                           _ EwDec_max: inout Double,
                           _ path_max: inout [Int],
                           _ newVertPos_max: inout Mat<Double>,
                           _ energyChanges_max: inout Pair<Double, Double>) {
        
        var splitInterior = splitInterior
        let filterExp_b: Double = 0.8
        let filterMult_b: Double = 1.0 // TODO: Better use ratio
        
        var bestCandVerts: [Int] = []
        if (!propogate) {
            let SD = SymDirichletEnergy()
            var divGradPerVert = Vec<Double>()
            SD.computeDivGradPerVert(self, &divGradPerVert)
            
            var candVerts_b: [Double: Int] = [:]
            var candVerts_in: [Double: Int] = [:]
            if (splitInterior) {
                for vI in 0..<V_rest.rows {
                    if (vNeighbor[vI].count <= 2) {
                        // this vertex is impossible to be splitted further
                        continue
                    }
                    
                    if (!isBoundaryVert(vI)) {
                        var connectToBound: Bool = false
                        for nbVI in vNeighbor[vI] {
                            if (isBoundaryVert(nbVI)) {
                                connectToBound = true
                                break
                            }
                        }
                        if (!connectToBound) {
                            // don't split vertices connected to boundary here
                            candVerts_in[-divGradPerVert[vI] / vertWeight[vI]] = vI
                        }
                    }
                }
                inSplitTotalAmt = candVerts_in.count
            } else {
                for vI in 0..<V_rest.rows {
                    if (vNeighbor[vI].count <= 2) {
                        // this vertex is impossible to be splitted further
                        continue
                    }
                    
                    if (isBoundaryVert(vI)) {
                        candVerts_b[-divGradPerVert[vI] / vertWeight[vI]] = vI
                    }
                }
            }
            
            let sortedCandVerts_b = candVerts_b.sorted(by: {$0.key < $1.key })
            let sortedCandVerts_in = candVerts_in.sorted(by: {$0.key < $1.key })
            
            if (!splitInterior) {
                if (sortedCandVerts_b.isEmpty) {
                    EwDec_max = 0.0
                    return
                }
                var bestCandAmt_b = Int(pow(Double(sortedCandVerts_b.count), filterExp_b) * filterMult_b)
                if (bestCandAmt_b < 2){
                    bestCandAmt_b = 2
                }
                bestCandVerts.reserveCapacity(bestCandAmt_b)
                for (_, value) in sortedCandVerts_b {
                    bestCandVerts.append(value)
                    if (bestCandVerts.count >= bestCandAmt_b) {
                        break
                    }
                }
            } else {
                if (sortedCandVerts_in.isEmpty) {
                    EwDec_max = 0.0
                    return
                }
                var bestCandAmt_in = Int(pow(Double(sortedCandVerts_in.count), filterExp_in))
                if (bestCandAmt_in < 2) {
                    bestCandAmt_in = 2
                }
                bestCandVerts.reserveCapacity(bestCandVerts.count + bestCandAmt_in)
                for (_, value) in sortedCandVerts_in {
                    bestCandVerts.append(value)
                    if (bestCandVerts.count >= bestCandAmt_in) {
                        break
                    }
                }
            }
        } else {
            // see whether fracture could be propogated from each fracture tail
            if (curFracTail < 0) {
                if (curInteriorFracTails.first < 0) {
                    EwDec_max = -__DBL_MAX__
                    path_max.removeAll()
                    newVertPos_max.resize(0, 2)
                    return
                } else {
                    assert(curInteriorFracTails.second >= 0)
                    splitInterior = false
                    bestCandVerts.append(curInteriorFracTails.first)
                    bestCandVerts.append(curInteriorFracTails.second)
                }
            } else {
                splitInterior = false
                bestCandVerts.append(curFracTail)
            }
        }
        
        assert(!bestCandVerts.isEmpty)
        
        // evaluate local energy decrease
        print("evaluate vertex splits, \(bestCandVerts.count) candidate verts")
        // run in parallel
        EwDecs.conservativeResize(to: bestCandVerts.count)
        var operationType: Int = -1
        // query boundary splits
        if (!splitInterior) {
            if (propogate) {
                paths_p = .init(repeating: [], count: bestCandVerts.count)
                newVertPoses_p = .init(repeating: .init(), count: bestCandVerts.count)
                energyChanges_p = .init(repeating: .init(0, 0), count: bestCandVerts.count)
                DispatchQueue.concurrentPerform(iterations: bestCandVerts.count) { [self] candI in
                    EwDecs[candI] = computeLocalLDec(bestCandVerts[candI],
                                                     lambda_t,
                                                     &paths_p[candI],
                                                     &newVertPoses_p[candI],
                                                     &energyChanges_p[candI])
                }
                /*
                for candI in 0..<bestCandVerts.count {
                    EwDecs[candI] = computeLocalLDec(bestCandVerts[candI],
                                                     lambda_t,
                                                     &paths_p[candI],
                                                     &newVertPoses_p[candI],
                                                     &energyChanges_p[candI])
                }*/
            } else {
                operationType = 0
                paths_bSplit = .init(repeating: [], count: bestCandVerts.count)
                //newVertPoses_bSplit = .init(repeating: .init(), count: bestCandVerts.count)
                newVertPoses_bSplit.conservativeResize(to: bestCandVerts.count)
                //energyChanges_bSplit = .init(repeating: .init(0, 0), count: bestCandVerts.count)
                energyChanges_bSplit.conservativeResize(to: bestCandVerts.count)
                DispatchQueue.concurrentPerform(iterations: bestCandVerts.count) { [self] candI in
                    EwDecs[candI] = computeLocalLDec(bestCandVerts[candI],
                                                     lambda_t,
                                                     &paths_bSplit[candI],
                                                     &newVertPoses_bSplit[candI],
                                                     &energyChanges_bSplit[candI])
                }
                /*
                for candI in 0..<bestCandVerts.count {
                    EwDecs[candI] = computeLocalLDec(bestCandVerts[candI],
                                                     lambda_t,
                                                     &paths_bSplit[candI],
                                                     &newVertPoses_bSplit[candI],
                                                     &energyChanges_bSplit[candI])
                }*/
            }
        } else {
            assert(!propogate)
            operationType = 1
            // query interior splits
            paths_iSplit = .init(repeating: [], count: bestCandVerts.count)
            //newVertPoses_iSplit = .init(repeating: .init(), count: bestCandVerts.count)
            newVertPoses_iSplit.conservativeResize(to: bestCandVerts.count)
            //energyChanges_iSplit = .init(repeating: .init(0, 0), count: bestCandVerts.count)
            energyChanges_iSplit.conservativeResize(to: bestCandVerts.count)
            DispatchQueue.concurrentPerform(iterations: bestCandVerts.count) { [self] candI in
                EwDecs[candI] = computeLocalLDec(bestCandVerts[candI],
                                                 lambda_t,
                                                 &paths_iSplit[candI],
                                                 &newVertPoses_iSplit[candI],
                                                 &energyChanges_iSplit[candI])
                if (EwDecs[candI] != -__DBL_MAX__) {
                    EwDecs[candI] *= 0.5
                }
            }
            /*
            for candI in 0..<bestCandVerts.count {
                EwDecs[candI] = computeLocalLDec(bestCandVerts[candI],
                                                 lambda_t,
                                                 &paths_iSplit[candI],
                                                 &newVertPoses_iSplit[candI],
                                                 &energyChanges_iSplit[candI])
                if (EwDecs[candI] != -__DBL_MAX__) {
                    EwDecs[candI] *= 0.5
                }
            }*/
        }
        
        var candI_max = 0
        for candI in 1..<bestCandVerts.count {
            if (EwDecs[candI] > EwDecs[candI_max]) {
                candI_max = candI
            }
        }
        
        EwDec_max = EwDecs[candI_max]
        switch (operationType) {
        case -1:
            path_max = paths_p[candI_max]
            newVertPos_max = .init(newVertPoses_p[candI_max], newVertPoses_p[candI_max].rows, newVertPoses_p[candI_max].cols)
            energyChanges_max = energyChanges_p[candI_max]
            
        case 0:
            path_max = paths_bSplit[candI_max]
            newVertPos_max = .init(newVertPoses_bSplit[candI_max], newVertPoses_bSplit[candI_max].rows, newVertPoses_bSplit[candI_max].cols)
            energyChanges_max = energyChanges_bSplit[candI_max]
            
        case 1:
            path_max = paths_iSplit[candI_max]
            newVertPos_max = .init(newVertPoses_iSplit[candI_max], newVertPoses_iSplit[candI_max].rows, newVertPoses_iSplit[candI_max].cols)
            energyChanges_max = energyChanges_iSplit[candI_max]
            
        default:
            assert(false)
        }
    }
    
     func splitEdge(_ lambda_t: Double,
                          _ thres: Double = 0.0,
                          _ propogate: Bool = false,
                          _ splitInterior: Bool = false) -> Bool {
        var EwDec_max: Double = .zero
        var path_max: [Int] = []
        var newVertPos_max: Matd = .init()
        var energyChanges_max: Pair<Double, Double> = .init(0.0, 0.0)
        querySplit(lambda_t, propogate, splitInterior,
                   &EwDec_max, &path_max, &newVertPos_max,
                   &energyChanges_max)
        print("E_dec threshold = \(thres)")
        if (EwDec_max > thres) {
            if (!splitInterior) {
                // boundary split
                print("boundary split E_dec = \(EwDec_max)")
                splitEdgeOnBoundary(Pair(path_max[0], path_max[1]), newVertPos_max)
                // TODO: process fractail here!
                updateFeatures()
            } else {
                // interior split
                assert(!propogate)
                print("interior split E_dec = \(EwDec_max)")
                cutPath(path_max, true, 1, newVertPos_max)
                print("interior edge splitted")
                fracTail.insert(path_max[0])
                fracTail.insert(path_max[2])
                curInteriorFracTails.first = path_max[0]
                curInteriorFracTails.second = path_max[2]
                curFracTail = -1
            }
            return true
        } else {
            print("max E_dec = \(EwDec_max) thres \(thres)")
            return false
        }
    }
    
     func queryMerge(_ lambda: Double,
                           _ propogate: Bool,
                           _ localEwDec_max: inout Double,
                           _ path_max: inout [Int],
                           _ newVertPos_max: inout Mat<Double>,
                           _ energyChanges_max: inout Pair<Double, Double>) {
        // TODO: local index updates in mergeBoundaryEdge()
        // TODO: parallelize the query
        print("evaluate edgemerge, \(cohE.rows) cohesive edge pairs.")
        
        localEwDec_max = -__DBL_MAX__
        
        if (!propogate) {
            paths_merge.removeAll()
            newVertPoses_merge.removeAll()
            energyChanges_merge.removeAll()
        }
        
        for cohI in 0..<cohE.rows {
            var forkVI: Int = 0
            if (cohE[cohI, 0] == cohE[cohI, 2]) {
                forkVI = 0
            } else if (cohE[cohI, 1] == cohE[cohI, 3]) {
                forkVI = 1
            } else {
                // only consider "zipper bottom" edge pairs for now
                continue
            }
            
            if (propogate) {
                if (cohE[cohI, forkVI] != curFracTail) {
                    continue
                }
            }
            
            // find incident triangles for inversion check and local energy decrease evaluation
            var triangles: [Int] = []
            var firstVertIncTriAmt: Int = 0
            for mergeVI in [1, 3] {
                for nbVI in vNeighbor[cohE[cohI, mergeVI - forkVI]] {
                    if let value = edge2Tri[Pair(cohE[cohI, mergeVI - forkVI], nbVI)] {
                        triangles.append(value)
                    }
                }
                if (mergeVI == 1) {
                    firstVertIncTriAmt = triangles.count
                    assert(firstVertIncTriAmt >= 1)
                }
            }
            
            var mergedPos: RVec2d = (V.row(cohE[cohI, 1 - forkVI]) + V.row(cohE[cohI, 3 - forkVI])) / 2.0
            let backup0: RVec2d = V.row(cohE[cohI, 1 - forkVI])
            let backup1: RVec2d = V.row(cohE[cohI, 3 - forkVI])
            V.row(cohE[cohI, 1 - forkVI]) <<== mergedPos
            V.row(cohE[cohI, 3 - forkVI]) <<== mergedPos
            
            if (checkInversion(true, triangles)) {
                V.row(cohE[cohI, 1 - forkVI]) <<== backup0
                V.row(cohE[cohI, 3 - forkVI]) <<== backup1
            } else {
                // project mergedPos to feasible set via Relaxation method for linear inequalities
                
                // find inequality constraints by opposite edge in incident triangles
                var inequalityConsMtr = Matd()
                var inequalityConsVec = Vec<Double>()
                for triII in 0..<triangles.count {
                    let triI: Int = triangles[triII]
                    let vI_toMerge: Int = ((triII < firstVertIncTriAmt) ? cohE[cohI, 1 - forkVI] : cohE[cohI, 3 - forkVI])
                    for i in 0..<3 {
                        if (F[triI, i] == vI_toMerge) {
                            let v1: RVec2d = V.row(F[triI, (i + 1) % 3])
                            let v2: RVec2d = V.row(F[triI, (i + 2) % 3])
                            let coef = RVec2d(v2[1] - v1[1], v1[0] - v2[0])
                            inequalityConsMtr.conservativeResize(inequalityConsMtr.rows + 1, 2)
                            inequalityConsMtr.row(inequalityConsMtr.rows - 1) <<== (coef / coef.norm())
                            inequalityConsVec.conservativeResize(inequalityConsVec.count + 1)
                            inequalityConsVec[inequalityConsVec.count - 1] = (v1[0] * v2[1] - v1[1] * v2[0]) / coef.norm()
                            break
                        }
                    }
                }
                assert(inequalityConsMtr.rows == triangles.count)
                assert(inequalityConsVec.count == triangles.count)
                
                // Relaxation method for linear inequalities
                let maxIter: Int = 70
                let eps_IC: Double = 1.0e-6 * avgEdgeLen
                for _ in 0..<maxIter {
                    var maxRes: Double = -__DBL_MAX__
                    for consI in 0..<inequalityConsMtr.rows {
                        var res: Double = inequalityConsMtr.row(consI).dot(mergedPos) - inequalityConsVec[consI]
                        if (res > -eps_IC) {
                            // project
                            mergedPos -= ((res + eps_IC) * inequalityConsMtr.row(consI))
                        }
                        if (res > maxRes) {
                            maxRes = res
                        }
                    }
                    
                    if (maxRes < 0.0) {
                        // converged (non-inversion satisfied)
                        // NOTE: although this maxRes is 1 iteration begind, it is OK for a convergence check
                        break
                    }
                }
                
                V.row(cohE[cohI, 1 - forkVI]) <<== mergedPos
                V.row(cohE[cohI, 3 - forkVI]) <<== mergedPos
                var noInversion: Bool = checkInversion(true, triangles)
                V.row(cohE[cohI, 1 - forkVI]) <<== backup0
                V.row(cohE[cohI, 3 - forkVI]) <<== backup1
                
                if (!noInversion) {
                    // because propogation is not at E_SD stationary, so it's possible to have no feasible region
                    continue
                }
            }
            
            // optimize local distortion
            var path: [Int] = []
            if (forkVI != 0) {
                path.append(contentsOf: [cohE[cohI, 0], cohE[cohI, 1], cohE[cohI, 2]])
            } else {
                path.append(contentsOf: [cohE[cohI, 3], cohE[cohI, 2], cohE[cohI, 1]])
            }
            var newVertPos: Matd = .init()
            var energyChanges: Pair<Double, Double> = .init(0, 0)
            let localEwDec = computeLocalLDec(0, lambda, &path, &newVertPos, &energyChanges, triangles, mergedPos)
            
            if (!propogate) {
                paths_merge.append(path)
                newVertPoses_merge.append(newVertPos)
                energyChanges_merge.append(energyChanges)
            }
            
            if (localEwDec > localEwDec_max) {
                localEwDec_max = localEwDec
                newVertPos_max = .init(newVertPos, newVertPos.rows, newVertPos.cols)
                path_max = path
                energyChanges_max = energyChanges
            }
        }
    }
    
     func mergeEdge(_ lambda: Double,
                   _ EDecThres: Double,
                   _ propogate: Bool) -> Bool {
        var localEwDec_max: Double = .zero
        var path_max: [Int] = []
         var newVertPos_max: Matd = .init()
         var energyChanges_max: Pair<Double, Double> = .init(.zero, .zero)
        queryMerge(lambda, propogate, &localEwDec_max, &path_max, &newVertPos_max, &energyChanges_max)
        
        print("E_dec threshold = \(EDecThres)")
        
        if (localEwDec_max > EDecThres) {
            print("merge edge E_dec = \(localEwDec_max)")
            mergeBoundaryEdges(Pair(path_max[0], path_max[1]),
                               Pair(path_max[1], path_max[2]),
                               newVertPos_max.row(0))
            print("edge merged")
            
            computeFeatures() // TODO: only update locally
            return true
        } else {
            print("max E_dec = \(localEwDec_max) < thres \(EDecThres)")
            return false
        }
    }
    
     func splitOrMerge(_ lambda_t: Double,
                             _ EDecThres: Double,
                             _ propogate: Bool,
                             _ splitInterior: Bool,
                             _ isMerge: inout Bool) -> Bool {
        
        assert((!propogate), "propogation is supported separately for split and merge!")
        
        var EwDec_max: Double = .zero
        var path_max: [Int] = []
        var newVertPos_max: Matd = .init()
        isMerge = false
        var energyChanges_split: Pair<Double, Double> = .init(0.0, 0.0)
        var energyChanges_merge: Pair<Double, Double> = .init(0.0, 0.0)
        
        if (splitInterior) {
            querySplit(lambda_t, propogate, splitInterior, &EwDec_max, &path_max, &newVertPos_max, &energyChanges_split)
        } else {
            var EwDec_max_split: Double = .zero
            var EwDec_max_merge: Double = .zero
            var path_max_split: [Int] = []
            var path_max_merge: [Int] = []
            var newVertPos_max_split: Matd = .init()
            var newVertPos_max_merge: Matd = .init()
            querySplit(lambda_t, propogate, splitInterior, &EwDec_max_split, &path_max_split, &newVertPos_max_split, &energyChanges_split)
            queryMerge(lambda_t, propogate, &EwDec_max_merge, &path_max_merge, &newVertPos_max_merge, &energyChanges_merge)
            
            if (EwDec_max_merge > EwDec_max_split) {
                isMerge = true
                EwDec_max = EwDec_max_merge
                path_max = path_max_merge
                newVertPos_max = newVertPos_max_merge
            } else {
                EwDec_max = EwDec_max_split
                path_max = path_max_split
                newVertPos_max = newVertPos_max_split
            }
        }
        
        if (EwDec_max > EDecThres) {
            if (isMerge) {
                print("merge edge E_dec = \(EwDec_max)")
                mergeBoundaryEdges(Pair(path_max[0], path_max[1]),
                                   Pair(path_max[1], path_max[2]),
                                   newVertPos_max.row(0))
                print("edge merged")
                computeFeatures() // TODO: only update locally
            } else {
                if (!splitInterior) {
                    // boundary split
                    print("boundary split E_dec = \(EwDec_max)")
                    splitEdgeOnBoundary(Pair(path_max[0], path_max[1]), newVertPos_max)
                    print("boundary edge splitted")
                    // TODO: process fractail here!
                    updateFeatures()
                } else {
                    // interior split
                    print("Interior split E_dec = \(EwDec_max)")
                    cutPath(path_max, true, 1, newVertPos_max)
                    print("interior edge splitted")
                    fracTail.insert(path_max[0])
                    fracTail.insert(path_max[2])
                    curInteriorFracTails.first = path_max[0]
                    curInteriorFracTails.second = path_max[2]
                    curFracTail = -1
                }
            }
            return true
        } else {
            print("max E_dec = \(EwDec_max) < thres \(EDecThres)")
            return false
        }
    }
    
     func onePointCut(_ vI: Int = 0) {
        assert((vI >= 0) && (vI < V_rest.rows))
        
        var path: [Int] = .init(vNeighbor[vI])
        assert(path.count >= 3)
        path[1] = vI
        path.removeSubrange(3..<path.endIndex)
        
        var makeCoh: Bool = true
        if (!makeCoh) {
            for pI in 0..<(path.count - 1) {
                initSeamLen += (V_rest.row(path[pI]) - V_rest.row(path[pI + 1])).norm()
            }
        }
        
        cutPath(path, makeCoh)
        
        if (makeCoh) {
            initSeams = cohE
        }
    }
    
     func highCurveOnePointCut() {
        var gaussianCurv: [Double] = .init(repeating: 2.0 * .pi, count: V.rows)
        for triI in 0..<F.rows {
            let triVInd: RVec3i = F.row(triI)
            let v: [RVec3d] = [V_rest.row(triVInd[0]),
                               V_rest.row(triVInd[1]),
                               V_rest.row(triVInd[2])]
            for vI in 0..<3 {
                let vI_post = (vI + 1) % 3
                let vI_pre = (vI + 2) % 3
                let e0: RVec3d = v[vI_pre] - v[vI]
                let e1: RVec3d = v[vI_post] - v[vI]
                gaussianCurv[triVInd[vI]] -= acos(max(-1, min(1.0, e0.dot(e1) / e0.norm() / e1.norm())))
            }
        }
        
        for i in 0..<gaussianCurv.count {
            if (gaussianCurv[i] < 0) {
                gaussianCurv[i] = -gaussianCurv[i]
            }
        }
        
        var vI_maxGC: Int = 0
        for vI in 1..<gaussianCurv.count {
            if (gaussianCurv[vI] > gaussianCurv[vI_maxGC]) {
                vI_maxGC = vI
            }
        }
        
        assert(vNeighbor[vI_maxGC].count >= 3)
        var vJ_maxGC: Int = vNeighbor[vI_maxGC].first!
        for vINB in vNeighbor[vI_maxGC] {
            if (gaussianCurv[vINB] > gaussianCurv[vJ_maxGC]) {
                vJ_maxGC = vINB
            }
        }
        
        assert(vNeighbor[vJ_maxGC].count >= 3)
        var vK_maxGC: Int = -1
        var gc_vK: Double = -__DBL_MAX__
        for vJNB in vNeighbor[vJ_maxGC] {
            if (gaussianCurv[vJNB] > gc_vK) && (vJNB != vI_maxGC) {
                vK_maxGC = vJNB
                gc_vK = gaussianCurv[vJNB]
            }
        }
        
        var path: [Int] = [vI_maxGC, vJ_maxGC, vK_maxGC]
        
        var makeCoh: Bool = true
        if (!makeCoh) {
            for pI in 0..<(path.count - 1) {
                initSeamLen += (V_rest.row(path[pI]) - V_rest.row(path[pI + 1])).norm()
            }
        }
        
        cutPath(path, makeCoh)
        
        if (makeCoh) {
            initSeams = cohE
        }
    }
    
     func farthestPointCut() {
        assert(vNeighbor.count == V_rest.rows)
        
        var graph: [[Int:Double]] = .init(repeating: [:], count: vNeighbor.count)
        for vI in 0..<vNeighbor.count {
            for nbI in vNeighbor[vI] {
                if (nbI > vI) {
                    graph[vI][nbI] = (V_rest.row(vI) - V_rest.row(nbI)).norm()
                    graph[nbI][vI] = graph[vI][nbI]
                }
            }
        }
        
        var path: [Int]!
        getFarthestPointPath(graph, getFarthestPoint(graph, 0), &path)
        
        var makeCoh: Bool = true
        if (!makeCoh) {
            for pI in 0..<(path.count - 1) {
                initSeamLen += (V_rest.row(path[pI]) - V_rest.row(path[pI + 1])).norm()
            }
        }
        
        cutPath(path, makeCoh)
        
        if (makeCoh) {
            initSeams = cohE
        }
    }
    
    // Trimesh should be passed by reference
     func geomImgCut(_ data_findExtrema: inout TriMesh) {
        data_findExtrema = self
        // compute UV map for find extremal point (interior)
        let mapType: Int = 0 // 0: SD, 1: harmonic (uniform), 2: harmonic (cotangent), 3: harmonic (MVC)
        if (mapType != 0) {
            var bnd: Vec<Int> = Veci()
            boundary_loop(F, &bnd) // Find the open boundary
            assert(bnd.count != 0)
            
            // TODO: ensure it doesn't have multiple boundaries ? or multi-components?
            
            // Map the boundary to a circle, preserving edge proportions
            var bnd_uv: Matd!
            IglUtils.map_vertices_to_circle(V_rest, bnd, &bnd_uv)
            
            var UV_Tutte: Matd!
            
            switch mapType {
            case 1:
                // harmonic map with uniform weights
                var A: SparseMatrix<Double>!
                var M: SparseMatrix<Double> = .init()
                IglUtils.computeUniformLaplacian(F, &A)
                do {
                    UV_Tutte = try harmonic(A, M, bnd, bnd_uv, 1)
                } catch {
                    print(error.localizedDescription)
                }
                
            case 2:
                // Harmonic parameterization
                do {
                    UV_Tutte = try harmonic(V, F, bnd, bnd_uv, 1)
                } catch {
                    print(error.localizedDescription)
                }
                
            case 3:
                // Shape preserving mesh parameterization
                // (Harmonic map with MVC weights)
                var A: SparseMatrix<Double>!
                IglUtils.computeMVCMtr(V_rest, F, &A)
                IglUtils.fixedBoundaryParam_MVC(A, bnd, bnd_uv: bnd_uv, &UV_Tutte)
                
            default:
                assert(false, "Unknown map specified for finding Geometry Image cuts")
            }
            data_findExtrema.V = UV_Tutte
        }
        // pick the vertex with largest L2 stretch
        var vI_extremal: Int = -1
        var L2stretchPerElem: Vec<Double>!
        var vertScores: Vec<Double>
        data_findExtrema.computeL2StretchPerElem(L2StretchPerElem: &L2stretchPerElem)
        vertScores = .init(V_rest.rows)
        for triI in 0..<F.rows {
            for i in 0..<3 {
                if (vertScores[F[triI, i]] < L2stretchPerElem[triI]) {
                    vertScores[F[triI, i]] = L2stretchPerElem[triI]
                }
            }
        }
        var extremal: Double = 0.0
        for vI in 0..<vertScores.count {
            if (!isBoundaryVert(vI)) {
                if (extremal < vertScores[vI]) {
                    extremal = vertScores[vI]
                    vI_extremal = vI
                }
            }
        }
        assert(vI_extremal >= 0)
        
        // construct mesh graph
        assert(vNeighbor.count == V_rest.rows)
        var graph: [[Int:Double]] = .init(repeating: [:], count: vNeighbor.count)
        for vI in 0..<vNeighbor.count {
            for nbI in vNeighbor[vI] {
                if (nbI > vI) {
                    graph[vI][nbI] = (V_rest.row(vI) - V_rest.row(nbI)).norm()
                    graph[nbI][vI] = graph[vI][nbI]
                }
            }
        }
        
        // find closest point on boundary
        let nV: Int = graph.count
        var dist: [Double]!
        var parent: [Int]!
        dijkstra(graph, vI_extremal, &dist, &parent)
        var minDistToBound: Double = __DBL_MAX__
        var vI_minDistToBound: Int = -1
        for vI in 0..<nV {
            if (isBoundaryVert(vI)) {
                if (dist[vI] < minDistToBound) {
                    minDistToBound = dist[vI]
                    vI_minDistToBound = vI
                }
            }
        }
        assert(vI_minDistToBound >= 0, "No boundary on the mesh!")
        
        // find shortest path to closest point on boundary
        var path: [Int] = []
        while(vI_minDistToBound >= 0) {
            path.append(vI_minDistToBound)
            vI_minDistToBound = parent[vI_minDistToBound]
        }
        path.reverse()
        
        cutPath(path, true)
    }
    
     func cutPath(_ path: [Int],
                        _ makeCoh: Bool = false,
                        _ changePos: Int = 0,
                        _ newVertPos: Mat<Double> = Matd(),
                        _ allowCutThrough: Bool = true) {
        var path = path
        assert(path.count >= 2)
        if (changePos == 1) {
            assert(changePos == 1, "right now only support change 1") // !!! still only allow 1?
            assert(newVertPos.cols == 2)
            assert(changePos * 2 == newVertPos.rows)
        }
        
        for pI in 1..<(path.count - 1) {
            assert(!isBoundaryVert(path[pI]), "Boundary vertices detected on split path, please split multiple times!")
        }
        
        let isFromBound: Bool = isBoundaryVert(path[0])
        let isToBound: Bool = isBoundaryVert(path.last!)
        if (allowCutThrough && (isFromBound || isToBound)) {
            var cutThrough: Bool = false
            if (isFromBound && isToBound) {
                cutThrough = true
            } else if (isToBound) {
                path.reverse()
                // always start cut from boundary
            }
            
            for vI in 0..<(path.count - 1) {
                let vInd_s: Int = path[vI]
                let vInd_e: Int = path[vI + 1]
                assert(edge2Tri[Pair(vInd_s, vInd_e)] != nil)
                assert(edge2Tri[Pair(vInd_e, vInd_s)] != nil)
                var newVertPos: Matd
                if (cutThrough) {
                    newVertPos = .init(4, 2)
                    newVertPos <<== [V.row(vInd_s), V.row(vInd_s), V.row(vInd_e), V.row(vInd_e)]
                } else {
                    newVertPos = .init(2, 2)
                    newVertPos <<== [V.row(vInd_s), V.row(vInd_s)]
                }
                splitEdgeOnBoundary(Pair(vInd_s, vInd_e), newVertPos, true)
                updateFeatures()
            }
        } else {
            // path is interior
            assert(path.count >= 3)
            
            var tri_left: [Int] = []
            let vI = path[1]
            var vI_new = path[0]
            while (true) {
                let tri = edge2Tri[Pair(vI, vI_new)]
                assert(tri != nil)
                tri_left.append(tri!)
                let triVInd: RVec3i = F.row(tri!)
                for i in 0..<3 {
                    if (triVInd[i] != vI) && (triVInd[i] != vI_new) {
                        vI_new = triVInd[i]
                        break
                    }
                }
                
                if (vI_new == path[2]) {
                    break
                }
                if (vI_new == path[0]) {
                    assert(false, "not a valid path!")
                }
            }
            
            let nV: Int = V_rest.rows
            V_rest.conservativeResize(nV + 1, 3)
            V_rest.row(nV) <<== V_rest.row(path[1])
            vertWeight.conservativeResize(nV + 1)
            vertWeight[nV] = vertWeight[path[1]]
            V.conservativeResize(nV + 1, 2)
            if (changePos != 0) {
                V.row(nV) <<== newVertPos.block(0, 0, 1, 2)
                V.row(path[1]) <<== newVertPos.block(1, 0, 1, 2)
            } else {
                V.row(nV) <<== V.row(path[1])
            }
            for triI in tri_left {
                for vI in 0..<3 {
                    if (F[triI, vI] == path[1]) {
                        F[triI, vI] = nV
                        break
                    }
                }
            }
            if (makeCoh) {
                let nCoh: Int = cohE.rows
                cohE.conservativeResize(nCoh + 2, 4)
                cohE.row(nCoh) <<== [nV, path[0], path[1], path[0]]
                cohE.row(nCoh + 1) <<== [path[2], nV, path[2], path[1]]
            }
            
            computeFeatures() // TODO: only update locally
            
            for vI in 2..<(path.count - 1) {
                let vInd_s: Int = path[vI]
                let vInd_e = path[vI + 1]
                assert(edge2Tri[Pair(vInd_s, vInd_e)] != nil)
                assert(edge2Tri[Pair(vInd_e, vInd_s)] != nil)
                var newVertPos: Matd = .init(2, 2)
                newVertPos <<== [V.row(vInd_s), V.row(vInd_s)]
                splitEdgeOnBoundary(Pair(vInd_s, vInd_e), newVertPos, true, allowCutThrough) // !!! make coh?
                updateFeatures()
            }
        }
    }
    
     func computeSeamScore(seamScore: inout Vec<Double>!) {
        seamScore = .init(cohE.rows)
        for cohI in 0..<cohE.rows {
            if (boundaryEdge[cohI] != 0) {
                seamScore[cohI] = -1.0
            } else {
                seamScore[cohI] = max((V.row(cohE[cohI, 0]) - V.row(cohE[cohI, 2])).norm(),
                                      (V.row(cohE[cohI, 1]) - V.row(cohE[cohI, 3])).norm()) / avgEdgeLen
            }
        }
    }
    
     func computeBoundaryLen(boundaryLen: inout Double) {
        boundaryLen = 0.0
        for (edge, _) in edge2Tri {
            if edge2Tri[Pair(edge.second, edge.first)] == nil {
                boundaryLen += (V_rest.row(edge.second) - V_rest.row(edge.first)).norm()
            }
        }
    }
    
     func computeSeamSparsity(_ sparsity: inout Double, _ triSoup: Bool = false) {
        let thres: Double = 1.0e0-2
        sparsity = 0.0
        for cohI in 0..<cohE.rows {
            if (boundaryEdge[cohI] == 0) {
                if (!triSoup ||
                    ((V.row(cohE[cohI, 0]) - V.row(cohE[cohI, 2])).norm() / avgEdgeLen > thres) ||
                    ((V.row(cohE[cohI, 1]) - V.row(cohE[cohI, 3])).norm() / avgEdgeLen > thres)) {
                    sparsity += edgeLen[cohI]
                }
            }
        }
        sparsity += initSeamLen
    }
    
     func computeStandardStretch(_ stretch_l2: inout Double,
                                       _ stretch_inf: inout Double,
                                       _ stretch_shear: inout Double,
                                       _ compress_inf: inout Double) {
        stretch_l2 = 0.0
        stretch_inf = -__DBL_MAX__
        stretch_shear = 0.0
        compress_inf = __DBL_MAX__
        for triI in 0..<F.rows {
            let triVInd: Vec3i = F.row(triI)
            let x_3D: [Vec3d] = [
                V_rest.row(triVInd[0]),
                V_rest.row(triVInd[1]),
                V_rest.row(triVInd[2])
            ]
            let uv: [Vec2d] = [
                V.row(triVInd[0]),
                V.row(triVInd[1]),
                V.row(triVInd[2])
            ]
            var dg: Mat2d = .init()
            IglUtils.computeDeformationGradient(x_3D, uv, &dg)
            
            let a: Double = dg.block(0, 0, 2, 1).squaredNorm()
            let b: Double = Vec2d(dg.block(0, 0, 2, 1)).dot(Vec2d(dg.block(0, 1, 2, 1)))
            let c: Double = dg.block(0, 1, 2, 1).squaredNorm()
            let t0: Double = a + c
            let t1: Double = sqrt((a - c) * (a - c) + 4.0 * b * b)
            let tau: Double = sqrt((t0 + t1) / 2.0)
            let gamma: Double = sqrt((t0 - t1) / 2.0)
            
            stretch_l2 += (t0 / 2.0 * triArea[triI])
            
            if (tau > stretch_inf) {
                stretch_inf = tau
            }
            
            stretch_shear += b * b / a / c * triArea[triI]
            
            if (gamma < compress_inf) {
                compress_inf = gamma
            }
        }
        stretch_l2 /= surfaceArea
        stretch_l2 = sqrt(stretch_l2)
        stretch_shear /= surfaceArea
        stretch_shear = sqrt(stretch_shear)
        
        var surfaceArea_UV: Double = 0.0
        for triI in 0..<F.rows {
            let triVInd: Vec3i = F.row(triI)
            
            let U1: Vec2d = V.row(triVInd[0])
            let U2: Vec2d = V.row(triVInd[1])
            let U3: Vec2d = V.row(triVInd[2])
            
            let U2m1: Vec2d = U2 - U1
            let U3m1: Vec2d = U3 - U1
            
            surfaceArea_UV += 0.5 * (U2m1[0] * U3m1[1] - U2m1[1] * U3m1[0])
        }
        
        // area scaling:
        let scaleFactor = sqrt(surfaceArea_UV / surfaceArea)
        stretch_l2 *= scaleFactor
        stretch_inf *= scaleFactor
        compress_inf *= scaleFactor // not meaningful now...
        // stretch_shear won't be affected by area scaling
    }
    
     func computeL2StretchPerElem(L2StretchPerElem: inout Vec<Double>!) {
        L2StretchPerElem = .init(F.rows)
        for triI in 0..<F.rows {
            let triVInd: Vec3i = F.row(triI)
            let x_3D: [Vec3d] = [V_rest.row(triVInd[0]),
                                 V_rest.row(triVInd[1]),
                                 V_rest.row(triVInd[2])]
            let uv: [Vec2d] = [V.row(triVInd[0]),
                               V.row(triVInd[1]),
                               V.row(triVInd[2])]
            var dg: Mat2d!
            IglUtils.computeDeformationGradient(x_3D, uv, &dg)
            
            let a: Double = dg.block(0, 0, 2, 1).squaredNorm()
            let c: Double = dg.block(0, 1, 2, 1).squaredNorm()
            let t0: Double = a + c
            
            L2StretchPerElem[triI] = sqrt(t0 / 2.0)
        }
    }
    
    // Note: Data is std::ofstream in c++
     func outputStandardStretch(file: inout Data) {
        var stretch_l2: Double = 0.0
        var stretch_inf: Double = 0.0
        var stretch_shear: Double = 0.0
        var compress_inf: Double = 0.0
        computeStandardStretch(&stretch_l2, &stretch_inf, &stretch_shear, &compress_inf)
        let string = String(stretch_l2) + " " + String(stretch_inf) + " " + String(stretch_shear) + " " + String(compress_inf)
        file.append(Data(string.utf8))
    }
    
     func computeAbsGaussianCurv(_ absGaussianCurv: inout Double) {
        var weights: [Double] = .init(repeating: 0.0, count: F.rows)
        var gaussianCurv: [Double] = .init(repeating: 2.0 * .pi, count: V.rows)
        for triI in 0..<F.rows {
            let triVInd: RVec3i = F.row(triI)
            var v: [RVec3d] = [
                V_rest.row(triVInd[0]),
                V_rest.row(triVInd[1]),
                V_rest.row(triVInd[2])
            ]
            for vI in 0..<3 {
                let vI_post = (vI + 1) % 3
                let vI_pre = (vI + 2) % 3
                let e0: RVec3d = v[vI_pre] - v[vI]
                let e1: RVec3d = v[vI_post] - v[vI]
                gaussianCurv[triVInd[vI]] -= acos(max(-1.0, min(1.0, e0.dot(e1) / e0.norm() / e1.norm())))
                weights[triVInd[vI]] += triArea[triI]
            }
        }
        
        absGaussianCurv = 0.0
        for vI in 0..<V.rows {
            if (!isBoundaryVert(vI)) {
                absGaussianCurv += abs(gaussianCurv[vI]) * weights[vI]
            }
        }
        absGaussianCurv /= surfaceArea * 3.0
    }
    
     func initRigidUV() {
        V.resize(V_rest.rows, 2)
        for triI in 0..<F.rows {
            let triVInd: Vec3i = F.row(triI)
            
            let x_3D: [Vec3d] = [
                V_rest.row(triVInd[0]),
                V_rest.row(triVInd[1]),
                V_rest.row(triVInd[2])
            ]
            var x: [Vec2d]!
            IglUtils.mapTriangleTo2D(x_3D, &x)
            
            V.row(triVInd[0]) <<== x[0]
            V.row(triVInd[1]) <<== x[1]
            V.row(triVInd[2]) <<== x[2]
        }
    }
    
     func checkInversion(_ triI: Int, _ mute: Bool) -> Bool {
        assert(triI < F.rows)
        
        let eps: Double = 0.0
        
        let triVInd: Vec3i = F.row(triI)
        
        let e_u: [Vec2d] = [
            V.row(triVInd[1]) - V.row(triVInd[0]),
            V.row(triVInd[2]) - V.row(triVInd[0])
        ]
        
        let dbArea: Double = e_u[0][0] * e_u[1][1] - e_u[0][1] * e_u[1][0]
        if (dbArea < eps) {
            if (!mute) {
                print("***Element inversion detected: \(dbArea) < \(eps)")
                print("mesh triangle count: \(F.rows)")
            }
            return false
        } else {
            return true
        }
    }
    
     func checkInversion(_ mute: Bool = false,
                        _ triangles: [Int] = []) -> Bool {
        if (triangles.isEmpty) {
            for triI in 0..<F.rows {
                if (!checkInversion(triI, mute)) {
                    return false
                }
            }
        } else {
            for triI in triangles {
                if (!checkInversion(triI, mute)) {
                    return false
                }
            }
        }
        
        return true
    }
    
    // Helper functions
     func computeLabpaciamMtr() {
        let L: SparseMatrix<Double> = cotmatrix(V_rest, F)
        LaplacianMtr.resize(L.rows, L.cols)
        LaplacianMtr.setZero()
        LaplacianMtr.reserve(L.nonZeros)
        
        var triplets: [Triplet<Double>] = []
        triplets.reserveCapacity(L.nonZeros)
        
        for k in 0..<L.outerSize {
            for it in L.innerIterator(k) {
                if (!fixedVert.contains(it.row) && !fixedVert.contains(it.col)) {
                    triplets.append(.init(i: it.row, j: it.col, value: -it.value))
                }
            }
        }
        for fixedVI in fixedVert {
            triplets.append(.init(i: fixedVI, j: fixedVI, value: 1.0))
        }
        LaplacianMtr.setFromTriplets(triplets)
    }
    
     func findBoundaryEdge(_ vI: Int,
                          _ startEdge: Pair<Int, Int>,
                          _ boundaryEdge: inout Pair<Int, Int>) -> Bool {
        guard var triI = edge2Tri[startEdge] else { fatalError() }
        
        let proceed: Bool = startEdge.first != vI
        let vI_neighbor: Int = (proceed ? startEdge.first : startEdge.second)
        var vI_new: Int = vI_neighbor
        while (true) {
            let triVInd: RVec3i = F.row(triI)
            for i in 0..<3 {
                if ((triVInd[i] != vI) && (triVInd[i] != vI_new)) {
                    vI_new = triVInd[i]
                    break
                }
            }
            
            if (vI_new == vI_neighbor) {
                return false
            }
            
            let triI = edge2Tri[proceed ? Pair(vI_new, vI) : Pair(vI, vI_new)]
            if triI == nil {
                boundaryEdge.first = (proceed ? vI : vI_new)
                boundaryEdge.second = (proceed ? vI_new : vI)
                
                return true
            }
        }
    }
    
     func insideTri(_ triI: Int, _ pos: RVec2d) -> Bool {
        let triVInd: RVec3i = F.row(triI)
        let e01: RVec2d = V.row(triVInd[1]) - V.row(triVInd[0])
        let e02: RVec2d = V.row(triVInd[2]) - V.row(triVInd[0])
        let e0p: RVec2d = pos - V.row(triVInd[0])
        
        // represent ep0 using e01 and e02
        var A: Mat2d = .init([e01[0], e02[0], e01[1], e02[1]], [2, 2])
        
        let bValues: [Double] = e0p.values
        fatalError("To be implemented")
    }
    
     func insideUVRegion(_ triangles: [Int], _ pos: RVec2d) -> Bool {
        for triI in triangles {
            if (insideTri(triI, pos)) {
                return true
            }
        }
        return false
    }
    
    // toBound = false indicate counter-clockwise
    @discardableResult
     func isBoundaryVert(_ vI: Int,
                               _ vI_neighbor: Int,
                               _ tri_toSep: inout [Int],
                               _ boundaryEdge: inout Pair<Int, Int>,
                               _ toBound: Bool = true) -> Bool {
        tri_toSep = []
        var triI: Int! = edge2Tri[toBound ? Pair(vI_neighbor, vI) : Pair(vI, vI_neighbor)]
        
        if triI == nil {
            boundaryEdge.first = (toBound ? vI : vI_neighbor)
            boundaryEdge.second = (toBound ? vI_neighbor : vI)
            return true
        }
        
        var vI_new = vI_neighbor
        
        repeat {
            tri_toSep.append(triI)
            let triVInd: RVec3i = F.row(triI)
            for i in 0..<3 {
                if (triVInd[i] != vI) && (triVInd[i] != vI_new) {
                    vI_new = triVInd[i]
                    break
                }
            }
            
            if (vI_new == vI_neighbor) {
                return false
            }
            
            triI = edge2Tri[toBound ? Pair(vI_new, vI) : Pair(vI, vI_new)]
            
            if (triI == nil) {
                boundaryEdge.first = (toBound ? vI : vI_new)
                boundaryEdge.second = (toBound ? vI_new : vI)
                return true
            }
        } while (true)
    }
    
     func isBoundaryVert(_ vI: Int) -> Bool {
        assert(vNeighbor.count == V.rows)
        assert(vI < vNeighbor.count)
        
        for vI_neighbor in vNeighbor[vI] {
            if (edge2Tri[Pair(vI, vI_neighbor)] == nil || edge2Tri[Pair(vI_neighbor, vI)] == nil) {
                return true
            }
        }
        
        return false
    }
    
     func compute2DInwardNormal(_ vI: Int,
                                      _ normal: inout RVec2d?) {
        var incTris: [[Int]] = .init(repeating: [], count: 2)
        var boundaryEdge: [Pair<Int, Int>] = .init(repeating: Pair(0, 0), count: 2)
        if (!isBoundaryVert(vI, vNeighbor[vI].first!, &incTris[0], &boundaryEdge[0], false)) {
            return
        }
        isBoundaryVert(vI, vNeighbor[vI].first!, &incTris[1], &boundaryEdge[1], true)
        assert(!(incTris[0].isEmpty && incTris[1].isEmpty))
        
        let boundaryEdgeDir: [RVec2d] = [
            (V.row(boundaryEdge[0].first) - V.row(boundaryEdge[0].second)).normalized(),
            (V.row(boundaryEdge[1].second) - V.row(boundaryEdge[1].first)).normalized()
        ]
        normal = (boundaryEdgeDir[0] + boundaryEdgeDir[1]).normalized()
        if (boundaryEdgeDir[1][0] * normal![1] - boundaryEdgeDir[1][1] * normal![0] < 0.0) {
            normal! *= -1
        }
    }
    
     func splitEdgeOnBoundary(_ edge: Pair<Int, Int>,
                                    _ newVertPos: Mat<Double>,
                                    _ changeVertPos: Bool = true,
                                    _ allowCutThrough: Bool = true) {
        
        assert(vNeighbor.count == V.rows)
        guard let _ = edge2Tri[edge] else { fatalError() }
        guard let _ = edge2Tri[Pair(edge.second, edge.first)] else { fatalError() }
        
        var duplicateBoth: Bool = false
        var vI_boundary: Int = edge.first
        var vI_interior: Int = edge.second
        if (isBoundaryVert(edge.first)) {
            if (allowCutThrough && isBoundaryVert(edge.second)) {
                if (changeVertPos) {
                    assert(newVertPos.rows == 4)
                }
                duplicateBoth = true
            }
        } else {
            assert(isBoundaryVert(edge.second), "Input edge must attach mesh boundary!")
            
            vI_boundary = edge.second
            vI_interior = edge.first
        }
        
        fracTail.remove(vI_boundary)
        if (!duplicateBoth) {
            fracTail.insert(vI_interior)
            curFracTail = vI_interior
        } else {
            curFracTail = -1
        }
        curInteriorFracTails.second = -1
        curInteriorFracTails.first = -1
        
        // duplicate vI_boundary
        var tri_toSep: [[Int]] = .init(repeating: [], count: 2)
        var boundaryEdge: [Pair<Int, Int>] = .init(repeating: Pair(0, 0), count: 2)
        for toBound in 0..<2 {
            isBoundaryVert(vI_boundary, vI_interior, &tri_toSep[1], &boundaryEdge[1], toBound != 0)
            assert(!tri_toSep[1].isEmpty)
        }
        if (duplicateBoth) {
            isBoundaryVert(vI_interior, vI_boundary, &tri_toSep[0], &boundaryEdge[0], true)
            assert(!tri_toSep[0].isEmpty)
        }
        
        let nV = V_rest.rows
        V_rest.conservativeResize(nV + 1, 3)
        V_rest.row(nV) <<== V_rest.row(vI_boundary)
        vertWeight.conservativeResize(nV + 1)
        vertWeight[nV] = vertWeight[vI_boundary]
        V.conservativeResize(nV + 1, 2)
        if (changeVertPos) {
            V.row(nV) <<== newVertPos.block(1, 0, 1, 2)
            V.row(vI_boundary) <<== newVertPos.block(0, 0, 1, 2)
        } else {
            V.row(nV) <<== V.row(vI_boundary)
        }
        
        for triI in tri_toSep[1] {
            for vI in 0..<3 {
                if (F[triI, vI] == vI_boundary) {
                    // update triangle vertInd, edge2Tri and vNeighbor
                    let vI_post = F[triI, (vI + 1) % 3]
                    let vI_pre = F[triI, (vI + 2) % 3]
                    
                    F[triI, vI] = nV
                    
                    edge2Tri[Pair(vI_boundary, vI_post)] = nil
                    edge2Tri[Pair(nV, vI_post)] = triI
                    edge2Tri[Pair(vI_pre, vI_boundary)] = nil
                    edge2Tri[Pair(vI_pre, nV)] = triI
                    
                    vNeighbor[vI_pre].remove(vI_boundary)
                    vNeighbor[vI_pre].insert(nV)
                    vNeighbor[vI_post].remove(vI_boundary)
                    vNeighbor[vI_post].insert(nV)
                    vNeighbor[vI_boundary].remove(vI_pre)
                    vNeighbor[vI_boundary].remove(vI_post)
                    vNeighbor.conservativeResize(to: nV + 1)
                    vNeighbor[nV].insert(vI_pre)
                    vNeighbor[nV].insert(vI_post)
                    
                    break
                }
            }
        }
        
        vNeighbor[vI_boundary].insert(vI_interior)
        vNeighbor[vI_interior].insert(vI_boundary)
        
        // add cohesive edge pair and update cohEIndex
        let nCE: Int = cohE.rows
        cohE.conservativeResize(nCE + 1, 4)
        cohE.row(nCE) <<== [vI_interior, nV, vI_interior, vI_boundary]
        cohEIndex[Pair(vI_interior, nV)] = nCE
        cohEIndex[Pair(vI_boundary, vI_interior)] = -nCE - 1
        if let CEI = cohEIndex[boundaryEdge[1]] {
            if CEI >= 0 {
                cohE[CEI, 0] = nV
            } else {
                cohE[-CEI - 1, 3] = nV
            }
            cohEIndex[Pair(nV, boundaryEdge[1].second)] = CEI
            cohEIndex.removeValue(forKey: boundaryEdge[1])
        }
        
        if (duplicateBoth) {
            let nV: Int = V_rest.rows
            V_rest.conservativeResize(nV + 1, 3)
            V_rest.row(nV) <<== V_rest.row(vI_interior)
            vertWeight.conservativeResize(nV + 1)
            vertWeight[nV] = vertWeight[vI_interior]
            V.conservativeResize(nV + 1, 2)
            if changeVertPos {
                V.row(nV) <<== newVertPos.block(2, 0, 1, 2)
                V.row(vI_interior) <<== newVertPos.block(3, 0, 1, 2)
            } else {
                V.row(nV) <<== V.row(vI_interior)
            }
            
            for triI in tri_toSep[0] {
                for vI in 0..<3 {
                    if (F[triI, vI] == vI_interior) {
                        // update triangle vertInd, edge2Tri and Vneighbor
                        let vI_post = F[triI, (vI + 1) % 3]
                        let vI_pre = F[triI, (vI + 2) % 3]
                        
                        F[triI, vI] = nV
                        
                        edge2Tri[Pair(vI_interior, vI_post)] = nil
                        edge2Tri[Pair(nV, vI_post)] = triI
                        edge2Tri[Pair(vI_pre, vI_interior)] = nil
                        edge2Tri[Pair(vI_pre, nV)] = triI
                        
                        vNeighbor[vI_pre].remove(vI_interior)
                        vNeighbor[vI_pre].insert(nV)
                        vNeighbor[vI_post].remove(vI_interior)
                        vNeighbor[vI_post].insert(nV)
                        vNeighbor[vI_interior].remove(vI_pre)
                        vNeighbor[vI_interior].remove(vI_post)
                        vNeighbor.conservativeResize(to: nV + 1)
                        vNeighbor[nV].insert(vI_pre)
                        vNeighbor[nV].insert(vI_post)
                        
                        break
                    }
                }
            }
            
            // update cohesive edge pair and update cohEIndex
            cohE[nCE, 2] = nV
            cohEIndex.removeValue(forKey: Pair(vI_boundary, vI_interior))
            cohEIndex[Pair(vI_boundary, nV)] = -nCE - 1
            if let CEI = cohEIndex[boundaryEdge[0]] {
                if CEI >= 0 {
                    cohE[CEI, 0] = nV
                } else {
                    cohE[-CEI - 1, 3] = nV
                }
                cohEIndex[Pair(nV, boundaryEdge[0].second)] = CEI
                cohEIndex.removeValue(forKey: boundaryEdge[0])
            }
        }
     }
    
     func mergeBoundaryEdges(_ edge0: Pair<Int, Int>,
                            _ edge1: Pair<Int, Int>,
                            _ mergedPos: RVecd) {
        
        assert(edge0.second == edge1.first)
        assert(edge2Tri[Pair(edge0.second, edge0.first)] == nil)
        assert(edge2Tri[Pair(edge1.second, edge1.first)] == nil)
        
        fracTail.remove(edge0.second)
        fracTail.insert(edge0.first)
        curFracTail = edge0.first
        
        V.row(edge0.first) <<== mergedPos
        let vBackI: Int = V.rows - 1
        if (edge1.second < vBackI) {
            V_rest.row(edge1.second) <<== V_rest.row(vBackI)
            vertWeight[edge1.second] = vertWeight[vBackI]
            V.row(edge1.second) <<== V.row(vBackI)
            
            if fracTail.contains(vBackI) {
                fracTail.remove(vBackI)
                fracTail.insert(edge1.second)
            }
        } else {
            assert(edge1.second == vBackI)
        }
        V_rest.conservativeResize(vBackI, 3)
        vertWeight.conservativeResize(vBackI)
        V.conservativeResize(vBackI, 2)
        
        for triI in 0..<F.rows {
            for vI in 0..<3 {
                if (F[triI, vI] == edge1.second) {
                    F[triI, vI] = edge0.first
                    break
                }
            }
        }
        
        if (edge1.second < vBackI) {
            for triI in 0..<F.rows {
                for vI in 0..<3 {
                    if (F[triI, vI] == vBackI) {
                        F[triI, vI] = edge1.second
                    }
                }
            }
        }
        
        guard let cohInd = cohEIndex[edge0] else { fatalError() }
        
        let cohEBackI: Int = cohE.rows - 1
        if (cohInd >= 0) {
            if (cohInd < cohEBackI) {
                cohE.row(cohInd) <<== cohE.row(cohEBackI)
            } else {
                assert(cohInd == cohEBackI)
            }
        } else {
            if (-cohInd - 1 < cohEBackI) {
                cohE.row(-cohInd - 1) <<== cohE.row(cohEBackI)
            } else {
                assert(-cohInd - 1 == cohEBackI)
            }
        }
        cohE.conservativeResize(cohEBackI, 4)
        
        for cohI in 0..<cohE.rows {
            for pI in 0..<4 {
                if (cohE[cohI, pI] == edge1.second) {
                    cohE[cohI, pI] = edge0.first
                }
            }
        }
        if (edge1.second < vBackI) {
            for cohI in 0..<cohE.rows {
                for pI in 0..<4 {
                    if (cohE[cohI, pI] == vBackI) {
                        cohE[cohI, pI] = edge1.second
                    }
                }
            }
        }
        
        // closeup just interior splitted diamond
        // TODO: do it faster by knowing the edge in advance and locate using cohIndex
        for cohI in 0..<cohE.rows {
            if ((cohE[cohI, 0] == cohE[cohI, 2]) && (cohE[cohI, 1] == cohE[cohI, 3])) {
                fracTail.remove(cohE[cohI, 0])
                fracTail.remove(cohE[cohI, 1])
                curFracTail = -1
                
                if (cohI < cohE.rows - 1) {
                    cohE.row(cohI) <<== cohE.row(cohE.rows - 1)
                }
                cohE.conservativeResize(cohE.rows - 1, 4)
                break
            }
        }
        // TODO: locally update edge2Tri, vNeighbor, cohEIndex
    }
    
    // query vertex candidate for either split or merge
     func computeLocalLDec(_ vI: Int,
                                 _ lambda_t: Double,
                                 _ path_max: inout [Int],
                                 _ newVertPos_max: inout Mat<Double>,
                                 _ energyChanges_max: inout Pair<Double, Double>,
                                 _ incTris: [Int] = [],
                                 _ initMergedPos: RVec2d = RVec2d()) -> Double {
        
        if (!path_max.isEmpty) {
            // merge query
            assert(path_max.count >= 3)
            for pI in path_max {
                assert(isBoundaryVert(pI))
            }
            
            if (path_max.count == 3) {
                // zipper merge
                let numerator: Double = (V_rest.row(path_max[0]) - V_rest.row(path_max[1])).norm()
                let denominator: Double = virtualRadius * (vertWeight[path_max[0]] + vertWeight[path_max[1]])
                var seDec: Double = numerator / denominator / 2.0
                // closing up splitted diamond
                var closeup: Bool = false
                for nbVI in vNeighbor[path_max[0]] {
                    if (nbVI != path_max[1]) {
                        if (isBoundaryVert(nbVI)) {
                            if (vNeighbor[path_max[2]].contains(nbVI)) {
                                let numerator: Double = (V_rest.row(path_max[0]) - V_rest.row(nbVI)).norm()
                                let denominator: Double = virtualRadius * (vertWeight[path_max[0]] + vertWeight[nbVI])
                                seDec += numerator / denominator / 2.0
                                closeup = true
                                break
                            }
                        }
                    }
                }
                
                assert(incTris.count >= 2)
                var freeVert = Set<Int>()
                freeVert.insert(path_max[0])
                freeVert.insert(path_max[2])
                var mergeVert = [Int:Int]()
                mergeVert[path_max[0]] = path_max[2]
                mergeVert[path_max[2]] = path_max[0]
                var newVertPos = [Int:RVec2d]()
                let SDInc: Double = -computeLocalEdDec_merge(path_max, incTris, freeVert, &newVertPos, mergeVert, initMergedPos, closeup)
                energyChanges_max.first = SDInc
                energyChanges_max.second = -seDec
                if (SDInc == __DBL_MAX__) {
                    return -__DBL_MAX__
                } else {
                    let value = newVertPos[path_max[0]]
                    assert(value != nil)
                    newVertPos_max.resize(1, 2)
                    newVertPos_max.row(0) <<== value!
                    
                    return lambda_t * seDec - (1.0 - lambda_t) * SDInc
                }
            } else {
                assert(false, "currently not considering \"interior\" merge!")
            }
        }
        
        // split:
        var umbrella: [Int] = []
        var boundaryEdge: Pair<Int, Int> = Pair(0, 0)
        if (isBoundaryVert(vI, vNeighbor[vI].first!, &umbrella, &boundaryEdge, false)) {
            // boundary split
            var maxEwDec = -__DBL_MAX__
            energyChanges_max.first = __DBL_MAX__
            energyChanges_max.second = __DBL_MAX__
            if path_max.count > 2 {
                path_max.removeSubrange(2..<path_max.endIndex)
            }
            if path_max.isEmpty {
                path_max = .init(repeating: 0, count: 2)
            }
            for nbVI in vNeighbor[vI] {
                let edge = Pair(vI, nbVI)
                if (edge2Tri[edge] != nil && edge2Tri[Pair(nbVI, vI)] != nil) {
                    // interior edge
                    var newVertPosI = Matd()
                    let SDDec: Double = queryLocalEdDec_bSplit(edge, &newVertPosI)
                    
                    let numerator: Double = (V_rest.row(vI) - V_rest.row(nbVI)).norm()
                    let denominator: Double = virtualRadius * (vertWeight[vI] + vertWeight[nbVI]) / 2.0
                    let seInc: Double = numerator / denominator
                    let curEWDec = (1.0 - lambda_t) * SDDec - lambda_t * seInc
                    if (curEWDec > maxEwDec) {
                        maxEwDec = curEWDec
                        path_max[0] = vI
                        path_max[1] = nbVI
                        newVertPos_max = newVertPosI
                        energyChanges_max.first = -SDDec
                        energyChanges_max.second = seInc
                    }
                }
            }
            return maxEwDec
        } else {
            // interior split
            for nbVI in vNeighbor[vI] {
                if (isBoundaryVert(nbVI)) {
                    energyChanges_max.first = __DBL_MAX__
                    energyChanges_max.second = __DBL_MAX__
                    assert(false, "should have prevented this case outside")
                    return -__DBL_MAX__ // don't split vertices connected to boundary here
                }
            }
            
            if (umbrella.count > 10) {
                print("large degree vert, \(umbrella.count) incident tris")
            }
            
            //if path_max.count > 3 { path_max.removeSubrange(3..<path_max.endIndex)}
            path_max.conservativeResize(to: 3)
            var EwDec_max = -__DBL_MAX__
            energyChanges_max.first = __DBL_MAX__
            energyChanges_max.second = __DBL_MAX__
            var freeVert = Set<Int>()
            freeVert.insert(vI)
            var newVertPosMap = [Int:RVec2d]()
            var path = [Int].init(repeating: 0, count: 3)
            path[1] = vI
            for startI in 0..<(umbrella.count - 1) {
                for i in 0..<3 {
                    if (F[umbrella[startI], i] == vI) {
                        path[0] = F[umbrella[startI], (i + 1) % 3]
                        break
                    }
                }
                
                for endI in (startI + 1)..<umbrella.count {
                    for i in 0..<3 {
                        if (F[umbrella[endI], i] == vI) {
                            path[2] = F[umbrella[endI], (i + 1) % 3]
                            break
                        }
                    }
                    
                    var SDDec: Double = 0.0
                    var newVertPos: Matd = .init()
                    SDDec += computeLocalEdDec_inSplit(umbrella, freeVert, path, &newVertPos)
                    // TODO: share local mesh before split, also for boundary splits
                    
                    let expression1: Double = (V_rest.row(path[0]) - V_rest.row(path[1])).norm() * (vertWeight[path[0]] + vertWeight[path[1]])
                    let expression2: Double = (V_rest.row(path[1]) - V_rest.row(path[2])).norm() * (vertWeight[path[1]] + vertWeight[path[2]])
                    let seInc: Double = (expression1 + expression2) / virtualRadius / 2.0
                    let EwDec = (1.0 - lambda_t) * SDDec - lambda_t * seInc
                    if (EwDec > EwDec_max) {
                        EwDec_max = EwDec
                        newVertPos_max = newVertPos
                        path_max = path
                        energyChanges_max.first = -SDDec
                        energyChanges_max.second = seInc
                    }
                }
            }
            return EwDec_max
        }
    }
    
    // query interior incident edge of a boundary vertex candidate
     func queryLocalEdDec_bSplit(_ edge: Pair<Int, Int>,
                                       _ newVertPos: inout Mat<Double>) -> Double {
        
        assert(vNeighbor.count == V.rows)
        guard let _ = edge2Tri[edge] else { fatalError() }
        guard let _ = edge2Tri[Pair(edge.second, edge.first)] else { fatalError() }
        
        var vI_boundary: Int = edge.first
        var vI_interior: Int = edge.second
        var cutThrough: Bool = false
        if (isBoundaryVert(edge.first)) {
            if (isBoundaryVert(edge.second)) {
                cutThrough = true
            }
        } else {
            assert(isBoundaryVert(edge.second), "Input edge must attach boundary")
            
            vI_boundary = edge.second
            vI_interior = edge.first
        }
        
        if (cutThrough) {
            newVertPos.resize(4, 2)
        } else {
            newVertPos.resize(2, 2)
        }
        
        var freeVertGID = Set<Int>()
        freeVertGID.insert(vI_boundary)
        if (cutThrough) {
            freeVertGID.insert(vI_interior)
        }
        
        var tri_toSep = [Int]()
        var tri_toSep1 = [Int]()
        var boundaryEdge: Pair<Int, Int> = Pair(0, 0)
        isBoundaryVert(vI_boundary, vI_interior, &tri_toSep, &boundaryEdge, false)
        assert(!tri_toSep.isEmpty)
        isBoundaryVert(vI_boundary, vI_interior, &tri_toSep1, &boundaryEdge, true)
        assert(!tri_toSep1.isEmpty)
        tri_toSep.append(contentsOf: tri_toSep1)
        
        if (cutThrough) {
            for clockwide in 0..<2 {
                var tri_interior: [Int] = []
                var boundaryEdge_interior: Pair<Int, Int> = Pair(0, 0)
                isBoundaryVert(vI_interior, vI_boundary, &tri_interior, &boundaryEdge_interior, clockwide == 0 ? false : true)
                for triI in tri_interior {
                    var newTri: Bool = true
                    for triI_b in tri_toSep {
                        if (triI_b == triI) {
                            newTri = false
                            break
                        }
                    }
                    if (newTri) {
                        tri_toSep.append(triI)
                    }
                }
            }
        }
        
        var splitPath: [Int] = .init(repeating: 0, count: 2)
        splitPath[0] = vI_boundary
        splitPath[1] = vI_interior
        
        return computeLocalEdDec_bSplit(tri_toSep, freeVertGID, splitPath, &newVertPos)
    }
    
    // boundary split
     func computeLocalEdDec_bSplit(_ triangles: [Int],
                                         _ freeVert: Set<Int>,
                                         _ splitPath: [Int],
                                         _ newVertPos: inout Mat<Double>,
                                         _ maxIter: Int = 100) -> Double {
        
        assert(!triangles.isEmpty && !freeVert.isEmpty)
        
        // construct local mesh
        var localF = Mati(triangles.count, 3)
        var localV_rest = Matd()
        var localV = Matd()
        var fixedVert = Set<Int>()
        var globalVI2Local = [Int:Int]()
        var localTriI: Int = 0
        for triI in triangles {
            for vI in 0..<3 {
                let globalVI = F[triI, vI]
                if let local = globalVI2Local[globalVI] {
                    localF[localTriI, vI] = local
                } else {
                    let localVI: Int = localV_rest.rows
                    if (!freeVert.contains(globalVI)) {
                        fixedVert.insert(localVI)
                    }
                    localV_rest.conservativeResize(localVI + 1, 3)
                    localV_rest.row(localVI) <<== V_rest.row(globalVI)
                    localV.conservativeResize(localVI + 1, 2)
                    localV.row(localVI) <<== V.row(globalVI)
                    localF[localTriI, vI] = localVI
                    globalVI2Local[globalVI] = localVI
                }
            }
            localTriI += 1
        }
        
        let localMesh = TriMesh(localV_rest, localF, localV, Mati(), false)
        localMesh.resetFixedVert(fixedVert)
        
        // compute initial symmetric Dirichlet Energy value
        let SD = SymDirichletEnergy()
        var initE: Double = 0.0
        for triI in triangles {
            var energyValI: Double = 0.0
            SD.getEnergyValByElemID(self, triI, &energyValI)
            initE += energyValI
        }
        initE *= surfaceArea / localMesh.surfaceArea
        
        // split edge
        var UV_bnds = Matd()
        var E = Mati()
        var bnd = Veci()
        var cutThrough: Bool = false
        
        switch splitPath.count {
        case 0:
            // nothing to split
            assert(false, "currently we don't use this function without splitting!")
            break
            
        case 2: // boundary split
            assert(freeVert.contains(splitPath[0]))
            if (freeVert.contains(splitPath[1])) {
                cutThrough = true
            }
            
            // convert splitPath global index to local index
            var splitPath_local = [Int]()
            splitPath_local.reserveCapacity(splitPath.count)
            for pvI in splitPath {
                if let local = globalVI2Local[pvI] {
                    splitPath_local.append(local)
                }
            }
            
            // split
            localMesh.splitEdgeOnBoundary(Pair(splitPath_local[0], splitPath_local[1]), Matd(), false, cutThrough)
            
            if (scaffold != nil) {
                // separate the splitted vertices to leave room for airmesh
                var splittedV: [RVec2d] = [
                    localMesh.V.row(splitPath_local[0]),
                    localMesh.V.row(localMesh.V.rows - 1 - (cutThrough ? 1 : 0))
                ]
                var sepDir_oneV: [RVec2d?] = .init(repeating: nil, count: 2)
                localMesh.compute2DInwardNormal(splitPath_local[0], &sepDir_oneV[0])
                localMesh.compute2DInwardNormal(localMesh.V.rows - 1 - (cutThrough ? 1 : 0), &sepDir_oneV[1])
                var sepDir: [Vec<Double>] = [
                    .Zero(localMesh.V.rows * 2),
                    .Zero(localMesh.V.rows * 2)
                ]
                sepDir[0].block(splitPath_local[0] * 2, 0, 2, 1) <<== sepDir_oneV[0]!.transpose()
                sepDir[1].block((localMesh.V.rows - 1 - (cutThrough ? 1 : 0)) * 2, 0, 2, 1) <<== sepDir_oneV[1]!.transpose()
                let eps_sep: Double = (V.row(splitPath[1]) - V.row(splitPath[0])).squaredNorm() * 1.0e-4
                var curSqDist = (splittedV[0] - splittedV[1]).squaredNorm()
                
                while (curSqDist < eps_sep) {
                    for i in 0..<2 {
                        var stepSize_sep: Double = 1.0
                        SD.initStepSize(localMesh, sepDir[i], &stepSize_sep)
                        splittedV[i] += 0.1 * stepSize_sep * sepDir_oneV[i]!
                    }
                    localMesh.V.row(splitPath_local[0]) <<== splittedV[0]
                    localMesh.V.row(localMesh.V.rows - 1 - (cutThrough ? 1 : 0)) <<== splittedV[1]
                    
                    let lastSqDist: Double = curSqDist
                    curSqDist = (splittedV[0] - splittedV[1]).squaredNorm()
                    if (abs(curSqDist - lastSqDist) / lastSqDist < 1.0e-3) {
                        break
                    }
                }
                assert(localMesh.checkInversion())
                
                if (cutThrough) {
                    splittedV[0] = localMesh.V.row(splitPath_local[1])
                    splittedV[1] = localMesh.V.bottomRows(1)
                    localMesh.compute2DInwardNormal(splitPath_local[1], &sepDir_oneV[0])
                    localMesh.compute2DInwardNormal(localMesh.V.rows - 1, &sepDir_oneV[1])
                    sepDir[1] = .Zero(localMesh.V.rows * 2)
                    sepDir[0] = .Zero(localMesh.V.rows * 2)
                    sepDir[0].block(splitPath_local[1] * 2, 0, 2, 1) <<== sepDir_oneV[0]!.transpose()
                    sepDir[1].bottomRows(2) <<== sepDir_oneV[1]!.transpose()
                    var curSqDist: Double = (splittedV[0] - splittedV[1]).squaredNorm()
                    while (curSqDist < eps_sep) {
                        for i in 0..<2 {
                            var stepSize_sep: Double = 1.0
                            SD.initStepSize(localMesh, sepDir[i], &stepSize_sep)
                            splittedV[i] += 0.1 * stepSize_sep * sepDir_oneV[i]!
                        }
                        localMesh.V.row(splitPath_local[1]) <<== splittedV[0]
                        localMesh.V.bottomRows(1) <<== splittedV[1]
                        
                        let lastSqDist: Double = curSqDist
                        curSqDist = (splittedV[0] - splittedV[1]).squaredNorm()
                        if (abs(curSqDist - lastSqDist) / lastSqDist < 1.0e-3) {
                            break
                        }
                        // may update search dir, and accelerate
                    }
                    assert(localMesh.checkInversion())
                }
                
                // prepare local air mesh boundary
                var UV_temp = Matd()
                var bnd_temp = Vec<Int>()
                var loop_AMVI = Set<Int>()
                scaffold!.get1RingAirLoop(splitPath[0], &UV_temp, &E, &bnd_temp, &loop_AMVI)
                let loopVAmt_beforeSplit: Int = E.rows
                if (!cutThrough) {
                    E.bottomRows(1) <<== [loopVAmt_beforeSplit - 1, loopVAmt_beforeSplit]
                    E.conservativeResize(loopVAmt_beforeSplit + 2, 2)
                    E.bottomRows(2) <<== [loopVAmt_beforeSplit, loopVAmt_beforeSplit + 1, loopVAmt_beforeSplit + 1, 0]
                    UV_bnds.resize(loopVAmt_beforeSplit + 2, 2)
                    UV_bnds.bottomRows(loopVAmt_beforeSplit - 3) <<== UV_temp.bottomRows(loopVAmt_beforeSplit - 3)
                    // NOTE: former vertices will be filled with mesh coordinates while constructing the local airmesh
                    
                    bnd.resize(bnd_temp.count + 2)
                    bnd[0] = bnd_temp[0]
                    bnd[1] = localMesh.V.rows - 1 - (cutThrough ? 1 : 0)
                    bnd[2] = splitPath[1]
                    bnd.bottomRows(2) <<== bnd_temp.bottomRows(2)
                    for bndI in 0..<bnd.count {
                        if (bndI != 1) {
                            if let local = globalVI2Local[bnd[bndI]] {
                                bnd[bndI] = local
                            }
                        }
                    }
                } else {
                    var UV_temp1 = Matd()
                    var bnd_temp1 = Veci()
                    var E1 = Mati()
                    var loop1_AMVI = Set<Int>()
                    scaffold!.get1RingAirLoop(splitPath[1], &UV_temp1, &E1, &bnd_temp1, &loop1_AMVI)
                    // avoid generating air mesh with duplicated vertices
                    // NOTE: this also avoids forming tiny charts
                    for i in loop1_AMVI {
                        if (loop_AMVI.contains(i)) {
                            return -__DBL_MAX__
                        }
                    }
                    let loopVAmt1_beforeSplit: Int = E1.rows
                    
                    UV_bnds.resize(loopVAmt_beforeSplit + loopVAmt1_beforeSplit + 2, 2)
                    UV_bnds.bottomRows(UV_bnds.rows - 8) <<== [UV_temp1.bottomRows(loopVAmt1_beforeSplit - 3), UV_temp.bottomRows(loopVAmt_beforeSplit - 3)]
                    // NOTE: former vertices will be filled with mesh coordinates while constructing the local air mesh
                    
                    bnd.resize(8)
                    bnd[0] = bnd_temp[0]
                    bnd[1] = localMesh.V.rows - 2
                    bnd[2] = splitPath[1]
                    bnd[3] = bnd_temp1[2]
                    bnd[4] = bnd_temp1[0]
                    bnd[5] = localMesh.V.rows - 1
                    bnd[6] = splitPath[0]
                    bnd[7] = bnd_temp[2]
                    for bndI in 0..<bnd.count {
                        if (bndI != 1 && bndI != 5) {
                            if let local = globalVI2Local[bnd[bndI]] {
                                bnd[bndI] = local
                            }
                        }
                    }
                    
                    E.resize(UV_bnds.rows, 2)
                    E.row(0) <<== [0, 1]
                    E.row(1) <<== [1, 2]
                    E.row(2) <<== [2, 3]
                    if (loopVAmt1_beforeSplit - 3 == 0) {
                        E.row(3) <<== [3, 4]
                    } else {
                        E.row(3) <<== [3, 8]
                        for i in 0..<(loopVAmt1_beforeSplit - 3) {
                            E.row(4 + i) <<== [8 + i, 9 + i]
                        }
                        E[loopVAmt1_beforeSplit, 1] = 4
                    }
                    E.row(loopVAmt1_beforeSplit + 1) <<== [4, 5]
                    E.row(loopVAmt1_beforeSplit + 2) <<== [5, 6]
                    E.row(loopVAmt1_beforeSplit + 3) <<== [6, 7]
                    if (loopVAmt_beforeSplit - 3 == 0) {
                        E.row(loopVAmt1_beforeSplit + 4) <<== [7, 0]
                    } else {
                        E.row(loopVAmt1_beforeSplit + 4) <<== [7, loopVAmt1_beforeSplit + 5]
                        for i in 0..<(loopVAmt_beforeSplit - 3) {
                            E.row(loopVAmt1_beforeSplit + 5 + i) <<== [loopVAmt1_beforeSplit + 5 + i, loopVAmt1_beforeSplit + 6 + i]
                        }
                        E[loopVAmt1_beforeSplit + loopVAmt_beforeSplit + 1, 1] = 0
                    }
                }
            }
        case 3: // interior split
            // not processed here
            break
            
        default:
            assert(false, "invalid split path!")
            break
        }
        
        // conduct optimization on local mesh
        let energyTerms: [Energy] = [SD]
        let param: UnsafeMutablePointer<Double> = .allocate(capacity: 1)
        param.initialize(to: 1.0)
        let energyParams: [UnsafeMutablePointer<Double>] = [param]
        defer {
            param.deallocate()
        }
        
        let optimizer = Optimizer(localMesh, energyTerms, energyParams, 0, true, scaffold != nil, UV_bnds, E, bnd, true)
        optimizer.precompute()
        
        optimizer.setRelGL2Tol(1.0e-6)
        optimizer.solve(maxIter)
        
        var curE: Double = 0.0
        optimizer.computeEnergyVal(optimizer.getResult(), optimizer.getScaffold(), &curE, true)
        let eDec: Double = (initE - curE) * localMesh.surfaceArea / surfaceArea
        
        // get new vertex positions
        newVertPos.resize(2, 2)
        newVertPos <<== [optimizer.getResult().V.row(globalVI2Local[splitPath[0]]!),
                         optimizer.getResult().V.row(localMesh.V.rows - 1 - (cutThrough ? 1 : 0))]
        if (cutThrough) {
            newVertPos.conservativeResize(4, 2)
            newVertPos.row(2) <<== optimizer.getResult().V.bottomRows(1)
            newVertPos.row(3) <<== optimizer.getResult().V.row(globalVI2Local[splitPath[1]]!)
        }
        
        return eDec
    }
    
    // interior split
     func computeLocalEdDec_inSplit(_ triangles: [Int],
                                          _ freeVert: Set<Int>,
                                          _ path: [Int],
                                          _ newVertPos: inout Mat<Double>,
                                          _ maxIter: Int = 100) -> Double {
        
        assert(!triangles.isEmpty && !freeVert.isEmpty)
        
        // construct local mesh
        var localF: Mati = .init(triangles.count, 3)
        var localV_rest = Matd()
        var localV = Matd()
        var fixedVert = Set<Int>()
        var globalVI2local = [Int:Int]()
        var localTriI: Int = 0
        for triI in triangles {
            for vI in 0..<3 {
                let globalVI: Int = F[triI, vI]
                if let localVI = globalVI2local[globalVI] {
                    localF[localTriI, vI] = localVI
                } else {
                    let localVI: Int = localV_rest.rows
                    if (!freeVert.contains(globalVI)) {
                        fixedVert.insert(localVI)
                    }
                    localV_rest.conservativeResize(localVI + 1, 3)
                    localV_rest.row(localVI) <<== V_rest.row(globalVI)
                    localV.conservativeResize(localVI + 1, 2)
                    localV.row(localVI) <<== V.row(globalVI)
                    localF[localTriI, vI] = localVI
                    globalVI2local[globalVI] = localVI
                }
            }
            localTriI += 1
        }
        let localMesh = TriMesh(localV_rest, localF, localV, Mati(), false)
        localMesh.resetFixedVert(fixedVert)
        
        let SD = SymDirichletEnergy()
        var initE: Double = 0.0
        for triI in triangles {
            var energyValI: Double = 0.0
            SD.getEnergyValByElemID(self, triI, &energyValI)
            initE += energyValI
        }
        initE *= (surfaceArea / localMesh.surfaceArea)
        
        // convert split path global index to local index
        var path_local: [Int] = []
        path_local.reserveCapacity(path.count)
        for pvI in path {
            if let local = globalVI2local[pvI] {
                path_local.append(local)
            }
        }
        // split
        localMesh.cutPath(path_local, true, 0, Matd(), false)
        
        let isBijective: Bool = (scaffold != nil) // TODO: verify this line is correctly converted from C++
        
        // construct air mesh
        var UV_bnds: Matd = .init()
        var E: Mati = .init()
        var bnd: Veci = .init()
        if (isBijective) {
            // separate vertex
            // TODO: write into a function
            var splittedV: [RVec2d] = [
                localMesh.V.row(path_local[1]),
                localMesh.V.row(localMesh.V.rows - 1)
            ]
            var sepDir_oneV: [RVec2d?] = .init(repeating: nil, count: 2)
            localMesh.compute2DInwardNormal(path_local[1], &sepDir_oneV[0])
            localMesh.compute2DInwardNormal(localMesh.V.rows - 1, &sepDir_oneV[1])
            let sepDir: [Vec<Double>] = [
                .Zero(localMesh.V.rows * 2),
                .Zero(localMesh.V.rows * 2)
            ]
            sepDir[0].block(path_local[1] * 2, 0, 2, 1) <<== sepDir_oneV[0]!.transpose()
            sepDir[1].block((localMesh.V.rows - 1) * 2, 0, 2, 1) <<== sepDir_oneV[1]!.transpose()
            let eps_sep: Double = (V.row(path[1]) - V.row(path[0])).squaredNorm() * 1.0e-4
            var curSqDist: Double = (splittedV[0] - splittedV[1]).squaredNorm()
            while (curSqDist < eps_sep) {
                for i in 0..<2 {
                    var stepSize_sep: Double = 1.0
                    SD.initStepSize(localMesh, sepDir[i], &stepSize_sep)
                    splittedV[i] += 0.1 * stepSize_sep * sepDir_oneV[i]!
                }
                localMesh.V.row(path_local[1]) <<== splittedV[0]
                localMesh.V.row(localMesh.V.rows - 1) <<== splittedV[1]
                
                let lastSqDist: Double = curSqDist
                curSqDist = (splittedV[0] - splittedV[1]).squaredNorm()
                if (abs(curSqDist - lastSqDist) / lastSqDist < 1.0e-3) {
                    break
                }
            }
            
            // establish air mesh information
            UV_bnds.resize(4, 2)
            E = .init([0, 1, 1, 2, 2, 3, 3, 0], [4, 2])
            bnd = .init([path_local[2], path_local[1], path_local[0], localMesh.V.rows - 1])
        }
        
        // conduct optimization on local mesh
        let energyTerms: [Energy] = [SD]
        let param: UnsafeMutablePointer<Double> = .allocate(capacity: 1)
        param.initialize(to: 1.0)
        let energyParams: [UnsafeMutablePointer<Double>] = [param]
        defer {
            param.deallocate()
        }
        let optimizer = Optimizer(localMesh, energyTerms, energyParams, 0, true, isBijective, UV_bnds, E, bnd, true)
        optimizer.precompute()
        optimizer.setRelGL2Tol(1.0e-6)
        optimizer.solve(maxIter) // do not output, the other part
        var curE: Double = 0.0
        optimizer.computeEnergyVal(optimizer.getResult(), optimizer.getScaffold(), &curE, true)
        let eDec: Double = (initE - curE) * localMesh.surfaceArea / surfaceArea
        
        // get new vertex positions
        newVertPos = .init(2, 2)
        newVertPos.row(0) <<== optimizer.getResult().V.bottomRows(1)
        newVertPos.row(1) <<== optimizer.getResult().V.row(path_local[1])
        
        return eDec
    }
    
    // merge
     func computeLocalEdDec_merge(_ path: [Int],
                                        _ triangles: [Int],
                                        _ freeVert: Set<Int>,
                                        _ newVertPos: inout [Int: RVec2d],
                                        _ mergeVert: [Int: Int],
                                        _ initMergedPos: RVec2d,
                                        _ closeup: Bool = false,
                                        _ maxIter: Int = 100) -> Double {
        
        assert(!triangles.isEmpty && !freeVert.isEmpty)
        assert(!mergeVert.isEmpty)
        
        let isBijective: Bool = ((scaffold != nil) && (!closeup))
        
        // construct local mesh
        var localF = Mati(triangles.count, 3)
        var localV_rest = Matd()
        var localV = Matd()
        var fixedVert = Set<Int>()
        var globalVI2Local = [Int:Int]()
        var localTriI: Int = 0
        for triI in triangles {
            for vI in 0..<3 {
                let globalVI = F[triI, vI]
                if let merge = mergeVert[globalVI] {
                    // one of the vertices to be merged
                    let localVIFinder = globalVI2Local[globalVI]
                    let localVIFinder_mergePair = globalVI2Local[merge]
                    let selfAdded: Bool = (localVIFinder != nil)
                    let mergePairAdded: Bool = (localVIFinder_mergePair != nil)
                    if (selfAdded) {
                        assert(mergePairAdded)
                        localF[localTriI, vI] = localVIFinder!
                    } else {
                        assert(!mergePairAdded)
                        let localVI: Int = localV_rest.rows
                        if (!freeVert.contains(globalVI)) {
                            fixedVert.insert(localVI)
                        }
                        localV_rest.conservativeResize(localVI + 1, 3)
                        localV_rest.row(localVI) <<== V_rest.row(globalVI)
                        localV.conservativeResize(localVI + 1, 2)
                        localV.row(localVI) <<== initMergedPos
                        localF[localTriI, vI] = localVI
                        globalVI2Local[globalVI] = localVI
                        
                        globalVI2Local[merge] = localVI
                    }
                } else {
                    // normal vertices
                    if let local = globalVI2Local[globalVI] {
                        localF[localTriI, vI] = local
                    } else {
                        let localVI: Int = localV_rest.rows
                        if (!freeVert.contains(globalVI)) {
                            fixedVert.insert(localVI)
                        }
                        localV_rest.conservativeResize(localVI + 1, 3)
                        localV_rest.row(localVI) <<== V_rest.row(globalVI)
                        localV.conservativeResize(localVI + 1, 2)
                        localV.row(localVI) <<== V.row(globalVI)
                        localF[localTriI, vI] = localVI
                        globalVI2Local[globalVI] = localVI
                    }
                }
            }
            localTriI += 1
        }
        var localMesh = TriMesh(localV_rest, localF, localV, Mati(), false)
        localMesh.resetFixedVert(fixedVert)
        
        let SD = SymDirichletEnergy()
        var initE: Double = 0.0
        for triI in triangles {
            var energyValI: Double = 0.0
            SD.getEnergyValByElemID(self, triI, &energyValI)
            initE += energyValI
        }
        initE *= surfaceArea / localMesh.surfaceArea
        
        // construct air mesh
        var UV_bnds = Matd()
        var E = Mati()
        var bnd = Veci()
        
        if (isBijective) {
            if (!scaffold!.getCornerAirLoop(path, initMergedPos, &UV_bnds, &E, &bnd)) {
                // if initPos causes the composite loop to self-intersect, or the loop is totally inverted
                // (potentially violating bijectivity), abondon this query
                return -__DBL_MAX__
            }
            
            for bndI in 0..<bnd.count {
                if let local = globalVI2Local[bnd[bndI]] {
                    bnd[bndI] = local
                }
            }
        }
        
        // conduct optimization on local mesh
        let energyTerms: [Energy] = [SD]
        let param: UnsafeMutablePointer<Double> = .allocate(capacity: 1)
        param.initialize(to: 1.0)
        let energyParams: [UnsafeMutablePointer<Double>] = [param]
        defer {
            param.deallocate()
        }
        let optimizer = Optimizer(localMesh, energyTerms, energyParams, 0, true, isBijective, UV_bnds, E, bnd, true)
        optimizer.precompute()
        optimizer.setRelGL2Tol(1.0e-6)
        optimizer.solve(maxIter) // do not output, the other part
        var curE: Double = 0
        optimizer.computeEnergyVal(optimizer.getResult(), optimizer.getScaffold(), &curE, true)
        let eDec: Double = (initE - curE) * localMesh.surfaceArea / surfaceArea
        
        // get new Vertex positions
        newVertPos.removeAll()
        for vI_free in freeVert {
            newVertPos[vI_free] = optimizer.getResult().V.row(globalVI2Local[vI_free]!)
        }
        
        return eDec
    }
    /*
     func clone() -> TriMesh {
        let clone = TriMesh()
        clone.V_rest = .init(V_rest, V_rest.rows, V_rest.cols)
        clone.V = .init(V, V.rows, V.cols)
        clone.F = .init(F, F.rows, F.cols)
        clone.cohE = .init(cohE, cohE.rows, cohE.cols)
        clone.initSeams = .init(initSeams, initSeams.rows, initSeams.cols)
        
        clone.scaffold = scaffold
        clone.areaThres_AM = areaThres_AM
        
        clone.boundaryEdge = .init(boundaryEdge, boundaryEdge.rows, boundaryEdge.cols)
        clone.edgeLen = .init(edgeLen, edgeLen.rows, edgeLen.cols)
        clone.LaplacianMtr = LaplacianMtr
        clone.triArea = .init(triArea, triArea.rows, triArea.cols)
        clone.triNormal = .init(triNormal, triNormal.rows, triNormal.cols)
        clone.surfaceArea = surfaceArea
        clone.triAreaSq = .init(triAreaSq, triAreaSq.rows, triAreaSq.cols)
        clone.e0dote1 = .init(e0dote1, e0dote1.rows, e0dote1.cols)
        clone.e0SqLen = .init(e0SqLen, e0SqLen.rows, e0SqLen.cols)
        clone.e1SqLen = .init(e1SqLen, e1SqLen.rows, e1SqLen.cols)
        clone.e0SqLen_div_dbAreaSq = .init(e0SqLen_div_dbAreaSq, e0SqLen_div_dbAreaSq.rows, e0SqLen_div_dbAreaSq.cols)
        clone.e1SqLen_div_dbAreaSq = .init(e1SqLen_div_dbAreaSq, e1SqLen_div_dbAreaSq.rows, e1SqLen_div_dbAreaSq.cols)
        clone.e0dote1_div_dbAreaSq = .init(e0dote1_div_dbAreaSq, e0dote1_div_dbAreaSq.rows, e0dote1_div_dbAreaSq.cols)
        clone.avgEdgeLen = avgEdgeLen
        clone.virtualRadius = virtualRadius
        clone.validSplit = validSplit
        clone.bbox = .init(bbox, bbox.rows, bbox.cols)
        clone.vertWeight = .init(vertWeight, vertWeight.rows, vertWeight.cols)
        
        clone.edge2Tri = edge2Tri
        clone.vNeighbor = vNeighbor
        clone.cohEIndex = cohEIndex
        
        clone.fracTail = fracTail
        clone.curFracTail = curFracTail
        clone.curInteriorFracTails = curInteriorFracTails
        clone.initSeamLen = initSeamLen
        
        return clone
    }*/
}

private func initCylinder<MV: Matrix, MF: Matrix>
(_ r1_x: Double, _ r1_y: Double, _ r2_x: Double, _ r2_y: Double, _ height: Double, _ circle_res: Int, _ height_resolution: Int,
 _ V: inout MV,
 _ F: inout MF,
 _ uv_coords_per_face: UnsafeMutablePointer<Mat<Double>>? = nil,
 _ uv_coords_face_ids: UnsafeMutablePointer<Mat<Int>>? = nil)
where MV.Element == Double, MF.Element == Int
{
    let nvertices: Int = circle_res * (height_resolution + 1)
    let nfaces: Int = 2 * circle_res * height_resolution
    
    V.resize(nvertices, 3)
    if (uv_coords_per_face != nil) {
        uv_coords_per_face!.pointee.resize(nvertices, 2)
    }
    F.resize(nfaces, 3)
    for j in 0..<(height_resolution + 1) {
        for i in 0..<circle_res {
            let t: Double = Double(j) / Double(height_resolution)
            let h: Double = height * t
            let theta: Double = Double(i) * 2 * .pi / Double(circle_res)
            let r_x: Double = r1_x * t + r2_x * (1 - t)
            let r_y: Double = r1_y * t + r2_y * (1 - t)
            V.row(j * circle_res + i) <<== Vec3d(r_x * cos(theta), height - h, r_y * sin(theta))
            if (uv_coords_per_face != nil) {
                uv_coords_per_face!.pointee.row(j * circle_res + i) <<== Vec2d(r_x * cos(theta), r_y * sin(theta))
            }
            
            if (j < height_resolution) {
                let vl0: Int = j * circle_res + i
                let vl1: Int = j * circle_res + ((i + 1) % circle_res)
                let vu0: Int = (j + 1) * circle_res + 1
                let vu1: Int = (j + 1) * circle_res + ((i + 1) % circle_res)
                F.row(2 * (j * circle_res + i) + 0) <<== Vec3i(vl0, vl1, vu1)
                F.row(2 * (j * circle_res + i) + 1) <<== Vec3i(vu0, vl0, vu1)
            }
        }
    }
}

// A utility function to find the vertex with minimum distance value, from
// the set of vertices not yet included in shortest path tree
private func minDistance(_ dist: [Double], _ sptSet: [Bool]) -> Int {
    // Initialize min value
    var min: Double = __DBL_MAX__
    var min_index: Int = -1
    
    for v in 0..<dist.count {
        if (!sptSet[v] && dist[v] <= min) {
            min = dist[v]
            min_index = v
        }
    }
    
    return min_index
}

// Function that implements Dijkstra's single source shortest path algorithm
// for a graph represented using adjacency matrix representation
func dijkstra(_ graph: [[Int:Double]],
              _ src: Int,
              _ dist: inout [Double]!,
              _ parent: inout [Int]!) {
    
    let nV: Int = graph.count
    
    dist = .init(repeating: __DBL_MAX__, count: nV)
    
    // The output array. dist[i] will hold the shortest
    // distance from src to i
    
    var sptSet: [Bool] = .init(repeating: false, count: nV) // sptSet[i] will true if vertex i is included in shortest
    // path tree or shortest distance from src to i is finalized
    
    parent = .init(repeating: -1, count: nV)
    
    // Distance of source vertex from itself is always 0
    dist[src] = 0.0
    
    // Find shortest path for all vertices
    for _ in 0..<(nV - 1) {
        // Pick the minimum distance vertex from the set of vertices not
        // yet processed. u is always equal to src in first iteration
        var u: Int = minDistance(dist, sptSet)
        
        // Mark the picked vertex as processed
        sptSet[u] = true
        
        for (key, value) in graph[u] {
            // Update dist[v] only if is not in sptSet, there is an edge from
            // u to v, and total weight of path from src to v through u is
            // smaller than current value of dist[v]
            if ((!sptSet[key])
                && (dist[u] != __DBL_MAX__)
                && (dist[u] + value < dist[key])) {
                dist[key] = dist[u] + value
                parent[key] = u
            }
        }
    }
}

private func getFarthestPoint(_ graph: [[Int:Double]], _ src: Int) -> Int {
    let nV: Int = graph.count
    var dist: [Double]!
    var parent: [Int]!
    dijkstra(graph, src, &dist, &parent)
    
    var maxDist: Double = 0.0
    var vI_maxDist: Int = -1
    for vI in 0..<nV {
        if (dist[vI] > maxDist) {
            maxDist = dist[vI]
            vI_maxDist = vI
        }
    }
    assert(vI_maxDist >= 0)
    return vI_maxDist
}

private func getFarthestPointPath(_ graph: [[Int:Double]], _ src: Int, _ path: inout [Int]!) {
    let nV: Int = graph.count
    var dist: [Double]!
    var parent: [Int]!
    dijkstra(graph, src, &dist, &parent)
    
    var maxDist: Double = 0.0
    var vI_maxDist: Int = -1
    for vI in 0..<nV {
        if (dist[vI] > maxDist) {
            maxDist = dist[vI]
            vI_maxDist = vI
        }
    }
    assert(vI_maxDist >= 0)
    path = []
    while (vI_maxDist >= 0) {
        path.append(vI_maxDist)
        vI_maxDist = parent[vI_maxDist]
    }
    path.reverse()
}
