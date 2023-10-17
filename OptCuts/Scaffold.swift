//
//  Scaffold.swift
//  OptCuts
//
//  Created by Reza on 6/20/23.
//

import Foundation
import Matrix
import GeometryProcessing

 class Scaffold {
    // MARK: - Properties
    var airMesh: TriMesh! // tesselation of voided regions
    var bnd: Vec<Int> = .init()
    var localVI2Global: Vec<Int> = .init() // map between airMesh indices to augmented system indices
    var meshVI2AirMesh: [Int: Int] = [:] // the inverse map of bnd
    var wholeMeshSize: Int = 0 // augmented system size
    
    // MARK: - Initialization
     init() {
        
    }
    
     init(_ mesh: TriMesh,
                _ UV_bnds: Mat<Double> = Matd(),
                _ E: Mat<Int> = Mati(),
                _ p_bnd: Vec<Int> = Veci()) {
        
        airMesh = .init()
        
        var UV_bnds: Mat<Double> = UV_bnds
        var E: Mat<Int> = E
        
        assert(E.rows == UV_bnds.rows)
        
        var H = Matd()
        var fixAMBoundary: Bool = false
        var edgeLen_eps: Double = mesh.avgEdgeLen * 0.5 // NOTE: different from what's used in [Jiang et al. 2017]
        if (E.rows == 0) {
            var bnd_all: [[Int]] = []
            boundary_loop(mesh.F, &bnd_all)
            assert(bnd_all.count != 0)
            
            var curBndVAmt: Int = 0
            for bndI in 0..<bnd_all.count {
                UV_bnds.conservativeResize(curBndVAmt + bnd_all[bndI].count, 2)
                E.conservativeResize(curBndVAmt + bnd_all[bndI].count, 2)
                bnd.conservativeResize(curBndVAmt + bnd_all[bndI].count)
                for bndVI in 0..<(bnd_all[bndI].count - 1) {
                    E.row(curBndVAmt + bndVI) <<== [curBndVAmt + bndVI, curBndVAmt + bndVI + 1]
                    UV_bnds.row(curBndVAmt + bndVI) <<== mesh.V.row(bnd_all[bndI][bndVI])
                    bnd[curBndVAmt + bndVI] = bnd_all[bndI][bndVI]
                    meshVI2AirMesh[bnd[curBndVAmt + bndVI]] = curBndVAmt + bndVI
                }
                E.row(curBndVAmt + bnd_all[bndI].count - 1) <<== [curBndVAmt + bnd_all[bndI].count - 1, curBndVAmt]
                UV_bnds.row(curBndVAmt + bnd_all[bndI].count - 1) <<== mesh.V.row(bnd_all[bndI].last!)
                bnd[curBndVAmt + bnd_all[bndI].count - 1] = bnd_all[bndI].last!
                meshVI2AirMesh[bnd[curBndVAmt + bnd_all[bndI].count - 1]] = curBndVAmt + bnd_all[bndI].count - 1
                curBndVAmt = E.rows
            }
            
            let scaleFactor: Double = 3.0 // >= 2.0 is recommended
            let bandAmt: Int = 2 // >= 2 is recommended
            let segLen: Double = mesh.avgEdgeLen * pow(scaleFactor, Double(bandAmt))
            let margin: Double = (mesh.avgEdgeLen - segLen) / (1.0 - scaleFactor) // TODO: use boundary edge lengths of current UV map?
            let minX: Double = mesh.V.col(0).minCoeff() - margin
            let maxX: Double = mesh.V.col(0).maxCoeff() + margin
            let minY: Double = mesh.V.col(1).minCoeff() - margin
            let maxY: Double = mesh.V.col(1).maxCoeff() + margin
            // segment the bounding box:
            let segAmtX: Int = Int(ceil((maxX - minX) / segLen))
            let segAmtY: Int = Int(ceil((maxY - minY) / segLen))
            E.conservativeResize(E.rows + (segAmtX + segAmtY) * 2, 2)
            for segI in (bnd.count)..<(E.rows - 1) {
                E.row(segI) <<== [segI, segI + 1]
            }
            E.row(E.rows - 1) <<== [E.rows - 1, bnd.count]
            let stepX: Double = (maxX - minX) / Double(segAmtX)
            let stepY: Double = (maxY - minY) / Double(segAmtY)
            UV_bnds.conservativeResize(UV_bnds.rows + (segAmtX + segAmtY) * 2, 2)
            for segI in 0..<segAmtX {
                UV_bnds.row(bnd.count + segI) <<== [minX + Double(segI) * stepX, minY]
            }
            for segI in 0..<segAmtY {
                UV_bnds.row(bnd.count + segAmtX + segI) <<== [maxX, minY + Double(segI) * stepY]
            }
            for segI in 0..<segAmtX {
                UV_bnds.row(bnd.count + segAmtX + segAmtY + segI) <<== [maxX - Double(segI) * stepX, maxY]
            }
            for segI in 0..<segAmtY {
                UV_bnds.row(bnd.count + 2 * segAmtX + segAmtY + segI) <<== [minX, maxY - Double(segI) * stepY]
            }
            
            // compute connected component of mesh
            var compI_V = Veci()
            vertex_components(mesh.F, &compI_V)
            // mark holes
            var processedComp = Set<Int>()
            for vI in 0..<compI_V.count {
                if (!processedComp.contains(compI_V[vI])) {
                    H.conservativeResize(H.rows + 1, 2)
                    
                    var incTris = [Int]()
                    var bEdge: Pair<Int, Int> = Pair(0, 0)
                    if (mesh.isBoundaryVert(vI, mesh.vNeighbor[vI].first!, &incTris, &bEdge, false)) {
                        // push vI a libble bit inside mesh and then add into H
                        // get all incident triangles
                        var temp: [Int] = []
                        mesh.isBoundaryVert(vI, mesh.vNeighbor[vI].first!, &temp, &bEdge, true)
                        incTris.append(contentsOf: temp)
                        // construct local mesh
                        var localF = Mati(incTris.count, 3)
                        var localV_rest = Matd()
                        var localV = Matd()
                        var globalVI2local: [Int:Int] = [:]
                        var localTriI: Int = 0
                        for triI in incTris {
                            for vI in 0..<3 {
                                let globalVI = mesh.F[triI, vI]
                                if let value = globalVI2local[globalVI] {
                                    localF[localTriI, vI] = value
                                } else {
                                    let localVI: Int = localV_rest.rows
                                    localV_rest.conservativeResize(localVI + 1, 3)
                                    localV_rest.row(localVI) <<== mesh.V_rest.row(globalVI)
                                    localV.conservativeResize(localVI + 1, 2)
                                    localV.row(localVI) <<== mesh.V.row(globalVI)
                                    localF[localTriI, vI] = localVI
                                    globalVI2local[globalVI] = localVI
                                }
                            }
                            localTriI += 1
                        }
                        let localMesh = TriMesh(localV_rest, localF, localV, Mati(), false)
                        // compute inward normal
                        var sepDir_oneV: RVec2d!
                        mesh.compute2DInwardNormal(vI, &sepDir_oneV)
                        let sepDir: Vec<Double> = .Zero(localMesh.V.rows * 2)
                        sepDir.block(globalVI2local[vI]! * 2, 0, 2, 1) <<== sepDir_oneV.transpose()
                        var stepSize_sep: Double = 1.0
                        let SD = SymDirichletEnergy()
                        SD.initStepSize(localMesh, sepDir, &stepSize_sep)
                        H.bottomRows(1) <<== mesh.V.row(vI) + 0.5 * stepSize_sep * sepDir_oneV
                    } else {
                        H.bottomRows(1) <<== mesh.V.row(vI)
                    }
                    
                    processedComp.insert(compI_V[vI])
                }
            }
        } else {
            assert(p_bnd.rows > 0)
            
            edgeLen_eps *= 0.1
            fixAMBoundary = true
            bnd = .init(p_bnd, p_bnd.rows, p_bnd.cols)
            for bndI in 0..<bnd.count {
                UV_bnds.row(bndI) <<== mesh.V.row(bnd[bndI])
                meshVI2AirMesh[bnd[bndI]] = bndI
            }
            
            //[NOTE] this option is for optimization on local stencil
            // where there shouldn't be holes on the one-ring,
            // so no processing for H
        }
        triangulate(UV_bnds, E, H, "qYQ", &airMesh.V, &airMesh.F)
        // "Y" for no Steiner points on mesh boundary
        // "q" for high quality mesh generation
        // "Q" for qiet mode (no output)
        
        airMesh.V_rest.resize(airMesh.V.rows, 3)
        airMesh.V_rest <<== [airMesh.V, Mat<Double>.Zero(airMesh.V.rows, 1)]
        airMesh.areaThres_AM = sqrt(3.0) / 4.0 * edgeLen_eps * edgeLen_eps // for preventing degenerate air mesh triangles
        airMesh.computeFeatures()
        
        localVI2Global = .init(bnd, bnd.rows, bnd.cols)
        localVI2Global.conservativeResize(airMesh.V.rows)
        for vI in bnd.count..<airMesh.V.rows {
            localVI2Global[vI] = mesh.V.rows + vI - bnd.count
        }
        wholeMeshSize = mesh.V.rows + airMesh.V.rows - bnd.count
        
        // mark fixed boundary vertices if any in air mesh
        var count: Int = 0
        for meshFixedVI in mesh.fixedVert {
            if let vertInd = meshVI2AirMesh[meshFixedVI] {
                count += 1
                airMesh.fixedVert.insert(vertInd)
            }
        }
        
        if (fixAMBoundary) {
            // fix air mesh boundary vertices
            for vI in bnd.count..<UV_bnds.rows {
                airMesh.fixedVert.insert(vI)
            }
        }
        // add bijectivity to optimization on local stencil: interior split?
        // filter out operations that will cause overlap initially
        // fix bounding box?
    }
    
    // augment mesh gradient with air mesh gradient with parameter w_scaf
     func augmentGradient(_ gradient: inout Vec<Double>,
                         _ gradient_scaf: Vec<Double>,
                         _ w_scaf: Double) {
        
        assert(gradient.count / 2 + airMesh.V.rows - bnd.count == wholeMeshSize)
        assert(w_scaf >= 0.0)
        
        let systemSize0: Int = gradient.count
        gradient.conservativeResize(wholeMeshSize * 2)
        if (w_scaf == 0.0) {
            // for line search
            gradient.bottomRows(gradient.count - systemSize0) <<== Vec<Double>.Zero(gradient.count - systemSize0)
        } else {
            for vI in 0..<bnd.count {
                gradient[bnd[vI] * 2] += w_scaf * gradient_scaf[vI * 2]
                gradient[bnd[vI] * 2 + 1] += w_scaf * gradient_scaf[vI * 2 + 1]
            }
            if (gradient_scaf.count - bnd.count * 2 != 0) {
                let tmp: Vec<Double> = w_scaf * gradient_scaf.bottomRows(gradient_scaf.count - bnd.count * 2)
                gradient.bottomRows(gradient.count - systemSize0) <<== tmp
            }
        }
    }
    
    // augment mesh proxy matrix with air mesh proxy matrix with parameter w_scaf
     func augmentProxyMatrix(_ I: inout Vec<Int>,
                            _ J: inout Vec<Int>,
                            _ V: inout Vec<Double>,
                            _ I_scaf: Vec<Int>,
                            _ J_scaf: Vec<Int>,
                            _ V_scaf: Vec<Double>,
                            _ w_scaf: Double) {
        
        assert(w_scaf > 0.0)
        let tupleSize0: Int = I.count
        
        V.conservativeResize(tupleSize0 + V_scaf.count)
        V.bottomRows(V_scaf.count) <<== w_scaf * V_scaf
        I.conservativeResize(tupleSize0 + I_scaf.count)
        J.conservativeResize(tupleSize0 + J_scaf.count)
        
        for tupleI in 0..<I_scaf.count {
            I[tupleI + tupleSize0] = localVI2Global[I_scaf[tupleI] / 2] * 2 + I_scaf[tupleI] % 2
            J[tupleI + tupleSize0] = localVI2Global[J_scaf[tupleI] / 2] * 2 + J_scaf[tupleI] % 2
        }
    }
    
    // when using dense representation
     func augmentProxyMatrix(_ P: inout Mat<Double>,
                                   _ P_scaf: Mat<Double>,
                                   _ w_scaf: Double) {
        
        assert(w_scaf > 0.0)
         
        let P_oldRows: Int = P.rows
        P.conservativeResize(wholeMeshSize * 2, wholeMeshSize * 2)
        P.rightCols(P.rows - P_oldRows).setZero()
        P.block(P_oldRows, 0, P.rows - P_oldRows, P_oldRows).setZero()
        for localVI in 0..<localVI2Global.count {
            let _2globalVI: Int = localVI2Global[localVI] * 2
            for localVJ in 0..<localVI2Global.count {
                let _2globalVJ: Int = localVI2Global[localVJ] * 2
                let tmp: Matd = w_scaf * P_scaf.block(localVI * 2, localVJ * 2, 2, 2)
                P.block(_2globalVI, _2globalVJ, 2, 2) += tmp
            }
        }
    }
    
    // extract air mesh searchDir from augmented searchDir
     func wholeSearchDir2airMesh(_ searchDir: Vec<Double>,
                                       _ searchDir_airMesh: inout Vec<Double>) {
        
        assert(searchDir.count / 2 == wholeMeshSize)
        
        searchDir_airMesh.resize(airMesh.V.rows * 2)
        for vI in 0..<bnd.count {
            searchDir_airMesh[vI * 2] = searchDir[bnd[vI] * 2]
            searchDir_airMesh[vI * 2 + 1] = searchDir[bnd[vI] * 2 + 1]
        }
        
        let airSystemSize = (airMesh.V.rows - bnd.count) * 2
        searchDir_airMesh.bottomRows(airSystemSize) <<== searchDir.bottomRows(airSystemSize)
    }
    
    // stepForward air mesh using augmented searchDir
     func stepForward(_ V0: Mat<Double>,
                     _ searchDir: Vec<Double>,
                     _ stepSize: Double) {
        
        assert(searchDir.count / 2 == wholeMeshSize)
        assert(V0.rows == airMesh.V.rows)
        
        var searchDir_airMesh = Vec<Double>()
        wholeSearchDir2airMesh(searchDir, &searchDir_airMesh)
        for vI in 0..<airMesh.V.rows {
            airMesh.V[vI, 0] = V0[vI, 0] + stepSize * searchDir_airMesh[vI * 2]
            airMesh.V[vI, 1] = V0[vI, 1] + stepSize * searchDir_airMesh[vI * 2 + 1]
        }
    }
    
     func mergeVNeighbor(_ vNeighbor_mesh: [Set<Int>],
                        _ vNeighbor: inout [Set<Int>]) {
         vNeighbor.removeAll(keepingCapacity: true)
         vNeighbor.append(contentsOf: vNeighbor_mesh)
        //vNeighbor = vNeighbor_mesh
        vNeighbor.conservativeResize(to: wholeMeshSize)
        // TODO: Make parallel
        
         
        for scafVI in 0..<airMesh.vNeighbor.count {
            for nb_scafVI in airMesh.vNeighbor[scafVI] {
                vNeighbor[localVI2Global[scafVI]].insert(localVI2Global[nb_scafVI])
            }
        }
    }
    
     func mergeFixedV(_ fixedV_mesh: Set<Int>,
                     _ fixedV: inout Set<Int>) {
        
        fixedV = fixedV_mesh
        for fixedV_scafVI in airMesh.fixedVert {
            fixedV.insert(localVI2Global[fixedV_scafVI])
        }
    }
    
    // for rendering purposes
     func augmentUVwithAirMesh(_ UV: inout Mat<Double>,
                              _ scale: Double) {
        
        let meshSize: Int = UV.rows
        UV.conservativeResize(wholeMeshSize, 2)
        let tmp: Matd = scale * airMesh.V.bottomRows(airMesh.V.rows - bnd.count)
        UV.bottomRows(wholeMeshSize - meshSize) <<== tmp
    }
    
     func augmentFwithAirMesh(F: inout Mat<Int>) {
        let meshFAmt = F.rows
        F.conservativeResize(meshFAmt + airMesh.F.rows, 3)
        for fI in meshFAmt..<F.rows {
            for vI in 0..<3 {
                F[fI, vI] = localVI2Global[airMesh.F[fI - meshFAmt, vI]]
            }
        }
    }
    
     func augmentFColorwithAirMesh(FColor: inout Mat<Double>) {
        FColor.conservativeResize(FColor.rows + airMesh.F.rows, 3)
        FColor.bottomRows(airMesh.F.rows) <<== Matd.Ones(airMesh.F.rows, 3)
    }
    
    // get 1-ring aurmesh loop for scaffolding optimization on local stencils
     func get1RingAirLoop(_ vI: Int,
                                _ UV: inout Mat<Double>,
                                _ E: inout Mat<Int>,
                                _ bnd: inout Vec<Int>,
                                _ loop_AMVI: inout Set<Int>) {
        
        guard let value = meshVI2AirMesh[vI] else { fatalError() }
        assert(airMesh.isBoundaryVert(value))
        
        var umbrella0 = [Int]()
        var umbrella = [Int]()
        var boundaryEdge: Pair<Int, Int> = Pair(0, 0)
        let nbV0I: Int = airMesh.vNeighbor[value].first!
        airMesh.isBoundaryVert(value, nbV0I, &umbrella0, &boundaryEdge, false)
        airMesh.isBoundaryVert(value, nbV0I, &umbrella, &boundaryEdge, true)
        if (!umbrella.isEmpty) {
            umbrella.reverse()
        }
        if (!umbrella0.isEmpty) {
            umbrella.append(contentsOf: umbrella0)
        }
        assert(!umbrella.isEmpty)
        
        loop_AMVI.removeAll()
        UV.resize(umbrella.count + 2, 2)
        E.resize(umbrella.count + 2, 2)
        bnd.resize(3)
        for tI in 0..<umbrella.count {
            let triI = umbrella[tI]
            for vI in 0..<3 {
                if (airMesh.F[triI, vI] == value) {
                    if (tI == 0) {
                        UV.row(1) <<== airMesh.V.row(value)
                        loop_AMVI.insert(value)
                        E.row(0) <<== [0, 1]
                        bnd[1] = self.bnd[value]
                        bnd[2] = self.bnd[airMesh.F[triI, (vI + 1) % 3]]
                    }
                    UV.row(tI + 2) <<== airMesh.V.row(airMesh.F[triI, (vI + 1) % 3])
                    loop_AMVI.insert(airMesh.F[triI, (vI + 1) % 3])
                    E.row(tI + 1) <<== [tI + 1, tI + 2]
                    if (tI + 1 == umbrella.count) {
                        UV.row(0) <<== airMesh.V.row(airMesh.F[triI, (vI + 2) % 3])
                        loop_AMVI.insert(airMesh.F[triI, (vI + 2) % 3])
                        E.row(tI + 2) <<== [tI + 2, 0]
                        bnd[0] = self.bnd[airMesh.F[triI, (vI + 2) % 3]]
                    }
                    break
                }
            }
        }
    }
    
     func getCornerAirLoop(_ corner_mesh: [Int],
                                 _ mergedPos: RVec2d,
                                 _ UV: inout Mat<Double>,
                                 _ E: inout Mat<Int>,
                                 _ bnd: inout Vec<Int>) -> Bool {
        
        assert(corner_mesh.count == 3)
        
        // convert mesh vertex index to air mesh vertex index
        var corner = [Int].init(repeating: 0, count: 3)
        for i in 0..<3 {
            if let value = meshVI2AirMesh[corner_mesh[i]] {
                corner[i] = value
            }
        }
        
        // get incident triangles
        var incTris = Set<Int>()
        var boundaryEdge: Pair<Int, Int> = Pair(0, 0)
        for vI in 0..<3 {
            var incTris_temp = [[Int]].init(repeating: [], count: 2)
            airMesh.isBoundaryVert(corner[vI], airMesh.vNeighbor[corner[vI]].first!, &incTris_temp[0], &boundaryEdge, false)
            incTris.formUnion(incTris_temp[0])
            
            airMesh.isBoundaryVert(corner[vI], airMesh.vNeighbor[corner[vI]].first!, &incTris_temp[1], &boundaryEdge, true)
            incTris.formUnion(incTris_temp[1])
            
            assert(incTris_temp[0].count + incTris_temp[1].count > 0)
        }
        
        var F_inc = Mati(incTris.count, 3)
        var fI: Int = 0
        for triI in incTris {
            F_inc.row(fI) <<== airMesh.F.row(triI)
            fI += 1
        }
        
        // compute outer loop and ensure no vertex duplication
        var loop = [Int]()
        boundary_loop(F_inc, &loop)
        var testDuplication: Set<Int> = .init(loop)
        
        if (testDuplication.count != loop.count) {
            assert(false, "vertex duplication found in loop!")
        }
        
        // eliminate merged vertices from loop
        var startI: Int = -1
        for vI in 0..<loop.count {
            if (loop[vI] == corner[0]) {
                if ((loop[(vI - 1 + loop.count) % loop.count] == corner[1]) &&
                    (loop[(vI - 2 + loop.count) % loop.count] == corner[2])) {
                    startI = vI
                    break
                } else {
                    assert(false, "corner not found on air mesh boundary loop")
                }
            }
        }
        assert(startI >= 0)
        var delI: Int = (startI - 2 + loop.count) % loop.count
        if (delI + 1 == loop.count) {
            loop.remove(at: loop.startIndex + delI)
            loop.remove(at: loop.startIndex)
            startI = 0
        } else {
            loop.remove(at: loop.startIndex + delI)
            loop.remove(at: loop.startIndex + delI)
        }
        if (startI < 3) {
            startI = 0
        } else {
            startI -= 2
        }
        
        // sort loop
        var loop_merged = [Int]()
        let beginI: Int = (startI - 1 + loop.count) % loop.count
        loop_merged.append(contentsOf: loop[(loop.startIndex + beginI)..<loop.endIndex])
        loop_merged.append(contentsOf: loop[loop.startIndex..<(loop.startIndex + beginI)])
        
        // construct air mesh data
        bnd.resize(3)
        for i in 0..<3 {
            assert(loop_merged[i] < self.bnd.count)
            bnd[i] = self.bnd[loop_merged[i]]
        }
        
        E.resize(loop_merged.count, 2)
        UV.resize(loop_merged.count, 2)
        for i in 0..<(loop_merged.count - 1) {
            UV.row(i) <<== airMesh.V.row(loop_merged[i])
            E.row(i) <<== [i, i + 1]
        }
        E.bottomRows(1) <<== [loop_merged.count - 1, 0]
        UV.bottomRows(1) <<== airMesh.V.row(loop_merged.last!)
        UV.row(1) <<== mergedPos
        
        // check loop boundary intersection
        for eJ in 2..<(E.rows - 1) {
            if (IglUtils.Test2DsegmentSegment(UV.row(E[0, 0]), UV.row(E[0, 1]),
                                              UV.row(E[eJ, 0]), UV.row(E[eJ, 1]))) {
                return false
            }
        }
        for eI in 1..<(E.rows - 2) {
            for eJ in (eI + 2)..<E.rows {
                if (IglUtils.Test2DsegmentSegment(UV.row(E[eI, 0]), UV.row(E[eI, 1]),
                                                  UV.row(E[eJ, 0]), UV.row(E[eJ, 1]))) {
                    return false
                }
            }
        }
        
        // ensure the loop is not totally inverted
        var rotAngle: Double = 0.0
        for eI in 0..<(E.rows - 1) {
            rotAngle += IglUtils.computeRotAngle(UV.row(E[eI, 1]) - UV.row(E[eI, 0]),
                                                 UV.row(E[eI + 1, 1]) - UV.row(E[eI + 1, 0]))
        }
        rotAngle += IglUtils.computeRotAngle(UV.row(E[E.rows - 1, 1]) - UV.row(E[E.rows - 1, 0]),
                                             UV.row(E[0, 1]) - UV.row(E[0, 0]))
        if (rotAngle > 0.0) {
            return true
        } else {
            return false
        }
    }
}
