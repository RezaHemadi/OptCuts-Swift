//
//  parameterize.swift
//  OptCuts
//
//  Created by Reza on 10/12/23.
//

import Foundation
import Matrix
import MetalKit
import GeometryProcessing

public func parameterize(device: MTLDevice, inputV: Matd, inputF: Mati, inputN: Matd, p_method: OptCuts.MethodType = .MT_OPTCUTS, p_lambda: Double?, p_upperBound: Double?, bijective: Bool = true, p_initCut: Int = 0, scaleUV: Bool) throws -> Data {
    start = DispatchTime.now()
    E_se_bestFeasible = Double.greatestFiniteMagnitude
    
    offlineMode = true
    headlessMode = true
    print("Headless mode")
    
    V = .init(inputV, inputV.rows, inputV.cols)
    N = .init(inputN, inputN.rows, inputN.cols)
    F = .init(inputF, inputF.rows, inputF.cols)
    
    vertAmt_input = V.rows
    var B = Vec<Int>()
    let isVertexManifold = is_vertex_manifold(F, &B)
    let isEdgeManifold = is_edge_manifold(F)
    let isManifold = isVertexManifold && isEdgeManifold
    if (!isManifold) {
        if !isVertexManifold {
            throw OptCutsError.nonManifoldVertex
        } else {
            throw OptCutsError.nonManifoldEdge
        }
    }
    
    // Set lambda
    lambda_init = 0.999
    if let lambda = p_lambda {
        lambda_init = lambda
        if lambda < 0.0 || lambda >= 1.0 {
            print("overrite invalid lambda \(lambda) to 0.999")
            lambda_init = 0.999
        }
    }
    
    methodType = p_method
    
    var startDS: String = ""
    switch methodType {
    case .MT_OPTCUTS_NODUAL:
        startDS = "OptCuts_noDual"
        
    case .MT_OPTCUTS:
        startDS = "OptCuts"
        
    case .MT_EBCUTS:
        startDS = "EBCuts"
        bijectiveParam = false
        
    case .MT_DISTMIN:
        lambda_init = 0.0
        startDS = "DistMin"
    }
    
    upperBound = p_upperBound ?? 4.1
    if (upperBound == 0.0) {
        var useScriptedBound: Bool = false
    } else {
        if (upperBound <= 4.0) {
            print("input b_d <= 4.0! use 4.1 instead.")
            upperBound = 4.1
        } else {
            print("use b_d = \(upperBound)")
        }
    }
    
    bijectiveParam = bijective
    print("bijectivity \(bijectiveParam ? "On" : "OFF")")
    
    initCutOption = p_initCut
    
    switch initCutOption {
    case 0:
        print("random 2-edge initial cut for closed surface")
        
    case 1:
        print("farthest 2-point initial cut for closed surface")
        
    default:
        print("input initial cut option invalid, use default")
        print("random 2-edge initial cut for closed surface")
        initCutOption = 0
    }
    
    var folderTail: String = ""
    
    if (UV.rows != 0) {
        // with input UV
        
    } else {
        // no input UV
        // * Harmonic map for initialization
        var bnd = Veci()
        boundary_loop(F, &bnd) // Find the open boundary
        if (bnd.count != 0) {
            // disk-topology
            
            // ODO: what if it has multiple boundaries or multi-components?
            
            // Map the boundary to a circle, preserving edge proportions
            var bnd_uv: Matd!
            IglUtils.map_vertices_to_circle(V, bnd, &bnd_uv)
            
            var UV_Tutte = Matd()
            
            // Harmonic map with uniform weights
            if (bnd.count == V.rows) {
                UV_Tutte.resize(V.rows, 2)
                for bndVI in 0..<bnd_uv.rows {
                    UV_Tutte.row(bnd[bndVI]) <<== bnd_uv.row(bndVI)
                }
            } else {
                var A: SparseMatrix<Double>!
                var M = SparseMatrix<Double>()
                IglUtils.computeUniformLaplacian(F, &A)
                do {
                    UV_Tutte = try harmonic(A, M, bnd, bnd_uv, 1)
                } catch {
                    print(error.localizedDescription)
                    fatalError("harmnic map generation failed!")
                }
            }
            
            triSoup.append(.init(V, F, UV_Tutte, Mati(), false))
        } else {
            // closed surface
            let genus: Int = 1 - euler_characteristic(V, F) / 2
            if (genus != 0) {
                print("Input surface genus = " + String(genus) + " or has multiple connected components!")
                
                var cuts: [[Int]] = [[]]
                cut_to_disk(F, &cuts)
                
                // record cohesive edge information,
                // transfer information format for cut_mesh
                var temp = TriMesh(V, F, Matd(), Mati(), false)
                var cutFlags = Mati(F.rows, 3)
                var cohEdgeRecord = Mati()
                cutFlags.setZero()
                for seamI in cuts {
                    for segI in 0..<(seamI.count - 1) {
                        var edge = Pair<Int, Int>(seamI[segI], seamI[segI + 1])
                        guard let tri = temp.edge2Tri[edge] else { fatalError() }
                        var i: Int = 0
                        for _ in 0..<3 {
                            if (temp.F[tri, i] == edge.first) {
                                cutFlags[tri, i] = 1
                                break
                            }
                            i += 1
                        }
                        
                        let cohERI: Int = cohEdgeRecord.rows
                        cohEdgeRecord.conservativeResize(cohERI + 1, 4)
                        cohEdgeRecord[cohERI, 0] = tri
                        cohEdgeRecord[cohERI, 1] = i
                        
                        edge.second = seamI[segI]
                        edge.first = seamI[segI + 1]
                        guard let tri = temp.edge2Tri[edge] else { fatalError() }
                        for i in 0..<3 {
                            if (temp.F[tri, i] == edge.first) {
                                cutFlags[tri, i] = 1
                                break
                            }
                        }
                        
                        cohEdgeRecord[cohERI, 2] = tri
                        cohEdgeRecord[cohERI, 3] = i
                    }
                }
                
                var Vcut = Matd()
                var Fcut = Mati()
                cut_mesh(temp.V_rest, temp.F, cutFlags, &Vcut, &Fcut)
                
                V = Vcut
                F = Fcut
                
                boundary_loop(F, &bnd) // Find the open boundary
                assert(bnd.count != 0)
                
                var bnd_uv: Matd!
                IglUtils.map_vertices_to_circle(V, bnd, &bnd_uv)
                
                var UV_Tutte: Matd!
                var A: SparseMatrix<Double>!
                var M = SparseMatrix<Double>()
                IglUtils.computeUniformLaplacian(F, &A)
                do {
                    UV_Tutte = try harmonic(A, M, bnd, bnd_uv, 1)
                } catch {
                    print(error.localizedDescription)
                    exit(EXIT_FAILURE)
                }
                
                var ptr = TriMesh(V, F, UV_Tutte, Mati(), false)
                ptr.buildCohEfromRecord(cohEdgeRecord)
                triSoup.append(ptr)
            } else {
                var temp: TriMesh = TriMesh(V, F, Matd(), Mati(), false)
                
                switch initCutOption {
                case 0:
                    temp.onePointCut()
                    rand1PinitCut = true
                    
                case 1:
                    temp.farthestPointCut()
                    
                default:
                    assert(false)
                }
                
                boundary_loop(temp.F, &bnd)
                assert(bnd.count != 0)
                
                var bnd_uv: Matd!
                IglUtils.map_vertices_to_circle(temp.V_rest, bnd, &bnd_uv)
                
                var A: SparseMatrix<Double>!
                var M = SparseMatrix<Double>()
                IglUtils.computeUniformLaplacian(temp.F, &A)
                
                var UV_Tutte: Matd!
                do {
                    UV_Tutte = try harmonic(A, M, bnd, bnd_uv, 1)
                } catch {
                    print(error.localizedDescription)
                }
                triSoup.append(.init(V, F, UV_Tutte, temp.F, false, temp.initSeamLen))
                
                // try initialize one-point cut with different vertices
                // until no inversion is detected
                var splitVI: Int = 0
                while (!triSoup.last!.checkInversion(true)) {
                    print("element inversion detected during UV initialization due to rounding errors, trying another vertex...")
                    
                    temp = TriMesh(V, F, Matd(), Mati(), false)
                    splitVI += 1
                    temp.onePointCut(splitVI)
                    
                    boundary_loop(temp.F, &bnd)
                    assert(bnd.count != 0)
                    
                    IglUtils.map_vertices_to_circle(temp.V_rest, bnd, &bnd_uv)
                    
                    IglUtils.computeUniformLaplacian(temp.F, &A)
                    
                    do {
                        UV_Tutte = try harmonic(A, M, bnd, bnd_uv, 1)
                    } catch {
                        print(error.localizedDescription)
                    }
                    
                    triSoup.removeLast()
                    triSoup.append(.init(V, F, UV_Tutte, temp.F, false, temp.initSeamLen))
                }
            }
        }
    }
    
    // * Our approach
    texScale = 10.0 / (triSoup[0].bbox.row(1) - triSoup[0].bbox.row(0)).maxCoeff()
    let newParam: UnsafeMutablePointer<Double> = .allocate(capacity: 1)
    newParam.initialize(to: 1.0 - lambda_init)
    energyParams.append(newParam)
    energyTerms.append(SymDirichletEnergy())
    
    optimizer = try Optimizer(triSoup[0], energyTerms, energyParams, 0, false, bijectiveParam && !rand1PinitCut) // for random one point initial cut, don't need air meshes in the beginning since it's impossible for a quad to itself
    
    optimizer.precompute()
    triSoup.append(optimizer.getResult())
    triSoup_backup = optimizer.getResult().clone()
    triSoup.append(optimizer.getData_findExtrema()) // for visualizing UV map for finding extrema
    
    if (lambda_init > 0.0) {
        // fracture mode
        fractureMode = true
    }
    
    // fatalError("1575") // if you want to specify seams
    if (headlessMode) {
        while (!outerLoopFinished) {
            preDrawFunc()
            postDrawFunc()
        }
    }
    
    let elapsed = DispatchTime.now().uptimeNanoseconds - start.uptimeNanoseconds
    let secs = Double(elapsed) / 1_000_000_000
    print("optimization took \(secs) seconds.")
    
    let triMesh = triSoup[channel_result]
    let data = triMesh.saveAsMesh(F0: F, scaleUV: scaleUV)
    cleanup()
    return data
}

public func parameterize(device: MTLDevice, inputURL: URL, p_method: OptCuts.MethodType = .MT_OPTCUTS, p_lambda: Double?, p_upperBound: Double?, bijective: Bool = true, p_initCut: Int = 0) throws -> Data {
    start = DispatchTime.now()
    E_se_bestFeasible = Double.greatestFiniteMagnitude
    
    offlineMode = true
    headlessMode = true
    print("Headless mode")
    
    // Load Mesh
    var V_Float = Mat<Float>()
    var N_Float = Mat<Float>()
    var F_UInt32 = Mat<UInt32>()
    let allocator = MTKMeshBufferAllocator(device: device)
    GeometryLoader.LoadMesh(allocator: allocator,
                            device: device,
                            url: inputURL,
                            vertices: &V_Float, normals: &N_Float, faces: &F_UInt32)
    V = V_Float.castToDouble()
    N = N_Float.castToDouble()
    F = F_UInt32.castToInt()
    
    vertAmt_input = V.rows
    var B = Vec<Int>()
    let isVertexManifold = is_vertex_manifold(F, &B)
    let isEdgeManifold = is_edge_manifold(F)
    let isManifold = isVertexManifold && isEdgeManifold
    if (!isManifold) {
        if !isVertexManifold {
            throw OptCutsError.nonManifoldVertex
        } else {
            throw OptCutsError.nonManifoldEdge
        }
    }
    
    // Set lambda
    lambda_init = 0.999
    if let lambda = p_lambda {
        lambda_init = lambda
        if lambda < 0.0 || lambda >= 1.0 {
            print("overrite invalid lambda \(lambda) to 0.999")
            lambda_init = 0.999
        }
    }
    
    methodType = p_method
    
    var startDS: String = ""
    switch methodType {
    case .MT_OPTCUTS_NODUAL:
        startDS = "OptCuts_noDual"
        
    case .MT_OPTCUTS:
        startDS = "OptCuts"
        
    case .MT_EBCUTS:
        startDS = "EBCuts"
        bijectiveParam = false
        
    case .MT_DISTMIN:
        lambda_init = 0.0
        startDS = "DistMin"
    }
    
    upperBound = p_upperBound ?? 4.1
    if (upperBound == 0.0) {
        var useScriptedBound: Bool = false
    } else {
        if (upperBound <= 4.0) {
            print("input b_d <= 4.0! use 4.1 instead.")
            upperBound = 4.1
        } else {
            print("use b_d = \(upperBound)")
        }
    }
    
    bijectiveParam = bijective
    print("bijectivity \(bijectiveParam ? "On" : "OFF")")
    
    initCutOption = p_initCut
    
    switch initCutOption {
    case 0:
        print("random 2-edge initial cut for closed surface")
        
    case 1:
        print("farthest 2-point initial cut for closed surface")
        
    default:
        print("input initial cut option invalid, use default")
        print("random 2-edge initial cut for closed surface")
        initCutOption = 0
    }
    
    var folderTail: String = ""
    
    if (UV.rows != 0) {
        // with input UV
        
    } else {
        // no input UV
        // * Harmonic map for initialization
        var bnd = Veci()
        boundary_loop(F, &bnd) // Find the open boundary
        if (bnd.count != 0) {
            // disk-topology
            
            // ODO: what if it has multiple boundaries or multi-components?
            
            // Map the boundary to a circle, preserving edge proportions
            var bnd_uv: Matd!
            IglUtils.map_vertices_to_circle(V, bnd, &bnd_uv)
            
            var UV_Tutte = Matd()
            
            // Harmonic map with uniform weights
            if (bnd.count == V.rows) {
                UV_Tutte.resize(V.rows, 2)
                for bndVI in 0..<bnd_uv.rows {
                    UV_Tutte.row(bnd[bndVI]) <<== bnd_uv.row(bndVI)
                }
            } else {
                var A: SparseMatrix<Double>!
                var M = SparseMatrix<Double>()
                IglUtils.computeUniformLaplacian(F, &A)
                do {
                    UV_Tutte = try harmonic(A, M, bnd, bnd_uv, 1)
                } catch {
                    print(error.localizedDescription)
                    fatalError("harmnic map generation failed!")
                }
            }
            
            triSoup.append(.init(V, F, UV_Tutte, Mati(), false))
        } else {
            // closed surface
            let genus: Int = 1 - euler_characteristic(V, F) / 2
            if (genus != 0) {
                print("Input surface genus = " + String(genus) + " or has multiple connected components!")
                
                var cuts: [[Int]] = [[]]
                cut_to_disk(F, &cuts)
                
                // record cohesive edge information,
                // transfer information format for cut_mesh
                var temp = TriMesh(V, F, Matd(), Mati(), false)
                var cutFlags = Mati(F.rows, 3)
                var cohEdgeRecord = Mati()
                cutFlags.setZero()
                for seamI in cuts {
                    for segI in 0..<(seamI.count - 1) {
                        var edge = Pair<Int, Int>(seamI[segI], seamI[segI + 1])
                        guard let tri = temp.edge2Tri[edge] else { fatalError() }
                        var i: Int = 0
                        for _ in 0..<3 {
                            if (temp.F[tri, i] == edge.first) {
                                cutFlags[tri, i] = 1
                                break
                            }
                            i += 1
                        }
                        
                        let cohERI: Int = cohEdgeRecord.rows
                        cohEdgeRecord.conservativeResize(cohERI + 1, 4)
                        cohEdgeRecord[cohERI, 0] = tri
                        cohEdgeRecord[cohERI, 1] = i
                        
                        edge.second = seamI[segI]
                        edge.first = seamI[segI + 1]
                        guard let tri = temp.edge2Tri[edge] else { fatalError() }
                        for i in 0..<3 {
                            if (temp.F[tri, i] == edge.first) {
                                cutFlags[tri, i] = 1
                                break
                            }
                        }
                        
                        cohEdgeRecord[cohERI, 2] = tri
                        cohEdgeRecord[cohERI, 3] = i
                    }
                }
                
                var Vcut = Matd()
                var Fcut = Mati()
                cut_mesh(temp.V_rest, temp.F, cutFlags, &Vcut, &Fcut)
                do {
                    try GeometryLoader.WriteToFile(outputFolderPath.appending(component: "_disk.obj"),
                                                   allocator: allocator,
                                                   vertexData: Vcut.castToFloat(),
                                                   faceData: Fcut.castToUInt32())
                } catch {
                    print(error.localizedDescription)
                }
                
                V = Vcut
                F = Fcut
                
                boundary_loop(F, &bnd) // Find the open boundary
                assert(bnd.count != 0)
                
                var bnd_uv: Matd!
                IglUtils.map_vertices_to_circle(V, bnd, &bnd_uv)
                
                var UV_Tutte: Matd!
                var A: SparseMatrix<Double>!
                var M = SparseMatrix<Double>()
                IglUtils.computeUniformLaplacian(F, &A)
                do {
                    UV_Tutte = try harmonic(A, M, bnd, bnd_uv, 1)
                } catch {
                    print(error.localizedDescription)
                    exit(EXIT_FAILURE)
                }
                
                var ptr = TriMesh(V, F, UV_Tutte, Mati(), false)
                ptr.buildCohEfromRecord(cohEdgeRecord)
                triSoup.append(ptr)
                outputFolderPath.append(component: inputURL.lastPathComponent + "_HighGenus" + String(lambda_init) + "_" + startDS + folderTail)
            } else {
                var temp: TriMesh = TriMesh(V, F, Matd(), Mati(), false)
                
                switch initCutOption {
                case 0:
                    temp.onePointCut()
                    rand1PinitCut = true
                    
                case 1:
                    temp.farthestPointCut()
                    
                default:
                    assert(false)
                }
                
                boundary_loop(temp.F, &bnd)
                assert(bnd.count != 0)
                
                var bnd_uv: Matd!
                IglUtils.map_vertices_to_circle(temp.V_rest, bnd, &bnd_uv)
                
                var A: SparseMatrix<Double>!
                var M = SparseMatrix<Double>()
                IglUtils.computeUniformLaplacian(temp.F, &A)
                
                var UV_Tutte: Matd!
                do {
                    UV_Tutte = try harmonic(A, M, bnd, bnd_uv, 1)
                } catch {
                    print(error.localizedDescription)
                }
                triSoup.append(.init(V, F, UV_Tutte, temp.F, false, temp.initSeamLen))
                
                // try initialize one-point cut with different vertices
                // until no inversion is detected
                var splitVI: Int = 0
                while (!triSoup.last!.checkInversion(true)) {
                    print("element inversion detected during UV initialization due to rounding errors, trying another vertex...")
                    
                    temp = TriMesh(V, F, Matd(), Mati(), false)
                    splitVI += 1
                    temp.onePointCut(splitVI)
                    
                    boundary_loop(temp.F, &bnd)
                    assert(bnd.count != 0)
                    
                    IglUtils.map_vertices_to_circle(temp.V_rest, bnd, &bnd_uv)
                    
                    IglUtils.computeUniformLaplacian(temp.F, &A)
                    
                    do {
                        UV_Tutte = try harmonic(A, M, bnd, bnd_uv, 1)
                    } catch {
                        print(error.localizedDescription)
                    }
                    
                    triSoup.removeLast()
                    triSoup.append(.init(V, F, UV_Tutte, temp.F, false, temp.initSeamLen))
                }
                
                outputFolderPath.append(component: inputURL.lastPathComponent +
                "_Tutte_" + String(lambda_init) + "_" + startDS + folderTail)
            }
        }
    }
    
    // * Our approach
    texScale = 10.0 / (triSoup[0].bbox.row(1) - triSoup[0].bbox.row(0)).maxCoeff()
    let newParam: UnsafeMutablePointer<Double> = .allocate(capacity: 1)
    newParam.initialize(to: 1.0 - lambda_init)
    energyParams.append(newParam)
    energyTerms.append(SymDirichletEnergy())
    
    optimizer = try Optimizer(triSoup[0], energyTerms, energyParams, 0, false, bijectiveParam && !rand1PinitCut) // for random one point initial cut, don't need air meshes in the beginning since it's impossible for a quad to itself
    
    optimizer.precompute()
    triSoup.append(optimizer.getResult())
    triSoup_backup = optimizer.getResult().clone()
    triSoup.append(optimizer.getData_findExtrema()) // for visualizing UV map for finding extrema
    
    if (lambda_init > 0.0) {
        // fracture mode
        fractureMode = true
    }
    
    // fatalError("1575") // if you want to specify seams
    if (headlessMode) {
        while (!outerLoopFinished) {
            preDrawFunc()
            postDrawFunc()
        }
    }
    
    let elapsed = DispatchTime.now().uptimeNanoseconds - start.uptimeNanoseconds
    let secs = Double(elapsed) / 1_000_000_000
    print("optimization took \(secs) seconds.")
    
    let triMesh = triSoup[channel_result]
    let data = triMesh.saveAsMesh(F0: F, scaleUV: true)
    cleanup()
    return data
}
