//
//  Optimizer.swift
//  OptCuts
//
//  Created by Reza on 6/20/23.
//

import Foundation
import Matrix
import GeometryProcessing

private let HessianCoefMaxCount: Int = 400_000

/// A class for solving an optimization problem
 class Optimizer {
    // MARK: - Properties
     private var I: Vec<Int> = .init()
     private var J: Vec<Int> = .init()
     private var V: Vec<Double> = .init()
    
    // these data are references
    private let data0: TriMesh // Initial guess
    private let energyTerms: [Energy]  // E_0, E_1, E_2, ...
    private let energyParams: [UnsafeMutablePointer<Double>] // a_0, a_1, a_2, ...
    // E = \Sigma_i a_i E_i
    
    // Owned data
    private var useDense: Bool
    private var propogateFracture: Int
    private var fractureInitiated: Bool = false
    private var allowEDecRelTol: Bool
    private let mute: Bool
    private let paradisoThreadAmt: Bool
    private var needRefactorize: Bool
    private var globalIterNum: Int
    private var topoIter: Int
    private var relGL2Tol: Double
    private var energyParamSum: Double
    private var result: TriMesh! // intermediate results of each iteration
    private var data_findExtrema: TriMesh! // intermediate results for deciding  the cuts in each topology step
    private var scaffolding: Bool // whether to enable bijectivity parameterization
    private let w_scaf: Double
    private var scaffold: Scaffold = .init() // air meshes to enforce bijectivity
    private var I_mtr = Vec<Int>() // triplet representation
    private var J_mtr = Vec<Int>() // triplet representation
    private var V_mtr = Vec<Double>()
    private var Hessian: Mat<Double> = .init() // when using dense representation
    // cholesky solver for solving the linear system for search directions
    private var linSysSolver: LinSysSolver<Vec<Int>, Vec<Double>>!
    //private let denseSolver: Eigen::LDLT<Eigen::MatrixXd>
    private var gradient = Vec<Double>() // energy gradient computed in each iteration
    private var searchDir = Vec<Double>() // search direction computed in each iteration
    private var lastEnergyVal: Double = .zero
    private var lastEDec: Double!
    private var targetGres: Double!
    private var gradient_ET: [Vec<Double>]
    private var gradient_scaffold = Vec<Double>()
    private var energyVal_ET: [Double]
    private var energyVal_scaffold: Double = .zero
    
    private let UV_bnds_scaffold: Mat<Double>
    private let E_scaffold: Mat<Int>
    private let bnd_scaffold: Vec<Int>
    private var vNeighbor_withScaf: [Set<Int>] = []
    private var fixedV_withScaf: Set<Int> = .init()
    private var denseSolver: LDLT!
    
    private var lastPropogate: Bool = false
    
    // std::ostringstream buffer_energyValPerIter
    // std::ostringstream buffer_gradientPerIter
    // std::ofstream file_energyValPerIter
    // std::ofstream file_gradientPerIter
    
    // MARK: - Initialization
     init(_ p_data0: TriMesh,
                _ p_energyTerms: [Energy],
                _ p_energyParams: [UnsafeMutablePointer<Double>],
                _ p_propogateFracture: Int = 1,
                _ p_mute: Bool = false,
                _ p_scaffolding: Bool = false,
                _ UV_bnds: Mat<Double> = Matd(),
                _ E: Mat<Int> = Mati(),
                _ bnd: Vec<Int> = Veci(),
                _ p_useDense: Bool = false) {
        
        data0 = p_data0
        energyTerms = p_energyTerms
        energyParams = p_energyParams
        
        assert(energyTerms.count == energyParams.count)
        
        useDense = p_useDense
        
        energyParamSum = 0.0
        for ePI in energyParams {
            energyParamSum += ePI.pointee
        }
        
        gradient_ET = .init(repeating: Vec<Double>(), count: energyTerms.count)
        energyVal_ET = .init(repeating: 0.0, count: energyTerms.count)
        
        allowEDecRelTol = true
        propogateFracture = p_propogateFracture
        mute = p_mute
        
        if (!data0.checkInversion()) {
            exit(-1)
        }
        
        globalIterNum = 0
        relGL2Tol = 1.0e-12
        topoIter = 0
        
        needRefactorize = false
        for energyTermI in energyTerms {
            if (energyTermI.getNeedRefactorize()) {
                needRefactorize = true
                break
            }
        }
        
        paradisoThreadAmt = true
        
        scaffolding = p_scaffolding
        UV_bnds_scaffold = .init(UV_bnds, UV_bnds.rows, UV_bnds.cols)
        E_scaffold = .init(E, E.rows, E.cols)
        bnd_scaffold = .init(bnd, bnd.rows, bnd.cols)
        w_scaf = energyParams[0].pointee * 0.01
         
         if !p_useDense {
             linSysSolver = EigenLibSolver<Vec<Int>, Vec<Double>>()
             I.reserveCapacity(HessianCoefMaxCount)
             J.reserveCapacity(HessianCoefMaxCount)
             V.reserveCapacity(HessianCoefMaxCount)
         }
    }
    
    // MARK: - Methods
    
    /// Precompute preconditioning matrix and factorize for fast solve, prepare initial guess
     func precompute() {
        result = data0.clone()
        if scaffolding {
            scaffold = Scaffold(result, UV_bnds_scaffold, E_scaffold, bnd_scaffold)
            result.scaffold = scaffold
            scaffold.mergeVNeighbor(result.vNeighbor, &vNeighbor_withScaf)
            scaffold.mergeFixedV(result.fixedVert, &fixedV_withScaf)
        }
        
        computeHessian(result, scaffold)
        
        if (useDense) {
            if (!needRefactorize) {
                //denseSolver = LDLT(Hessian)
            }
        } else {
            linSysSolver.set_type(4, -2)
            linSysSolver.set_pattern(scaffolding ? vNeighbor_withScaf : result.vNeighbor,
                                     scaffolding ? fixedV_withScaf : result.fixedVert)
            linSysSolver.update_a(I_mtr, J_mtr, V_mtr)
            linSysSolver.analyze_pattern()
            if (!needRefactorize) {
                do {
                    linSysSolver.factorize()
                } catch {
                    print(error.localizedDescription)
                    exit(-1)
                }
            }
        }
        
        lastEDec = 0.0
        data_findExtrema = data0.clone()
        updateTargetGRes()
        computeEnergyVal(result, scaffold, &lastEnergyVal)
    }
    
    /// solve the optimization problem that minimizes E using a hill-climbing method,
    /// the final result will be in result
    @discardableResult
     func solve(_ maxIter: Int = 100) -> Int {
        for _ in 0..<maxIter {
            computeGradient(result, scaffold, &gradient)
            let sqn_g: Double = gradient.squaredNorm()
            if (!mute) {
                print("||gradient||^2 = \(sqn_g), targetGRes = \(targetGres!)")
            }
            if (sqn_g < targetGres) {
                // converged
                lastEDec = 0.0
                globalIterNum += 1
                return 1
            } else {
                if (solve_oneStep()) {
                    globalIterNum += 1
                    return 1
                }
            }
            globalIterNum += 1
            
            if (propogateFracture > 0) {
                if (!createFracture(lastEDec, propogateFracture)) {
                    // always perform the one decreasing E_w more
                    if (scaffolding) {
                        scaffold = Scaffold(result, UV_bnds_scaffold, E_scaffold, bnd_scaffold)
                        result.scaffold = scaffold
                        scaffold.mergeVNeighbor(result.vNeighbor, &vNeighbor_withScaf)
                        scaffold.mergeFixedV(result.fixedVert, &fixedV_withScaf)
                    }
                    
                    if (lastPropogate) {
                        lastPropogate = false
                        return 2 // for saving screenshots
                    }
                } else {
                    lastPropogate = true
                }
            } else {
                if (scaffolding) {
                    scaffold = Scaffold(result, UV_bnds_scaffold, E_scaffold, bnd_scaffold)
                    result.scaffold = scaffold
                    scaffold.mergeVNeighbor(result.vNeighbor, &vNeighbor_withScaf)
                    scaffold.mergeFixedV(result.fixedVert, &fixedV_withScaf)
                }
            }
        }
        return 0
    }
    
     func updatePrecondMtrAndFactorize() {
        if (needRefactorize) {
            // don't need to call this function
            return
        }
        
        if (!mute) {
            print("precompute proxy/Hessian matrix and factorize...")
        }
        computeHessian(result, scaffold)
        
        if (useDense) {
            denseSolver = LDLT(Hessian)
        } else {
            linSysSolver.update_a(I_mtr, J_mtr, V_mtr)
            linSysSolver.factorize()
        }
    }
    
     func updateEnergyData(_ updateEVal: Bool = true,
                                 _ updateGradient: Bool = true,
                                 _ updateHessian: Bool = true) {
        
        energyParamSum = 0.0
        for ePI in energyParams {
            energyParamSum += ePI.pointee
        }
        updateTargetGRes()
        
        if (updateEVal) {
            // compute energy and output
            computeEnergyVal(result, scaffold, &lastEnergyVal)
        }
        
        if (updateGradient) {
            // compute gradient and output
            computeGradient(result, scaffold, &gradient)
            if (gradient.squaredNorm() < targetGres) {
                print("||g||^2 = \(gradient.squaredNorm()) << after fracture initialization!")
            }
        }
        
        if (updateHessian) {
            // for the changing hessian
            if (!mute) {
                print("recompute proxy/Hessian matrix and factorize...")
            }
            computeHessian(result, scaffold)
            
            if (useDense) {
                if (!needRefactorize) {
                    denseSolver = LDLT(Hessian)
                }
            } else {
                linSysSolver.set_pattern(scaffolding ? vNeighbor_withScaf : result.vNeighbor,
                                         scaffolding ? fixedV_withScaf : result.fixedVert)
                linSysSolver.update_a(I_mtr, J_mtr, V_mtr)
                linSysSolver.analyze_pattern()
                if (!needRefactorize) {
                    linSysSolver.factorize()
                }
            }
        }
    }
    
     func createFracture(_ stressThres: Double,
                               _ propType: Int,
                               _ allowPropogate: Bool = true,
                               _ allowInSplit: Bool = false) -> Bool {
        var allowPropogate = allowPropogate
        
        if (propType == 0) {
            topoIter += 1
        }
        var changed: Bool = false
        var isMerge: Bool = false
        
        switch (methodType) {
        case .MT_OPTCUTS_NODUAL, .MT_OPTCUTS:
            data_findExtrema = result.clone()
            switch (propType) {
            case 0: // initiation
                changed = result.splitOrMerge(1.0 - energyParams[0].pointee, stressThres, false, allowInSplit, &isMerge)
                break
                
            case 1: // propogate split
                changed = result.splitEdge(1.0 - energyParams[0].pointee, stressThres, true, allowInSplit)
                break
                
            case 2: // propogate merge
                changed  = result.mergeEdge(1.0 - energyParams[0].pointee, stressThres, true)
                isMerge = true
                break
                
            default:
                assert(false)
                break
            }
            break
            
        case .MT_EBCUTS:
            result.geomImgCut(&data_findExtrema)
            allowPropogate = false
            changed = true
            break
            
        default:
            assert(false, "Fracture forbidden for current method type!")
            break
        }
        
        if (changed) {
            if (scaffolding) {
                scaffold = Scaffold(result, UV_bnds_scaffold, E_scaffold, bnd_scaffold)
                result.scaffold = scaffold
                scaffold.mergeVNeighbor(result.vNeighbor, &vNeighbor_withScaf)
                scaffold.mergeFixedV(result.fixedVert, &fixedV_withScaf)
            }
            
            updateEnergyData(true, false, true)
            fractureInitiated = true
            if ((!mute) && (propType == 0)) {
                writeEnergyValToFile(false)
            }
            
            if (allowPropogate && (propType == 0)) {
                propogateFracture = 1 + (isMerge ? 1 : 0)
            }
        }
        
        return changed
    }
    
    @discardableResult
     func createFracture(_ opType: Int,
                               _ path: [Int],
                               _ newVertPos: Mat<Double>,
                               _ allowPropogate: Bool) -> Bool {
        
        assert(methodType == .MT_OPTCUTS)
        
        topoIter += 1
        
        var isMerge: Bool = false
        data_findExtrema = result.clone()
        
        switch (opType) {
        case 0: // boundary split
            print("boundary split without querying again")
            result.splitEdgeOnBoundary(Pair(path[0], path[1]), newVertPos)
            // TODO: process fractail here!
            result.updateFeatures()
            break
            
        case 1: // interior split
            print("Interior split without querying again")
            result.cutPath(path, true, 1, newVertPos)
            result.fracTail.insert(path[0])
            result.fracTail.insert(path[2])
            result.curInteriorFracTails.first = path[0]
            result.curInteriorFracTails.second = path[2]
            result.curFracTail = -1
            break
            
        case 2: // merge
            print("corner edge merged without querying again")
            result.mergeBoundaryEdges(Pair(path[0], path[1]),
                                      Pair(path[1], path[2]), newVertPos.row(0))
            result.computeFeatures() // TODO: only update locally
            isMerge = true
            break
            
        default:
            assert(false)
            break
        }
        
        if (scaffolding) {
            scaffold = Scaffold(result, UV_bnds_scaffold, E_scaffold, bnd_scaffold)
            result.scaffold = scaffold
            scaffold.mergeVNeighbor(result.vNeighbor, &vNeighbor_withScaf)
            scaffold.mergeFixedV(result.fixedVert, &fixedV_withScaf)
        }
        
        updateEnergyData(true, false, true)
        fractureInitiated = true
        if (!mute) {
            writeEnergyValToFile(false)
        }
        
        if (allowPropogate) {
            propogateFracture = 1 + (isMerge ? 1 : 0)
        }
        
        return true
    }
    
     func setConfig(_ config: TriMesh,
                   _ iterNum: Int,
                   _ p_topoIter: Int) {
        topoIter = p_topoIter
        globalIterNum = iterNum
        result = config.clone()
        if (scaffolding) {
            scaffold = Scaffold(result, UV_bnds_scaffold, E_scaffold, bnd_scaffold)
            result.scaffold = scaffold
            scaffold.mergeVNeighbor(result.vNeighbor, &vNeighbor_withScaf)
            scaffold.mergeFixedV(result.fixedVert, &fixedV_withScaf)
        }
        
        updateEnergyData()
    }
    
     func setPropogateFracture(_ p_prop: Bool) {
        propogateFracture = p_prop ? 1 : 0
    }
    
     func setScaffolding(_ p_scaffolding: Bool) {
        scaffolding = p_scaffolding
        if (scaffolding) {
            scaffold = Scaffold(result, UV_bnds_scaffold, E_scaffold, bnd_scaffold)
            result.scaffold = scaffold
            scaffold.mergeVNeighbor(result.vNeighbor, &vNeighbor_withScaf)
            scaffold.mergeFixedV(result.fixedVert, &fixedV_withScaf)
        }
    }
    
     func setUseDense(_ p_useDense: Bool = true) {
        useDense = p_useDense
    }
    
     func computeLastEnergyVal() {
        computeEnergyVal(result, scaffold, &lastEnergyVal)
    }
    
     func getGradientVisual(_ arrowVec: inout Mat<Double>) {
        assert(result.V.rows * 2 == gradient.count)
        
        arrowVec.resize(result.V.rows, result.V.cols)
        for vI in 0..<result.V.rows {
            arrowVec[vI, 0] = gradient[vI * 2]
            arrowVec[vI, 1] = gradient[vI * 2 + 1]
            arrowVec.row(vI).normalize()
        }
        arrowVec *= avg_edge_length(result.V, result.F)
    }
    
    // Note: return by reference in original c++ code
     func getResult() -> TriMesh {
        return result
    }
    
    // Note: return by reference in original c++ code
     func getScaffold() -> Scaffold {
        return scaffold
    }
    
    // Note: return by reference in original c++ code
     func getAirMesh() -> TriMesh {
        return scaffold.airMesh
    }
    
     func isScaffolding() -> Bool {
        return scaffolding
    }
    
    // Note: return by reference in original c++ code
     func getData_findExtrema() -> TriMesh {
        return data_findExtrema
    }
    
     func getIterNum() -> Int {
        return globalIterNum
    }
    
     func getTopoIter() -> Int {
        return topoIter
    }
    
     func setRelGL2Tol(_ p_relTol: Double) {
        assert(p_relTol > 0.0)
        relGL2Tol = p_relTol
        updateTargetGRes()
    }
    
     func setAllowEDecRelTol(_ p_allowEDecRelTol: Bool) {
        allowEDecRelTol = p_allowEDecRelTol
    }
    
     func flushEnergyFileOutput() {
        //fatalError("to be implemented")
    }
    
     func flushGradFileOutput() {
        //fatalError("to be implemented")
    }
    
     func clearEnergyFileOutputBuffer() {
        fatalError("to be implemented")
    }
    
     func clearGradFileOutputBuffer() {
        fatalError("to be implemented")
    }
    
    // MARK: - Helper Methods
    
    // solve for new configuration in the next iteration
    // Note: must compute current gradient first
     func solve_oneStep() -> Bool {
        if (needRefactorize) {
            // for the changing hessian
            if (!mute) {
                print("recompute proxy/Hessian matrix...")
            }
            if (!fractureInitiated) {
                computeHessian(result, scaffold)
            }
            if (!mute) {
                print("factorizing proxy/Hessian matrix...")
            }
            
            if (!fractureInitiated) {
                if (!useDense) {
                    if (scaffolding) {
                        linSysSolver.set_pattern(scaffolding ? vNeighbor_withScaf : result.vNeighbor,
                                                 scaffolding ? fixedV_withScaf : result.fixedVert)
                        linSysSolver.update_a(I_mtr, J_mtr, V_mtr)
                        linSysSolver.analyze_pattern()
                    } else {
                        linSysSolver.update_a(I_mtr, J_mtr, V_mtr)
                    }
                }
            }
            do {
                if (useDense) {
                    denseSolver = LDLT(Hessian)
                } else {
                    linSysSolver.factorize()
                }
            } catch {
                print(error.localizedDescription)
                exit(-1)
            }
        }
        
        var minusG: Vec<Double> = -gradient
        if (useDense) {
            denseSolver.solve(minusG, &searchDir)
        } else {
            linSysSolver.solve(&minusG, &searchDir)
        }
        
        fractureInitiated = false
        
        let stopped: Bool = lineSearch()
        
        return stopped
    }
    
     func lineSearch() -> Bool {
        var stopped: Bool = false
        var stepSize: Double = 1.0
        initStepSize(result, &stepSize)
        stepSize *= 0.99 // producing degenerated element is not allowed
        if (!mute) {
            print("stepSize: \(stepSize)")
        }
        
        var lastEnergyVal_scaffold: Double = 0.0
        let resultV0: Matd = .init(result.V, result.V.rows, result.V.cols)
        var scaffoldV0: Matd = .init()
        if (scaffolding) {
            scaffoldV0 = .init(scaffold.airMesh.V, scaffold.airMesh.V.rows, scaffold.airMesh.V.cols)
            computeEnergyVal(result, scaffold, &lastEnergyVal) // this update is necessary since scaffold changes
            lastEnergyVal_scaffold = energyVal_scaffold
        }
        stepForward(resultV0, scaffoldV0, &result, scaffold, stepSize)
        var testingE: Double = 0.0
        computeEnergyVal(result, scaffold, &testingE)
        
        while (testingE > lastEnergyVal) { // ensure energy decrese
            stepSize /= 2.0
            if (stepSize == 0.0) {
                stopped = true
                break
            }
            
            stepForward(resultV0, scaffoldV0, &result, scaffold, stepSize)
            computeEnergyVal(result, scaffold, &testingE)
        }
        if (!mute) {
            //print("\(stepSize) \"(armijo)\"")
        }
        /*
        while (!result.checkInversion() || (scaffolding && !scaffold.airMesh.checkInversion() )) {
            assert(false, "element inversion after armijo shouldn't happen!")
            
            stepSize /= 2.0
            if (stepSize == 0.0) {
                assert(false, "line search failed!")
                stopped = true
                break
            }
            
            stepForward(resultV0, scaffoldV0, &result, scaffold, stepSize)
            computeEnergyVal(result, scaffold, &testingE)
        }*/
        
        lastEDec = lastEnergyVal - testingE
        if (scaffolding) {
            lastEDec += (-lastEnergyVal_scaffold + energyVal_scaffold)
        }
        if (allowEDecRelTol && (lastEDec / lastEnergyVal < 1.0e-6 * stepSize) && (stepSize > 1.0e-3)) {
            // avoid stapping in hard situations
            stopped = true
        }
        lastEnergyVal = testingE
        
        if (!mute) {
            print("stepLen = \((stepSize * searchDir).squaredNorm())")
            print("E_cur_smooth = \(testingE - energyVal_scaffold)")
            
            if (!stopped) {
                writeEnergyValToFile(false)
            }
        }
        
        return stopped
    }
    
     func stepForward(_ dataV0: Mat<Double>,
                     _ scaffoldV0: Mat<Double>,
                     _ data: inout TriMesh,
                     _ scaffoldData: Scaffold,
                     _ stepSize: Double) {
        
        assert(dataV0.rows == data.V.rows)
        if (scaffolding) {
            assert(data.V.rows + scaffoldData.airMesh.V.rows - scaffoldData.bnd.count == searchDir.count / 2)
        } else {
            assert(data.V.rows * 2 == searchDir.count)
        }
        //assert(data.V.rows == result.V.rows)
        
        /*
        for vI in 0..<data.V.rows {
            data.V[vI, 0] = dataV0[vI, 0] + stepSize * searchDir[vI * 2]
            data.V[vI, 1] = dataV0[vI, 1] + stepSize * searchDir[vI * 2 + 1]
        }*/
         withUnsafeMutablePointer(to: &data.V) { vPtr in
             DispatchQueue.concurrentPerform(iterations: vPtr.pointee.rows) { vI in
                 vPtr.pointee[vI, 0] = dataV0[vI, 0] + stepSize * searchDir[vI * 2]
                 vPtr.pointee[vI, 1] = dataV0[vI, 1] + stepSize * searchDir[vI * 2 + 1]
             }
         }
         
        if (scaffolding) {
            scaffoldData.stepForward(scaffoldV0, searchDir, stepSize)
        }
    }
    
     func updateTargetGRes() {
        targetGres = energyParamSum * Double(data0.V_rest.rows - data0.fixedVert.count) / Double(data0.V_rest.rows) * relGL2Tol
    }
    
     func computeEnergyVal(_ data: TriMesh,
                          _ scaffoldData: Scaffold,
                          _ energyVal: inout Double,
                          _ excludeScaffold: Bool = false) {
        
        energyTerms[0].computeEnergyVal(data, &energyVal_ET[0])
        energyVal = energyParams[0].pointee * energyVal_ET[0]
        for eI in 1..<energyTerms.count {
            energyTerms[eI].computeEnergyVal(data, &energyVal_ET[eI])
            energyVal += energyParams[eI].pointee * energyVal_ET[eI]
        }
        
        if (scaffolding && !excludeScaffold) {
            let SD = SymDirichletEnergy()
            SD.computeEnergyVal(scaffoldData.airMesh, &energyVal_scaffold, true)
            energyVal_scaffold *= (w_scaf / Double(scaffold.airMesh.F.rows))
            energyVal += energyVal_scaffold
        } else {
            energyVal_scaffold = 0.0
        }
    }
    
    func computeGradient(_ data: TriMesh,
                         _ scaffoldData: Scaffold,
                         _ gradient: inout Vec<Double>,
                         _ excludeScaffold: Bool = false) {
        
        energyTerms[0].computeGradient(data, &gradient_ET[0])
        gradient = energyParams[0].pointee * gradient_ET[0]
        for eI in 1..<energyTerms.count {
            energyTerms[eI].computeGradient(data, &gradient_ET[eI])
            gradient += energyParams[eI].pointee * gradient_ET[eI]
        }
        
        if (scaffolding) {
            let SD = SymDirichletEnergy()
            SD.computeGradient(scaffoldData.airMesh, &gradient_scaffold, true)
            scaffoldData.augmentGradient(&gradient, gradient_scaffold, (excludeScaffold ? 0.0 : (w_scaf / Double(scaffold.airMesh.F.rows))))
        }
    }
    
     func computeHessian(_ data: TriMesh,
                        _ scaffoldData: Scaffold,
                               _ loud: Bool = false) {
        
        if (useDense) {
            energyTerms[0].computeHessian(data, &Hessian)
            Hessian *= energyParams[0].pointee
            
            for eI in 1..<energyTerms.count {
                var HessianI = Matd()
                energyTerms[eI].computeHessian(data, &HessianI)
                Hessian += (energyParams[eI].pointee * HessianI)
            }
            
            if (scaffolding) {
                let SD = SymDirichletEnergy()
                var Hessian_scaf = Matd()
                SD.computeHessian(scaffoldData.airMesh, &Hessian_scaf, true)
                scaffoldData.augmentProxyMatrix(&Hessian, Hessian_scaf, w_scaf / Double(scaffold.airMesh.F.rows))
                
            }
        } else {
            I_mtr.resize(0)
            J_mtr.resize(0)
            V_mtr.resize(0)
            for eI in 0..<energyTerms.count {
                I.resize(0)
                J.resize(0)
                V.resize(0)
                /*
                var I = Vec<Int>()
                I.reserveCapacity(HessianCoefMaxCount)
                var J = Vec<Int>()
                J.reserveCapacity(HessianCoefMaxCount)
                var V = Vec<Double>()
                V.reserveCapacity(HessianCoefMaxCount)*/
                energyTerms[eI].computeHessian(data, &V, &I, &J)
                V *= energyParams[eI].pointee
                /*
                I_mtr.conservativeResize(I_mtr.count + I.count)
                I_mtr.bottomRows(I.count) <<== I
                J_mtr.conservativeResize(J_mtr.count + J.count)
                J_mtr.bottomRows(J.count) <<== J
                V_mtr.conservativeResize(V_mtr.count + V.count)
                V_mtr.bottomRows(V.count) <<== V*/
                
                I_mtr = .init(I, I.rows, I.cols)
                J_mtr = .init(J, J.rows, J.cols)
                V_mtr = .init(V, V.rows, V.cols)
            }
            if (scaffolding) {
                let SD = SymDirichletEnergy()
                I.resize(0)
                J.resize(0)
                V.resize(0)
                /*
                var I = Vec<Int>()
                I.reserveCapacity(HessianCoefMaxCount)
                var J = Vec<Int>()
                J.reserveCapacity(HessianCoefMaxCount)
                var V = Vec<Double>()
                V.reserveCapacity(HessianCoefMaxCount)*/
                SD.computeHessian(scaffoldData.airMesh, &V, &I, &J, true)
                scaffoldData.augmentProxyMatrix(&I_mtr, &J_mtr, &V_mtr, I, J, V, w_scaf / Double(scaffold.airMesh.F.rows))
            }
        }
    }
    
     func initStepSize(_ data: TriMesh,
                             _ stepSize: inout Double) {
        
        for eI in 0..<energyTerms.count {
            energyTerms[eI].initStepSize(data, searchDir, &stepSize)
        }
        if (scaffolding) {
            var searchDir_scaffold = Vec<Double>()
            scaffold.wholeSearchDir2airMesh(searchDir, &searchDir_scaffold)
            let SD = SymDirichletEnergy()
            SD.initStepSize(scaffold.airMesh, searchDir_scaffold, &stepSize)
        }
    }
    
     func writeEnergyValToFile(_ flush: Bool) {
        //fatalError("to be implemented")
    }
    
     func writeGradL2NormToFile(flush: Bool) {
        //fatalError("to be implemented")
    }
    
     func getLastEnergyVal(_ exculdeScaffold: Bool = false) -> Double {
        return ((exculdeScaffold && scaffolding) ? (lastEnergyVal - energyVal_scaffold) : lastEnergyVal)
    }
}
