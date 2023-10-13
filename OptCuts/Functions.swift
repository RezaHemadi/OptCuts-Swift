//
//  Functions.swift
//  OptCuts
//
//  Created by Reza on 10/13/23.
//

import Foundation
import Matrix
import GeometryProcessing
import MetalKit

func proceedOptimization(_ proceedNum: Int = 1) {
    for _ in 0..<proceedNum {
        guard (converged == 0) else { break }
        print("Iteration \(iterNum):")
        converged = optimizer.solve(1)
        iterNum = optimizer.getIterNum()
    }
}

func updateViewerData_meshEdges() {
    fatalError("To be implemented")
}

func updateViewData_seam(_ V: inout Matd, _ F: inout Matd, _ UV: inout Matd) {
    fatalError("To be implemented")
}

func updateViewerData_distortion() {
    fatalError("To be implemented")
}

func updateViewerData() {
    
}

func saveScreenshot(_ filePath: String, _ scale: Double = 1.0, _ writeGif: Bool = false, _ writePNG: Bool = true) {
    fatalError("To be implemented")
}

func saveInfoForPresent(_ fileName: String = "info.text") {
    fatalError("To be implemented")
}

func toggleOptimization() {
    optimization_on = !optimization_on
    if (optimization_on) {
        if (converged != 0) {
            optimization_on = false
            print("optimization converged.")
        } else {
            if (!headlessMode && iterNum == 0) {
                // Start gif begin
                // save screenshop
            }
            print("start/resume optimization, press again to pause.")
        }
    } else {
        print("pause optimization, press again to resume.")
    }
}

func computeOptPicked(_ energyChanges0: [Pair<Double, Double>],
                      _ energyChanges1: [Pair<Double, Double>],
                      _ lambda: Double) -> Int {
    assert(!energyChanges0.isEmpty)
    assert(!energyChanges1.isEmpty)
    assert(lambda >= 0.0 && lambda <= 1.0)
    
    var minEChange0: Double = __DBL_MAX__
    
    for ecI in 0..<energyChanges0.count {
        if ((energyChanges0[ecI].first == __DBL_MAX__) || (energyChanges0[ecI].second == __DBL_MAX__)) {
            continue
        }
        let EwChange: Double = energyChanges0[ecI].first * (1.0 - lambda) + energyChanges0[ecI].second * lambda
        if (EwChange < minEChange0) {
            minEChange0 = EwChange
        }
    }
    
    var minEChange1 = __DBL_MAX__
    for ecI in 0..<energyChanges1.count {
        if (energyChanges1[ecI].first == __DBL_MAX__ || energyChanges1[ecI].second == __DBL_MAX__) {
            continue
        }
        let EwChange: Double = energyChanges1[ecI].first * (1.0 - lambda) + energyChanges1[ecI].second * lambda
        if (EwChange < minEChange1) {
            minEChange1 = EwChange
        }
    }
    
    assert(minEChange0 != __DBL_MAX__ || minEChange1 != __DBL_MAX__)
    
    return (minEChange0 > minEChange1 ? 1 : 0)
}

func computeBestCand(_ energyChanges: [Pair<Double, Double>], _ lambda: Double, _ bestEChange: inout Double) -> Int {
    assert(lambda >= 0.0 && lambda <= 1.0)
    
    bestEChange = __DBL_MAX__
    var id_minEChange: Int = -1
    for ecI in 0..<energyChanges.count {
        if (energyChanges[ecI].first == __DBL_MAX__ || energyChanges[ecI].second == __DBL_MAX__) {
            continue
        }
        let EwChange: Double = energyChanges[ecI].first * (1.0 - lambda) + energyChanges[ecI].second * lambda
        if (EwChange < bestEChange) {
            bestEChange = EwChange
            id_minEChange = ecI
        }
    }
    
    return id_minEChange
}

func checkCand(_ energyChanges: [Pair<Double, Double>]) -> Bool {
    for candI in energyChanges {
        if (candI.first < 0.0 || candI.second < 0.0) {
            return true
        }
    }
    
    var minEChange: Double = __DBL_MAX__
    for candI in energyChanges {
        if (candI.first < minEChange) {
            minEChange = candI.first
        }
        if (candI.second < minEChange) {
            minEChange = candI.second
        }
    }
    
    print("candidates not valid, minEChange: \(minEChange)")
    return false
}

func updateLambda(_ measure_bound: Double, _ lambda_SD: UnsafeMutablePointer<Double> = energyParams[0], _ kappa: Double = 1.0, _ kappa2: Double = 1.0) -> Double {
    var lambda_SD = lambda_SD.pointee
    lambda_SD = max(0.0, kappa * (measure_bound - (upperBound - convTol_upperBound / 2.0)) + kappa2 * lambda_SD / (1.0 - lambda_SD))
    return lambda_SD / (1.0 + lambda_SD)
}

func updateLambda_stationaryV(_ cancelMomentum: Bool = true, _ checkConvergence: Bool = false) -> Bool {
    let edgeLengths: Matd = edge_lengths(vertices: triSoup[channel_result].V_rest, faces: triSoup[channel_result].F)
    let eps_E_se: Double = 1.0e-3 * edgeLengths.minCoeff() / triSoup[channel_result].virtualRadius
    
    // measurement and energy value computation
    let E_SD: Double = optimizer.getLastEnergyVal(true) / energyParams[0].pointee
    // energyParams[0] = 0.01
    var E_se: Double = .zero
    triSoup[channel_result].computeSeamSparsity(&E_se)
    E_se /= triSoup[channel_result].virtualRadius
    var stretch_l2: Double = .zero
    var stretch_inf: Double = .zero
    var stretch_shear: Double = .zero
    var compress_inf: Double = .zero
    triSoup[channel_result].computeStandardStretch(&stretch_l2, &stretch_inf, &stretch_shear, &compress_inf)
    let measure_bound: Double = E_SD
    let eps_lambda: Double = min(1.0e-3, abs(updateLambda(measure_bound) - energyParams[0].pointee))
    
    // TODO: ? Stop when first violates bounds from feasible, don't go to best feasible. check after each merge whether distortion is violated
    // oscillation detection
    if (iterNum != lastStationaryIterNum) {
        // not a roll back config
        let lambda: Double = 1.0 - energyParams[0].pointee
        var oscillate: Bool = false
        if configs_stationaryV.allSatisfy({$0.key < E_se}) {
            // all less than E_se
            if (!configs_stationaryV.isEmpty) {
                // use largest element
                let largest = configs_stationaryV.max(by: { $0.key < $1.key})
                if (abs(largest!.key - E_se) < eps_E_se) {
                    for lambdaI in largest!.value {
                        if (abs(lambdaI.first - lambda) < eps_lambda &&
                            abs(lambdaI.second - E_SD) < eps_E_se) {
                            oscillate = true
                            break
                        }
                    }
                }
            }
        } else if configs_stationaryV.allSatisfy({!($0.key < E_se)}) {
            // all not less than E_se
            let low = configs_stationaryV.sorted(by: {$0.key < $1.key}).first(where: { $0.key >= E_se})!
            if (abs(low.key - E_se) < eps_E_se) {
                for lambdaI in low.value {
                    if (abs(lambdaI.first - lambda) < eps_lambda &&
                        abs(lambdaI.second - E_SD) < eps_E_se) {
                        oscillate = true
                        break
                    }
                }
            }
        } else {
            let sorted = configs_stationaryV.sorted(by: {$0.key < $1.key})
            let low = sorted.firstIndex(where: {E_se <= $0.key})!
            let prev = low - 1
            
            if (abs(sorted[low].key - E_se) < eps_E_se) {
                for lambdaI in sorted[low].value {
                    if (abs(lambdaI.first - lambda) < eps_lambda &&
                        abs(lambdaI.second - E_SD) < eps_E_se) {
                        oscillate = true
                        break
                    }
                }
            }
            if (!oscillate && (abs(sorted[prev].key - E_se) < eps_E_se)) {
                for lambdaI in sorted[prev].value {
                    if (abs(lambdaI.first - lambda) < eps_lambda &&
                        abs(lambdaI.second - E_SD) < eps_E_se) {
                        oscillate = true
                        break
                    }
                }
            }
        }
        
        // record best feasible UV map
        if (measure_bound <= upperBound) && (E_se < E_se_bestFeasible) {
            iterNum_bestFeasible = iterNum
            triSoup_bestFeasible = triSoup[channel_result].clone()
            E_se_bestFeasible = E_se
        }
        
        if (oscillate && (iterNum_bestFeasible >= 0)) {
            // arrive at the best feasible config again
            if (iterNum_bestFeasible != iterNum) {
                optimizer.setConfig(triSoup_bestFeasible, iterNum, optimizer.getTopoIter())
            }
            return false
        } else {
            if configs_stationaryV[E_se] != nil {
                configs_stationaryV[E_se]!.append(Pair(lambda, E_SD))
            } else {
                configs_stationaryV[E_se] = [Pair(lambda, E_SD)]
            }
        }
    }
    
    lastStationaryIterNum = iterNum
    
    // convergence check
    if (checkConvergence) {
        if (measure_bound <= upperBound) {
            if (!saved) {
                saved = true
            }
            
            if (measure_bound >= upperBound - convTol_upperBound) {
                if (iterNum_bestFeasible != iterNum) {
                    assert(iterNum_bestFeasible >= 0)
                    optimizer.setConfig(triSoup_bestFeasible, iterNum, optimizer.getTopoIter())
                }
                return false
            }
        }
    }
    
    // lambda update (dual update)
    energyParams[0].pointee = updateLambda(measure_bound)
    // energyParams[0] = 0.33
    // TODO: needs to be careful on lambda update space
    
    // critical lambda scheme
    if (checkConvergence) {
        // update lambda untill feasible update on T might be triggered
        if (measure_bound > upperBound) {
            // need to cut further, increase energyParams[0]
            
            if (!energyChanges_merge.isEmpty &&
                computeOptPicked(energyChanges_bSplit, energyChanges_merge, 1.0 - energyParams[0].pointee) == 1) {
                // still picking merge
                repeat {
                    energyParams[0].pointee = updateLambda(measure_bound)
                } while ((computeOptPicked(energyChanges_bSplit, energyChanges_merge, 1.0 - energyParams[0].pointee) == 1))
            }
            
            if (!checkCand(energyChanges_iSplit) && !checkCand(energyChanges_bSplit)) {
                // if filtering too stroing
                reQuery = true
            } else {
                var eDec_b: Double = .zero
                var eDec_i: Double = .zero
                assert(!(energyChanges_bSplit.isEmpty && energyChanges_iSplit.isEmpty))
                var id_pickingBSplit: Int = computeBestCand(energyChanges_bSplit, 1.0 - energyParams[0].pointee, &eDec_b)
                var id_pickingISplit: Int = computeBestCand(energyChanges_iSplit, 1.0 - energyParams[0].pointee, &eDec_i)
                while ((eDec_b > 0.0) && (eDec_i > 0.0)) {
                    energyParams[0].pointee = updateLambda(measure_bound)
                    id_pickingBSplit = computeBestCand(energyChanges_bSplit, 1.0 - energyParams[0].pointee, &eDec_b)
                    id_pickingISplit = computeBestCand(energyChanges_iSplit, 1.0 - energyParams[0].pointee, &eDec_i)
                }
                if (eDec_b <= 0.0) {
                    opType_queried = 0
                    path_queried = paths_bSplit[id_pickingBSplit]
                    newVertPos_queried = newVertPoses_bSplit[id_pickingBSplit]
                } else {
                    opType_queried = 1
                    path_queried = paths_iSplit[id_pickingISplit]
                    newVertPos_queried = newVertPoses_iSplit[id_pickingISplit]
                }
            }
        } else {
            var noOp: Bool = true
            for ecI in energyChanges_merge {
                if (ecI.first != __DBL_MAX__) {
                    noOp = false
                    break
                }
            }
            if (noOp) {
                energyParams[0].pointee = 1.0 - eps_lambda
                optimizer.updateEnergyData(true, false, false)
                if (iterNum_bestFeasible != iterNum) {
                    optimizer.setConfig(triSoup_bestFeasible, iterNum, optimizer.getTopoIter())
                }
                return false
            }
            
            // !!! also account for iSplit for this switch?
            if (computeOptPicked(energyChanges_bSplit, energyChanges_merge, 1.0 - energyParams[0].pointee) == 0) {
                // still picking split
                repeat {
                    energyParams[0].pointee = updateLambda(measure_bound)
                } while (computeOptPicked(energyChanges_bSplit, energyChanges_merge, 1.0 - energyParams[0].pointee) == 0)
            }
            
            var eDec_m: Double = .zero
            assert(!energyChanges_merge.isEmpty)
            var id_pickingMerge: Int = computeBestCand(energyChanges_merge, 1.0 - energyParams[0].pointee, &eDec_m)
            while (eDec_m > 0.0) {
                energyParams[0].pointee = updateLambda(measure_bound)
                id_pickingMerge = computeBestCand(energyChanges_merge, 1.0 - energyParams[0].pointee, &eDec_m)
            }
            opType_queried = 2
            path_queried = paths_merge[id_pickingMerge]
            newVertPos_queried = newVertPoses_merge[id_pickingMerge]
        }
    }
    
    // lambda value result sanity check
    if (energyParams[0].pointee > 1.0 - eps_lambda) {
        energyParams[0].pointee = 1.0 - eps_lambda
    }
    if (energyParams[0].pointee < eps_lambda) {
        energyParams[0].pointee = eps_lambda
    }
    
    optimizer.updateEnergyData(true, false, false)
    
    return true
}

func converge_preDrawFunc() {
    infoName = "finalResult"
    
    if (!bijectiveParam) {
        // perform exact solve
        optimizer.setAllowEDecRelTol(false)
        converged = 0
        optimizer.setPropogateFracture(false)
        while (converged == 0) {
            proceedOptimization(1000)
        }
    }
    
    updateViewerData()
    optimizer.flushEnergyFileOutput()
    optimizer.flushGradFileOutput()
    
    optimization_on = false
    print("optimization converged, with secpast")
    outerLoopFinished = true
}

@discardableResult
func preDrawFunc() -> Bool {
    if (optimization_on) {
        if (offlineMode) {
            while(converged == 0) {
                proceedOptimization()
            }
            
        } else {
            proceedOptimization()
        }
        
        updateViewerData()
        
        if (converged != 0) {
            saveInfo_postDraw = true
            
            var stretch_l2: Double = .zero
            var stretch_inf: Double = .zero
            var stretch_shear: Double = .zero
            var compress_inf: Double = .zero
            triSoup[channel_result].computeStandardStretch(&stretch_l2, &stretch_inf, &stretch_shear, &compress_inf)
            let measure_bound: Double = optimizer.getLastEnergyVal(true) / energyParams[0].pointee
            
            switch (methodType) {
            case .MT_EBCUTS:
                if (measure_bound <= upperBound) {
                    infoName = "finalResult"
                    // perform exact solve
                    optimizer.setAllowEDecRelTol(false)
                    converged = 0
                    while (converged == 0) {
                        proceedOptimization(1000)
                    }
                    updateViewerData()
                    
                    optimization_on = false
                    print("optimization converged, with sec past")
                    outerLoopFinished = true
                } else {
                    infoName = String(iterNum)
                    
                    // continue to make geometry image cuts
                    assert(optimizer.createFracture(fracThres, 0, false))
                    converged = 0
                }
                optimizer.flushEnergyFileOutput()
                optimizer.flushGradFileOutput()
                break
                
            case .MT_OPTCUTS_NODUAL, .MT_OPTCUTS:
                infoName = String(iterNum)
                if (converged == 2) {
                    converged = 0
                    return false
                }
                
                if ((methodType == .MT_OPTCUTS) && (measure_bound <= upperBound)) {
                    if (!preDrawSaved) {
                        preDrawSaved = true
                    }
                }
                
                // if necessary, turn on scaffolding for random one point initial cut
                if (!optimizer.isScaffolding() && bijectiveParam && rand1PinitCut) {
                    optimizer.setScaffolding(true)
                }
                
                var E_se: Double = .zero
                triSoup[channel_result].computeSeamSparsity(&E_se)
                E_se /= triSoup[channel_result].virtualRadius
                let E_SD: Double = optimizer.getLastEnergyVal(true) / energyParams[0].pointee
                
                print("\(iterNum): \(E_SD) \(E_se) \(triSoup[channel_result].V_rest.rows)")
                optimizer.flushEnergyFileOutput()
                optimizer.flushGradFileOutput()
                
                // continue to split boundary
                if (methodType == .MT_OPTCUTS && !updateLambda_stationaryV()) {
                    // oscillation detected
                    converge_preDrawFunc()
                } else {
                    if (optimizer.createFracture(fracThres, 0, topoLineSearch)) {
                        converged = 0
                    } else {
                        // if no boundary op, try interior split if split is the current best boundary op
                        if ((measure_bound > upperBound) &&
                            optimizer.createFracture(fracThres, 0, topoLineSearch, true)) {
                            converged = 0
                        } else {
                            if ((methodType == .MT_OPTCUTS_NODUAL) || (!updateLambda_stationaryV(false, true))) {
                                // all converged
                                converge_preDrawFunc()
                            } else {
                                // split or merge after lambda update
                                if (reQuery) {
                                    filterExp_in += log(2.0) / log(Double(inSplitTotalAmt))
                                    filterExp_in = min(1.0, filterExp_in)
                                    
                                    while (!optimizer.createFracture(fracThres, 0, topoLineSearch, true)) {
                                        filterExp_in += log(2.0) / log(Double(inSplitTotalAmt))
                                        filterExp_in = min(1.0, filterExp_in)
                                    }
                                    reQuery = false
                                    // TODO: set filtering param back?
                                } else {
                                    optimizer.createFracture(opType_queried, path_queried, newVertPos_queried, topoLineSearch)
                                }
                                opType_queried = -1
                                converged = 0
                            }
                        }
                    }
                }
                break
                
            case .MT_DISTMIN:
                converge_preDrawFunc()
                break
            }
        }
    } else {
        if (isCapture3D && (capture3DI < 2)) {
            // change view accorgingly
            /*
            let rotDeg: Double = ((Double(capture3DI) < 8.0) ? ((.pi / 2.0) * (Double(capture3DI) / 2)) : .pi / 2.0)
            var rotAxis: Vec3<Float> = .init(0, 1.0, 0)
            if ((capture3DI / 2) == 4) {
                rotAxis = .init(1.0, 0, 0)
            } else if ((capture3DI / 2) == 5) {
                rotAxis = .init(-1.0, 0, 0)
            } */
            viewChannel = channel_result
            viewUV = false
            showSeam = true
            isLighting = false
            showTexture = ((capture3DI % 2) == 0 ? false : true)
            showDistortion = 2 - capture3DI % 2
            updateViewerData()
        }
    }
    return false
}

@discardableResult
func postDrawFunc() -> Bool {
    if (offlineMode && (iterNum == 0)) {
        toggleOptimization()
    }
    
    if (saveInfo_postDraw) {
        saveInfo_postDraw = false
    }
    
    if (outerLoopFinished) {
        if (offlineMode) {
            //exit(EXIT_SUCCESS)
        } else {
            outerLoopFinished = false
        }
    }
    
    return false
}
