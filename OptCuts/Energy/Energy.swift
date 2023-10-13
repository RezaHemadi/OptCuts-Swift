//
//  Energy.swift
//  OptCuts
//
//  Created by Reza on 6/20/23.
//

import Foundation
import Matrix
import GeometryProcessing

/// An abstract class for energy terms in the objective of an optimization problem
 class Energy {
    // MARK: - Properties
    private let needRefactorize: Bool
    
    // MARK: - Initialization
     init(p_needRefactorize: Bool) {
        self.needRefactorize = p_needRefactorize
    }
    
    // MARK: - Methods
     func getNeedRefactorize() -> Bool {
        return needRefactorize
    }
    
     func computeEnergyVal(_ data: TriMesh,
                          _ energyVal: inout Double,
                          _ uniformWeight: Bool = false) {
        
        var energyValPerElem: Vec<Double> = .init()
        getEnergyValPerElem(data, &energyValPerElem, uniformWeight)
        energyVal = energyValPerElem.sum()
    }
    
     func getEnergyValPerElem(_ data: TriMesh,
                             _ energyValPerElem: inout Vec<Double>,
                             _ uniformWeight: Bool = false) {
        fatalError("To be implemented")
    }
    
     func getEnergyValByElemID(_ data: TriMesh,
                              _ elemI: Int,
                              _ energyVal: inout Double,
                              _ uniformWeight: Bool = false) {
        fatalError("To be implemented")
    }
    
     func computeGradient(_ data: TriMesh,
                         _ gradient: inout Vec<Double>,
                         _ uniformWeight: Bool = false) {
        fatalError("To be implemented")
    }
    
     func computeHessian(_ data: TriMesh,
                        _ V: inout Vec<Double>,
                        _ I: inout Vec<Int>,
                        _ J: inout Vec<Int>,
                        _ uniformWeights: Bool = false) {
        fatalError("To be implemented")
    }
    
     func computeHessian(_ data: TriMesh,
                        _ Hessian: inout Mat<Double>,
                        _ uniformWeight: Bool = false) {
        fatalError("Implement in subclass")
    }
    
     func checkEnergyVal(_ data: TriMesh) {
        fatalError("To be implemented")
    }
    
    /// check with finite difference method, according to energyVal
     func checkGradient(data: TriMesh) {
        print("Checking energy gradient computation...")
        
        var energyVal0: Double = .zero
        computeEnergyVal(data, &energyVal0)
        let h: Double = 1.0e-8 * avg_edge_length(data.V, data.F)
        var perturbed = data
        var gradient_finiteDiff = Vec<Double>(data.V.rows * 2)
        for vI in 0..<data.V.rows {
            for dimI in 0..<2 {
                perturbed.V = .init(data.V, data.V.rows, data.V.cols)
                perturbed.V[vI, dimI] += h
                var energyVal_perturbed: Double = 0.0
                computeEnergyVal(perturbed, &energyVal_perturbed)
                gradient_finiteDiff[vI * 2 + dimI] = (energyVal_perturbed - energyVal0) / h
            }
            
            if ((vI + 1) % 100 == 0) {
                print("\(vI + 1) \(data.V.rows) vertices computed.")
            }
        }
        for fixedVI in data.fixedVert {
            gradient_finiteDiff[2 * fixedVI] = 0.0
            gradient_finiteDiff[2 * fixedVI + 1] = 00
        }
        
        var gradient_symmbolic = Vec<Double>()
        computeGradient(data, &gradient_symmbolic)
        
        var difVec = Vec<Double>()
        let dif_L2: Double = difVec.norm()
        let relErr: Double = dif_L2 / gradient_finiteDiff.norm()
        
        print("L2 dist = \(dif_L2), relErr = \(relErr)")
    }
    
    /// check with finite difference method, according to gradient
     func checkHessian(_ data: TriMesh,
                      _ useTriplet: Bool = false) {
        print("checking energy hessian computation...")
        
        var gradient0 = Vec<Double>()
        computeGradient(data, &gradient0)
        let h: Double = 1.0e-8 * avg_edge_length(data.V, data.F)
        var perturbed = data
        var hessian_finiteDiff = SparseMatrix<Double>()
        hessian_finiteDiff.resize(data.V.rows * 2, data.V.rows * 2)
        var triplets: [Tripletd] = []
        for vI in 0..<data.V.rows {
            if (data.fixedVert.contains(vI)) {
                triplets.append(.init(i: vI * 2, j: vI * 2, value: 1.0))
                triplets.append(.init(i: vI * 2 + 1, j: vI * 2 + 1, value: 1.0))
                continue
            }
            
            for dimI in 0..<2 {
                perturbed.V = .init(data.V, data.V.rows, data.V.cols)
                perturbed.V[vI, dimI] += h
                var gradient_perturbed = Vec<Double>()
                computeGradient(perturbed, &gradient_perturbed)
                var hessian_colI: Vec<Double> = (gradient_perturbed - gradient0) / h
                var colI: Int = vI * 2 + dimI
                for rowI in 0..<(data.V.rows * 2) {
                    if (data.fixedVert.contains(rowI / 2)) {
                        continue
                    }
                    
                    triplets.append(.init(i: rowI, j: colI, value: hessian_colI[rowI]))
                }
            }
            
            if ((vI + 1) % 100 == 0) {
                print("\(vI + 1) \(data.V.rows) vertices computed")
            }
        }
        
        hessian_finiteDiff.setFromTriplets(triplets)
        
        var hessian_symbolic = SparseMatrix<Double>()
        assert(useTriplet)
        
        var I = Vec<Int>()
        var J = Vec<Int>()
        var V = Vec<Double>()
        computeHessian(data, &V, &I, &J)
        triplets.removeAll()
        triplets.reserveCapacity(V.count)
        for entryI in 0..<V.count {
            triplets.append(.init(i: I[entryI], j: J[entryI], value: V[entryI]))
        }
        hessian_symbolic.resize(data.V.rows * 2, data.V.rows * 2)
        hessian_symbolic.setFromTriplets(triplets)
        
        var diffMtr: SparseMatrix<Double> = hessian_symbolic - hessian_finiteDiff
        let dif_L2: Double = diffMtr.norm()
        let relErr: Double = dif_L2 / hessian_finiteDiff.norm()
        
        print("L2 dist = \(dif_L2), relErr = \(relErr)")
    }
    
     func initStepSize(_ data: TriMesh,
                      _ searchDir: Vec<Double>,
                      _ stepSize: inout Double) {
        
    }
}
