//
//  EigenLibSolver.swift
//  OptCuts
//
//  Created by Reza on 6/20/23.
//

import Foundation
import Matrix

class EigenLibSolver<vectorTypeI: Vector, vectorTypeS: Vector>: LinSysSolver<vectorTypeI, vectorTypeS> {
    typealias Base = LinSysSolver<vectorTypeI, vectorTypeS>
    
    // MARK: - Properties
    private var useDense: Bool = false
    private var coefMtr_dense: Mat<Double> = .init()
    private var LDLT: LDLT = .init()
    private var coefMtr: SparseMatrix<Double> = .init()
    private var simplicalLDLT: SimplicalLDLT = .init()
    
    // MARK: - Methods
     override func set_type(_ threadAmt: Int,
                                  _ _mtype: Int,
                                  _ is_upper_half: Bool = false) {
        
        // TODO: support more matrix types, currently only SPD
        useDense = false
        // TODO: move to base class and support for CHOLMOD
    }
    
     override func set_pattern(_ vNeighbor: [Set<Int>],
                              _ fixedVert: Set<Int>) {
        if (useDense) {
            numRows = vNeighbor.count * 2
            coefMtr_dense.resize(numRows, numRows)
        } else {
            super.set_pattern(vNeighbor, fixedVert)
            
            // TODO: directly save into mtr
            coefMtr.resize(numRows, numRows)
            coefMtr.reserve(ja.count)
            ia.array() -= 1
            ja.array() -= 1
            coefMtr.setInnerIndices(ja.values)
            coefMtr.setOuterIndices(ia.values)
        }
    }
    
     override func set_pattern(_ mtr: SparseMatrix<Double>) { // Note: mtr must be SPD
        if (useDense) {
            coefMtr_dense = Matd(mtr)
        } else {
            coefMtr = mtr
        }
    }
    
     override func update_a(_ II: vectorTypeI,
                           _ JJ: vectorTypeI,
                           _ SS: vectorTypeS) where vectorTypeI.Element == Int, vectorTypeS.Element == Double {
        
        if (useDense) {
            coefMtr_dense.setZero()
            for i in 0..<II.count {
                coefMtr_dense[II[i], JJ[i]] += SS[i]
            }
        } else {
            super.update_a(II, JJ, SS)
            
            // TODO: directly save into coefMtr
            coefMtr.setValues(a.values)
        }
    }
    
    override func analyze_pattern() {
        if (!useDense) {
            //simplicalLDLT.analyzePattern(coefMtr)
            //simplicalLDLT = .init()
            //simplicalLDLT = .init()
            //simplicalLDLT.factorize(coefMtr)
            assert(simplicalLDLT.info == .success)
        }
    }
    
    override func factorize() -> Bool {
        var succeeded: Bool = false
        if (useDense) {
            LDLT.compute(coefMtr_dense)
            succeeded = (LDLT.info == .success)
        } else {
            simplicalLDLT.factorize(coefMtr)
            succeeded = (simplicalLDLT.info == .success)
        }
        
        assert(succeeded)
        return succeeded
    }
    
    override func solve(_ rhs: inout Vec<Double>,
                               _ result: inout Vec<Double>) {
        
        if useDense {
            LDLT.solve(rhs, &result)
            assert(LDLT.info == .success)
        } else {
            result = simplicalLDLT.solve(rhs)
            assert(simplicalLDLT.info == .success)
        }
    }
    
    override func coeffMtr(_ rowI: Int, _ colI: Int) -> Double {
        if (useDense) {
            return coefMtr_dense[rowI, colI]
        } else {
            return coefMtr[rowI, colI]
        }
    }
    
    override func setZero() {
        // TODO: useDense
        // TODO: directly manipulate valuesPtr without a
        super.setZero()
        coefMtr.setValues(a.values)
    }
    
    override func setCoeff(_ rowI: Int, _ colI: Int, _ val: Double) {
        // TODO: useDense
        // TODO: directly manipulate valuesPtr without a
        // TODO: faster 0(1) indices!!
        
        if (rowI <= colI) {
            assert(rowI < IJ2aI.count)
            if let aI = IJ2aI[rowI][colI] {
                a[aI] = val
                coefMtr.setValue(val, at: aI)
            }
        }
    }
    
    override func addCoeff(_ rowI: Int, _ colI: Int, _ val: Double) {
        
        // TODO: useDense
        // TODO: directly manipulate valuePtr without a
        // TODO: faster 0(1) indices!!
        if (rowI <= colI) {
            assert(rowI < IJ2aI.count)
            if let aI = IJ2aI[rowI][colI] {
                a[aI] += val
                coefMtr.addValue(val, at: aI)
            }
        }
    }
}
