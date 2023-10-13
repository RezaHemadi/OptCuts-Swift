//
//  LinSysSolver.swift
//  OptCuts
//
//  Created by Reza on 6/20/23.
//

import Foundation
import Matrix

 class LinSysSolver<vectorTypeI: Vector, vectorTypeS: Vector> {
    // MARK: - Properties
     var numRows: Int = 0
     var ia: Vec<Int> = .init()
     var ja: Vec<Int> = .init()
     var IJ2aI: [[Int:Int]] = []
     var a: Vec<Double> = .init()
    
    // MARK: - Methods
     func set_type(_ threadAmt: Int,
                  _ _mtype: Int,
                  _ is_uppder_half: Bool = false) {
        fatalError("implement in subclass")
    }
    
     func set_pattern(_ vNeighbor: [Set<Int>],
                     _ fixedVert: Set<Int>) {
        
        numRows = vNeighbor.count * DIM
        ia.resize(vNeighbor.count * DIM + 1)
        ia[0] = 1 // 1 + nnz above row i
        ja.resize(0) // colI of each element
        IJ2aI = .init(repeating: [:], count: vNeighbor.count * DIM) // map from matrix index to ja index
        for rowI in 0..<vNeighbor.count {
            if (!fixedVert.contains(rowI)) {
                let oldSize_ja: Int = ja.count
                IJ2aI[rowI * DIM][rowI * DIM] = oldSize_ja
                IJ2aI[rowI * DIM][rowI * DIM + 1] = oldSize_ja + 1
                if (DIM == 3) {
                    IJ2aI[rowI * DIM][rowI * DIM + 2] = oldSize_ja + 2
                }
                ja.conservativeResize(oldSize_ja + DIM)
                ja[oldSize_ja] = rowI * DIM + 1
                ja[oldSize_ja + 1] = rowI * DIM + 2
                if (DIM == 3) {
                    ja[oldSize_ja + 2] = rowI * DIM + 3
                }
                
                var nnz_rowI: Int = 1
                for colI in vNeighbor[rowI] {
                    if (!fixedVert.contains(colI)) {
                        if (colI > rowI) {
                            // only the lower-left part
                            // colI > rowI means upper-right, but we are preparing CSR here
                            // in a row-major maner and CHOLMOD is actually column-major
                            let oldSize_ja_temp: Int = ja.count
                            IJ2aI[rowI * DIM][colI * DIM] = oldSize_ja_temp
                            IJ2aI[rowI * DIM][colI * DIM + 1] = oldSize_ja_temp + 1
                            if (DIM == 3) {
                                IJ2aI[rowI * DIM][colI * DIM + 2] = oldSize_ja_temp + 2
                            }
                            ja.conservativeResize(oldSize_ja_temp + DIM)
                            ja[oldSize_ja_temp] = colI * DIM + 1
                            ja[oldSize_ja_temp + 1] = colI * DIM + 2
                            if (DIM == 3) {
                                ja[oldSize_ja_temp + 2] = colI * DIM + 3
                            }
                            nnz_rowI += 1
                        }
                    }
                }
                // another row for y,
                // excluding the left-bottom entry on the diagonal band
                IJ2aI[rowI * DIM + 1] = IJ2aI[rowI * DIM]
                for (key, value) in IJ2aI[rowI * DIM + 1] {
                    IJ2aI[rowI * DIM + 1][key] = value + (nnz_rowI * DIM - 1) // TODO: Verify code
                }
                ja.conservativeResize(ja.count + nnz_rowI * DIM - 1)
                ja.bottomRows(nnz_rowI * DIM - 1) <<== ja.block(oldSize_ja + 1, 0, nnz_rowI * DIM - 1, 1)
                
                if (DIM == 3) {
                    // third row for z
                    IJ2aI[rowI * DIM + 2] = IJ2aI[rowI * DIM + 1]
                    for (key, value) in IJ2aI[rowI * DIM + 2] {
                        IJ2aI[rowI * DIM + 2][key] = value + (nnz_rowI * DIM - 2)
                    }
                    ja.conservativeResize(ja.count + nnz_rowI * DIM - 2)
                    ja.bottomRows(nnz_rowI * DIM - 2) <<== ja.block(oldSize_ja + 2, 0, nnz_rowI * DIM - 2, 1)
                    
                    IJ2aI[rowI * DIM + 2].removeValue(forKey: rowI * DIM)
                    IJ2aI[rowI * DIM + 2].removeValue(forKey: rowI * DIM + 1)
                }
                IJ2aI[rowI * DIM + 1].removeValue(forKey: rowI * DIM)
                
                ia[rowI * DIM + 1] = ia[rowI * DIM] + nnz_rowI * DIM
                ia[rowI * DIM + 2] = ia[rowI * DIM + 1] + nnz_rowI * DIM - 1
                if (DIM == 3) {
                    ia[rowI * DIM + 3] = ia[rowI * DIM + 2] + nnz_rowI * DIM - 2
                }
            } else {
                let oldSize_ja: Int = ja.count
                IJ2aI[rowI * DIM][rowI * DIM] = oldSize_ja
                IJ2aI[rowI * DIM + 1][rowI * DIM + 1] = oldSize_ja + 1
                if (DIM == 3) {
                    IJ2aI[rowI * DIM + 2][rowI * DIM + 2] = oldSize_ja + 2
                }
                ja.conservativeResize(oldSize_ja + DIM)
                ja[oldSize_ja] = rowI * DIM + 1
                ja[oldSize_ja + 1] = rowI * DIM + 2
                if (DIM == 3) {
                    ja[oldSize_ja + 2] = rowI * DIM + 3
                }
                ia[rowI * DIM + 1] = ia[rowI * DIM] + 1
                ia[rowI * DIM + 2] = ia[rowI * DIM + 1] + 1
                if (DIM == 3) {
                    ia[rowI * DIM + 3] = ia[rowI * DIM + 2] + 1
                }
            }
        }
        a.resize(ja.count)
    }
    
     func set_pattern(_ mtr: SparseMatrix<Double>) { // Note: mtr must be SPD
        return
    }
    
     func update_a(_ II: vectorTypeI,
                  _ JJ: vectorTypeI,
                  _ SS: vectorTypeS) where vectorTypeI.Element == Int, vectorTypeS.Element == Double {
        
        // TODO: faster 0(1) indices!!
        assert(II.count == JJ.count)
        assert(II.count == SS.count)
        
        a.resize(ja.count)
        a.setZero()
        for tripletI in 0..<II.count {
            let i: Int = II[tripletI]
            let j: Int = JJ[tripletI]
            if (i <= j) {
                assert(i < IJ2aI.count)
                if let value = IJ2aI[i][j] {
                    a[value] += SS[tripletI]
                }
            }
        }
    }
    
     func update_a(mtr: SparseMatrix<Double>) {
        assert(false, "please implement in subclass!")
    }
    
     func analyze_pattern() {
        fatalError("implement in subclass")
    }
    
    @discardableResult
     func factorize() -> Bool {
        fatalError("implement in subclass")
    }
    
     func solve(_ rhs: inout Vec<Double>,
               _ result: inout Vec<Double>) {
        fatalError("implement in subclass!")
    }
    
     func multiply(_ x: Vec<Double>,
                  _ Ax: inout Vec<Double>) {
        
        assert(x.count == numRows)
        assert(IJ2aI.count == numRows)
        
        Ax.resize(numRows)
        Ax.setZero()
        for rowI in 0..<numRows {
            for (key, value) in IJ2aI[rowI] {
                Ax[rowI] += Double(value) * x[key]
                if (rowI != key) {
                    Ax[key] += Double(value) * x[rowI]
                }
            }
        }
    }
    
     func coeffMtr(_ rowI: Int, _ colI: Int) -> Double {
        var rowI = rowI
        var colI = colI
        
        if (rowI > colI) {
            // return only upper right part for symmetric matrix
            let temp: Int = rowI
            rowI = colI
            colI = temp
        }
        assert(rowI < IJ2aI.count)
        if let value = IJ2aI[rowI][colI] {
            return a[value]
        } else {
            return 0.0
        }
    }
    
     func getCoeffMtr(_ mtr: inout SparseMatrix<Double>) {
        mtr.resize(numRows, numRows)
        mtr.setZero()
        mtr.reserve(a.count * 2 - numRows)
        var triplets: [Tripletd] = []
        triplets.reserveCapacity(a.count * 2 - numRows)
        for rowI in 0..<numRows {
            for (key, value) in IJ2aI[rowI] {
                triplets.append(.init(i: rowI, j: key, value: a[value]))
                if (rowI != key) {
                    triplets.append(.init(i: key, j: rowI, value: a[value]))
                }
            }
        }
        mtr.setFromTriplets(triplets)
    }
    
     func setCoeff(_ rowI: Int, _ colI: Int, _ val: Double) {
        // TODO: faster 0(1) indices !!
        
        if (rowI <= colI) {
            assert(rowI < IJ2aI.count)
            if let value = IJ2aI[rowI][colI] {
                a[value] = val
            }
        }
    }
    
     func setZero() {
        a.setZero()
    }
    
     func addCoeff(_ rowI: Int, _ colI: Int, _ val: Double) {
        // TODO: faster 0(1) indices!!
        
        if (rowI <= colI) {
            assert(rowI < IJ2aI.count)
            if let value = IJ2aI[rowI][colI] {
                a[value] += val
            }
        }
    }
    
     func getNumRows() -> Int {
        return numRows
    }
    
     func getNumNonzeroes() -> Int {
        return a.count
    }
    
     func getIJ2aI() -> [[Int: Int]] {
        return IJ2aI
    }
    
    // Note: in original c++ code this is returned by reference
     func get_ia() -> Vec<Int> {
        return ia
    }
    
    // Note: in original c++ code this is returned by reference
     func get_ja() -> Vec<Int> {
        return ja
    }
    
    // Note: in original c++ code this is returned by reference
     func get_a() -> Vec<Double> {
        return a
    }
}
