//
//  OPSparseMatrix.swift
//  OptCuts
//
//  Created by Reza on 10/14/23.
//

import Foundation
import Matrix
import Accelerate

class OPSparseMatrix {
    static private let AllocateCount: Int = 400_000
    // MARK: - Properties
    var innerIndices: UnsafeMutablePointer<Int32>!
    var outerStarts: UnsafeMutablePointer<Int>!
    var values: UnsafeMutablePointer<Double>!
    var rows: Int32!
    var cols: Int32!
    
    // MARK: - Initializers
    init() {
        innerIndices = .allocate(capacity: Self.AllocateCount)
    }
    
    deinit {
        innerIndices.deallocate()
    }
    
    subscript(i: Int, j: Int) -> Double {
        get {
            fatalError()
        }
    }
    
    // MARK: - Methods
    func setOuterIndices(_ indices: UnsafeMutablePointer<Int>, count: Int) {
        outerStarts = indices
    }
    
    func setInnerIndices(_ indices: UnsafeMutablePointer<Int>, count: Int) {
        assert(count <= Self.AllocateCount)
        // TODO: test concurrent
        /*
        DispatchQueue.concurrentPerform(iterations: count) { i in
            (innerIndices + i).initialize(to: Int32(indices[i]))
        }*/
        
        for i in 0..<count {
            //(innerIndices + i).initialize(to: Int32(indices[i]))
            innerIndices[i] = Int32(indices[i])
        }
    }
    
    func setValues(values: UnsafeMutablePointer<Double>, count: Int) {
        self.values = values
    }
    
    func resize(_ rows: Int, _ cols: Int) {
        self.rows = Int32(rows)
        self.cols = Int32(cols)
    }
    
    func setValue(_ val: Double, at: Int) {
        fatalError()
    }
    
    func addValue(_ val: Double, at: Int) {
        fatalError()
    }
}

extension SimplicalLDLT {
    func factorize(_ mtr: OPSparseMatrix) {
        //var rowIndices: [Int32] = mtr.innerIndices.map({ Int32($0) })
        //var columnStarts: [Int] = mtr.outerStarts
        
        let a = SparseMatrix_Double(structure: structure!,
                                    data: mtr.values)
        factorization = SparseFactor(SparseFactorizationLDLT, a)
        
        /*
        let llt: SparseOpaqueFactorization_Double = values.withUnsafeMutableBufferPointer { valuesPtr in
            let a = SparseMatrix_Double(structure: structure,
                                        data: valuesPtr.baseAddress!)
            return SparseFactor(SparseFactorizationLDLT, a)
        }
        
        factorization = llt*/
    }
}

extension SimplicalLDLT {
    func analyzePattern(_ mtr: OPSparseMatrix) {
        
        if structure == nil {
            var attributes = SparseAttributes_t()
            attributes.triangle = SparseLowerTriangle
            attributes.kind = SparseSymmetric
            
            structure = SparseMatrixStructure(rowCount: mtr.rows,
                                              columnCount: mtr.cols,
                                              columnStarts: mtr.outerStarts,
                                              rowIndices: mtr.innerIndices,
                                              attributes: attributes,
                                              blockSize: 1)
        } else {
            if structure!.rowCount != mtr.rows {
                structure!.rowCount = mtr.rows
            }
            if structure!.columnCount != mtr.cols {
                structure!.columnCount = mtr.cols
            }
            structure!.columnStarts = mtr.outerStarts
            structure!.rowIndices = mtr.innerIndices
        }
    }
}
