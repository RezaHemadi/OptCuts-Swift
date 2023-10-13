//
//  CHOLMODSolver.swift
//  OptCuts
//
//  Created by Reza on 6/20/23.
//

import Foundation
import Matrix

class CHOLMODSolver<vectorTypeI: Vector, vectorTypeS: Vector>: LinSysSolver<vectorTypeI, vectorTypeS> {
    typealias Base = LinSysSolver<vectorTypeI, vectorTypeS>
    
    // MARK: - Properties
    
    // MARK: - Initialization
    override init() {
        fatalError("To be implemented")
    }
    
    deinit {
        fatalError("To be implemented")
    }
    
    // MARK: - Methods
    func set_type(threadAmt: Int, _mtype: Int, is_upper_half: Bool = false) {
        fatalError("To be implemented")
    }
    
    func set_pattern(vNeighbor: [Set<Int>],
                     fixedVert: Set<Int>) {
        fatalError("To be implemented")
    }
    
    func set_pattern(mtr: inout SparseMatrix<Double>) { // Note: mtr must be SPD
        fatalError("To be implemented")
    }
    
    func update_a(II: vectorTypeI,
                  JJ: vectorTypeI,
                  SS: vectorTypeS) {
        fatalError("To be implemented")
    }
    
    override func update_a(mtr: SparseMatrix<Double>) {
        fatalError("To be implemented")
    }
    
    override func analyze_pattern() {
        fatalError("To be implemented")
    }
    
    override func factorize() -> Bool {
        fatalError("To be implemented")
    }
    
    func solve(rhs: inout Vec<Double>,
               result: inout Vec<Double>) {
        fatalError("To be implemented")
    }
    
    func multiply(x: Vec<Double>, Ax: inout Vec<Double>) {
        fatalError("To be implemented")
    }
    
    override func setZero() {
        fatalError("To be implemented")
    }
    
    func setCoeff(rowI: Int, colI: Int, val: Double) {
        fatalError("To be implemented")
    }
    
    func addCoeff(rowI: Int, colI: Int, val: Double) {
        fatalError("To be implemented")
    }
}
