//
//  OptCutsError.swift
//  OptCuts
//
//  Created by Reza on 1/11/24.
//

import Foundation

public enum OptCutsError: Error {
    case elementInversion(Int)
    case nonManifoldEdge
    case nonManifoldVertex
}
