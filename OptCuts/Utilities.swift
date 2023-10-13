//
//  Utilities.swift
//  OptCuts
//
//  Created by Reza on 8/30/23.
//

import Foundation
import Matrix

protocol DefaultInitializable {
    init()
}

extension Set: DefaultInitializable where Element: BinaryInteger {}
extension Int: DefaultInitializable {}
extension Double: DefaultInitializable {}
extension Mat<Double>: DefaultInitializable {}
extension Pair: DefaultInitializable where T1: DefaultInitializable, T2: DefaultInitializable {
    init() {
        first = .init()
        second = .init()
    }
}

extension Array {
    mutating func conservativeResize(to newCount: Int) where Element: DefaultInitializable {
        if newCount == count { return }
        assert(newCount >= 0)
        
        if (newCount < count) {
            // remove extra elements from the end of array
            removeSubrange(startIndex+newCount..<endIndex)
        } else {
            // add extra elements to the array
            append(contentsOf: Array<Element>.init(repeating: .init(), count: newCount - count))
        }
    }
}
