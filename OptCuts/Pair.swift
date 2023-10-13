//
//  Pair.swift
//  OptCuts
//
//  Created by Reza on 6/20/23.
//

import Foundation

struct Pair<T1: Hashable, T2: Hashable> {
    var first: T1
    var second: T2
    
    init(_ firstValue: T1, _ secondValue: T2) {
        self.first = firstValue
        self.second = secondValue
    }
}

extension Pair: Hashable {
    func hash(into hasher: inout Hasher) {
        hasher.combine(first)
        hasher.combine(second)
    }
}
