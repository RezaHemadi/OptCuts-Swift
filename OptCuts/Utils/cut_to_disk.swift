//
//  cut_to_disk.swift
//  OptCuts
//
//  Created by Reza on 6/26/23.
//

import Foundation
import Matrix

// Given a triangle mesh, computes a set of edge cuts sufficient to carve the
// mesh into a topological disk, without disconnecting any connected components.
// Nothing else about the cuts (including number, total length, or smoothness)
// is guaranteed to be optimal.
//
// Simply-connected components without boundary (topological spheres) are left
// untouched (delete any edge if you really want a disk).
// All other connected components are cut into disks. Meshes with boundary are
// supported; boundary edges will be included as cuts.
//
// The cut mesh itself can be materialized using cut_mesh().
//
// Implements the triangle-deletion approach described by Gu et al's
// "Geometry Images."
//
// Template Parameters:
//   Index  Integrable type large enough to represent the total number of faces
//     and edges in the surface represented by F, and all entries of F.
//
// Inputs:
//   F  #F by 3 list of the faces (must be triangles)
//
// Outputs:
//   cuts  List of cuts. Each cut is a sequence of vertex indices (where
//     pairs of consecutive vertices share a face), is simple, and is either
//     a closed loop (in which the first and last indices are identical) or
//     an open curve. Cuts are edge-disjoint.
//
 func cut_to_disk<MF: Matrix>(_ F: MF, _ cuts: inout [[Int]]) where MF.Element == Int {
    
    cuts.removeAll()
    let nfaces: Int = F.rows
    
    if (nfaces == 0) { return }
    
    var edges: [Pair<Int, Int> : [Int]] = [:]
    // build edges
    
    for i in 0..<nfaces {
        for j in 0..<3 {
            let v0: Int = F[i, j]
            let v1 = F[i, (j + 1) % 3]
            let e: Pair<Int, Int> = Pair(min(v0, v1), max(v0, v1))
            if edges[e] != nil {
                edges[e]!.append(i)
            } else {
                edges[e] = [i]
            }
        }
    }
    
    let nedges: Int = edges.count
    // TODO: Verify Eigen:Matrix<Index, -1, -1> is equivalent to Mat<Int>
    var edgeVerts = Mati(nedges, 2)
    var edgeFaces = Mati(nedges, 2)
    var faceEdges = Mati(nfaces, 3)
    var boundaryEdges = Set<Int>()
    var edgeidx = [Pair<Int, Int> : Int]()
    var idx: Int = 0
    for (key, value) in edges {
        edgeidx[key] = idx
        edgeVerts[idx, 0] = key.first
        edgeVerts[idx, 1] = key.second
        edgeFaces[idx, 0] = value[0]
        if (value.count > 1) {
            edgeFaces[idx, 1] = value[1]
        } else {
            edgeFaces[idx, 1] = -1
            boundaryEdges.insert(idx)
        }
        idx += 1
    }
    for i in 0..<nfaces {
        for j in 0..<3 {
            let v0: Int = F[i, j]
            let v1: Int = F[i, (j + 1) % 3]
            let e: Pair<Int, Int> = Pair(min(v0, v1), max(v0, v1))
            faceEdges[i, j] = edgeidx[e]!
        }
    }
    
    var deleted: [Bool] = .init(repeating: false, count: nfaces)
    var deletededges = Set<Int>()
    
    // loop over faces
    for face in 0..<nfaces {
        // stop at first undeleted face
        if (deleted[face]) {
            continue
        }
        deleted[face] = true
        var processEdges = [Int]()
        for i in 0..<3 {
            let e: Int = faceEdges[face, i]
            if (boundaryEdges.contains(e)) {
                continue
            }
            var ndeleted: Int = 0
            if (deleted[edgeFaces[e, 0]]) {
                ndeleted += 1
            }
            if (deleted[edgeFaces[e, 1]]) {
                ndeleted += 1
            }
            if (ndeleted == 1) {
                processEdges.append(e)
            }
        }
        // delete all faces adjacent to edges with exactly one adjacent face
        while (!processEdges.isEmpty) {
            let nexte: Int = processEdges.removeFirst()
            var todelete: Int = nfaces
            if (!deleted[edgeFaces[nexte, 0]]) {
                todelete = edgeFaces[nexte, 0]
            }
            if (!deleted[edgeFaces[nexte, 1]]) {
                todelete = edgeFaces[nexte, 1]
            }
            if (todelete != nfaces) {
                deletededges.insert(nexte)
                deleted[todelete] = true
                for i in 0..<3 {
                    let e: Int = faceEdges[todelete, i]
                    if (boundaryEdges.contains(e)) {
                        continue
                    }
                    var ndeleted: Int = 0
                    if (deleted[edgeFaces[e, 0]]) {
                        ndeleted += 1
                    }
                    if (deleted[edgeFaces[e, 1]]) {
                        ndeleted += 1
                    }
                    if (ndeleted == 1) {
                        processEdges.append(e)
                    }
                }
            }
        }
    }
    
    // accumulated non-deleted edges
    var leftedges: [Int] = []
    for i in 0..<nedges {
        if (!deletededges.contains(i)) {
            leftedges.append(i)
        }
    }
    
    deletededges.removeAll()
    
    // prune spines
    var spinevertedges: [Int : [Int]] = [:]
    for i in leftedges {
        if spinevertedges[edgeVerts[i, 0]] != nil {
            spinevertedges[edgeVerts[i, 0]]!.append(i)
        } else {
            spinevertedges[edgeVerts[i, 0]] = [i]
        }
        if spinevertedges[edgeVerts[i, 1]] != nil {
            spinevertedges[edgeVerts[i, 1]]!.append(i)
        } else {
            spinevertedges[edgeVerts[i, 1]] = [i]
        }
        
        var vertsProcess = [Int]()
        var spinevertnbs = [Int : Int]()
        
        for (key, value) in spinevertedges {
            spinevertnbs[key] = value.count
            if (value.count == 1) { vertsProcess.append(key) }
        }
        while (!vertsProcess.isEmpty) {
            let vert = vertsProcess.removeFirst()
            for e in spinevertedges[vert]! {
                if (!deletededges.contains(e)) {
                    deletededges.insert(e)
                    for j in 0..<2 {
                        spinevertnbs[edgeVerts[e, j]]! -= 1
                        if spinevertnbs[edgeVerts[e, j]]! == 1 {
                            vertsProcess.append(edgeVerts[e, j])
                        }
                    }
                }
            }
        }
        
        var loopedges = [Int]()
        for i in leftedges {
            if (!deletededges.contains(i)) {
                loopedges.append(i)
            }
        }
        
        let nloopedges = loopedges.count
        if (nloopedges == 0) { return }
        
        var loopvertedges = [Int : [Int]]()
        
        for e in loopedges {
            if loopvertedges[edgeVerts[e, 0]] != nil {
                loopvertedges[edgeVerts[e, 0]]!.append(e)
            } else {
                loopvertedges[edgeVerts[e, 0]] = [e]
            }
            if loopvertedges[edgeVerts[e, 1]] != nil {
                loopvertedges[edgeVerts[e, 1]]!.append(e)
            } else {
                loopvertedges[edgeVerts[e, 1]] = [e]
            }
        }
        
        var usededges = Set<Int>()
        
        for e in loopedges {
            // make a cycle or chain starting from this edge
            while (!usededges.contains(e)) {
                var cycleverts = [Int]()
                var cycleedges = [Int]()
                cycleverts.append(edgeVerts[e, 0])
                cycleverts.append(edgeVerts[e, 1])
                cycleedges.append(e)
                
                var cycleidx = [Int : Int]()
                cycleidx[cycleverts[0]] = 0
                cycleidx[cycleverts[1]] = 1
                
                var curvert: Int = edgeVerts[e, 1]
                var cure: Int = 0
                var foundcycle: Bool = false
                while (curvert != -1 && !foundcycle) {
                    var nextvert: Int = -1
                    var nexte: Int = -1
                    for cande in loopvertedges[curvert]! {
                        if (!usededges.contains(cande) && cande != cure) {
                            var vidx: Int = 0
                            if (curvert == edgeVerts[cande, vidx]) {
                                vidx = 1
                            }
                            nextvert = edgeVerts[cande, vidx]
                            nexte = cande
                            break
                        }
                    }
                    if (nextvert != -1) {
                        if let idx = cycleidx[nextvert] {
                            // we've hit outselves
                            var cut = [Int]()
                            for i in idx..<cycleverts.count {
                                cut.append(cycleverts[i])
                            }
                            cut.append(nextvert)
                            cuts.append(cut)
                            for i in idx..<cycleedges.count {
                                usededges.insert(cycleedges[i])
                            }
                            usededges.insert(nexte)
                            foundcycle = true
                        } else {
                            cycleidx[nextvert] = cycleverts.count
                            cycleverts.append(nextvert)
                            cycleedges.append(nexte)
                        }
                    }
                    curvert = nextvert
                    cure = nexte
                }
                if (!foundcycle) {
                    // we've hit a dead end. reverse and try the other direction
                    cycleverts.reverse()
                    cycleedges.reverse()
                    curvert = cycleverts.last!
                    cure = cycleedges.last!
                    
                    while (curvert != -1 && !foundcycle) {
                        var nextvert: Int = -1
                        var nexte: Int = -1
                        for cande in loopvertedges[curvert]! {
                            if (!usededges.contains(cande) && cande != cure) {
                                var vidx: Int = 0
                                if (curvert == edgeVerts[cande, vidx]) {
                                    vidx = 1
                                }
                                nextvert = edgeVerts[cande, vidx]
                                nexte = cande
                                break
                            }
                        }
                        if (nextvert != -1) {
                            if let idx = cycleidx[nextvert] {
                                // we've hit outseleves
                                var cut = [Int]()
                                for i in idx..<cycleverts.count {
                                    cut.append(cycleverts[i])
                                }
                                cut.append(nextvert)
                                cuts.append(cut)
                                
                                for i in idx..<cycleedges.count {
                                    usededges.insert(cycleedges[i])
                                }
                                usededges.insert(nexte)
                                foundcycle = true
                            } else {
                                cycleidx[nextvert] = cycleverts.count
                                cycleverts.append(nextvert)
                                cycleedges.append(nexte)
                            }
                        }
                        curvert = nextvert
                        cure = nexte
                    }
                    if (!foundcycle) {
                        // we've found a chain
                        var cut = [Int]()
                        for i in 0..<cycleverts.count {
                            cut.append(cycleverts[i])
                        }
                        cuts.append(cut)
                        for i in 0..<cycleedges.count {
                            usededges.insert(cycleedges[i])
                        }
                    }
                }
            }
        }
    }
}
