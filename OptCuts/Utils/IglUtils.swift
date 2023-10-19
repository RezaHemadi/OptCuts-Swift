//
//  IglUtils.swift
//  OptCuts
//
//  Created by Reza on 6/27/23.
//

import Foundation
import Matrix

extension Matrix where Element: Equatable {
    func verifySymmetry() {
        for i in 0..<rows {
            for j in 0..<cols {
                guard i != j else { continue }
                
                assert(self[i, j] == self[j, i])
            }
        }
    }
}

// A static class implementing basic geometry processing operations that are not provided in libigl
 class IglUtils {
     static func computeGraphLaplacian<MF: Matrix>(_ F: MF,
                                                  _ graphL: inout SparseMatrix<Double>!) where MF.Element == Int {
        
        // compute vertex adjacency
        let vertAmt: Int = F.maxCoeff() + 1
        var adjacentVertices = [Set<Int>].init(repeating: Set<Int>(), count: vertAmt)
        
        for rowI in 0..<F.rows {
            adjacentVertices[F[rowI, 0]].insert(F[rowI, 1])
            adjacentVertices[F[rowI, 1]].insert(F[rowI, 0])
            adjacentVertices[F[rowI, 1]].insert(F[rowI, 2])
            adjacentVertices[F[rowI, 2]].insert(F[rowI, 1])
            adjacentVertices[F[rowI, 2]].insert(F[rowI, 0])
            adjacentVertices[F[rowI, 0]].insert(F[rowI, 2])
        }
        
        graphL = SparseMatrix<Double>.init(vertAmt, vertAmt)
        graphL.reserve(vertAmt * 7)
        var triplets: [Tripletd] = []
        triplets.reserveCapacity(vertAmt * 7)
        for rowI in 0..<vertAmt {
            triplets.append(.init(i: rowI, j: rowI, value: -Double(adjacentVertices[rowI].count)))
            for neighborI in adjacentVertices[rowI] {
                triplets.append(.init(i: rowI, j: neighborI, value: 1.0))
            }
        }
        graphL.setFromTriplets(triplets)
    }
    
    /// graph laplacian with half-weighted boundary edgee, the computation is also faster
     static func computeUniformLaplacian<MF: Matrix>(_ F: MF,
                                                    _ graphL: inout SparseMatrix<Double>!) where MF.Element == Int {
        
        let vertAmt: Int = F.maxCoeff() + 1
        var triplet: [Tripletd?] = .init(repeating: nil, count: F.rows * 9)
        // TODO: Make parallel
        for rowI in 0..<F.rows {
            let startInd = rowI * 9
            
            triplet[startInd] = .init(i: F[rowI, 0], j: F[rowI, 1], value: 1.0)
            triplet[startInd + 1] = .init(i: F[rowI, 1], j: F[rowI, 0], value: 1.0)
            triplet[startInd + 2] = .init(i: F[rowI, 1], j: F[rowI, 2], value: 1.0)
            triplet[startInd + 3] = .init(i: F[rowI, 2], j: F[rowI, 1], value: 1.0)
            triplet[startInd + 4] = .init(i: F[rowI, 2], j: F[rowI, 0], value: 1.0)
            triplet[startInd + 5] = .init(i: F[rowI, 0], j: F[rowI, 2], value: 1.0)
            
            triplet[startInd + 6] = .init(i: F[rowI, 0], j: F[rowI, 0], value: -2.0)
            triplet[startInd + 7] = .init(i: F[rowI, 1], j: F[rowI, 1], value: -2.0)
            triplet[startInd + 8] = .init(i: F[rowI, 2], j: F[rowI, 2], value: -2.0)
        }
        graphL = SparseMatrix<Double>(vertAmt, vertAmt)
        graphL.setFromTriplets(triplet.compactMap({$0}))
    }
    
     static func computeMVCMtr<MV: Matrix, MF: Matrix>
    (_ V: MV, _ F: MF, _ MVCMtr: inout SparseMatrix<Double>!)
    where MV.Element == Double, MF.Element == Int {
        
        var HETan: [Pair<Int, Int> : Double] = [:]
        var thirdPoint: [Pair<Int, Int> : Int] = [:]
        var vvNeighbor: [Set<Int>] = .init(repeating: Set<Int>(), count: V.rows)
        
        for triI in 0..<F.rows {
            let v0I: Int = F[triI, 0]
            let v1I: Int = F[triI, 1]
            let v2I: Int = F[triI, 2]
            
            let e01: Vec3d = V.row(v1I) - V.row(v0I)
            let e12: Vec3d = V.row(v2I) - V.row(v1I)
            let e20: Vec3d = V.row(v0I) - V.row(v2I)
            let dot0102: Double = -e01.dot(e20)
            let dot1210: Double = -e12.dot(e01)
            let dot2021: Double = -e20.dot(e12)
            let cos0102: Double = dot0102 / (e01.norm() * e20.norm())
            let cos1210: Double = dot1210 / (e01.norm() * e12.norm())
            let cos2021: Double = dot2021 / (e12.norm() * e20.norm())
            
            HETan[Pair(v0I, v1I)] = sqrt(1.0 - cos0102 * cos0102) / (1.0 + cos0102)
            HETan[Pair(v1I, v2I)] = sqrt(1.0 - cos1210 * cos1210) / (1.0 + cos1210)
            HETan[Pair(v2I, v0I)] = sqrt(1.0 - cos2021 * cos2021) / (1.0 + cos2021)
            
            thirdPoint[Pair(v0I, v1I)] = v2I
            thirdPoint[Pair(v1I, v2I)] = v0I
            thirdPoint[Pair(v2I, v0I)] = v1I
            
            vvNeighbor[v0I].insert(v1I)
            vvNeighbor[v0I].insert(v2I)
            vvNeighbor[v1I].insert(v0I)
            vvNeighbor[v1I].insert(v2I)
            vvNeighbor[v2I].insert(v0I)
            vvNeighbor[v2I].insert(v1I)
        }
        
        MVCMtr = SparseMatrix<Double>.init(V.rows, V.rows)
        MVCMtr.reserve(V.rows * 7)
        var triplets: [Tripletd] = []
        triplets.reserveCapacity(V.rows * 7)
        for rowI in 0..<V.rows {
            for nbVI in vvNeighbor[rowI] {
                var weight: Double = getHETan(HETan, rowI, nbVI)
                if let point = thirdPoint[Pair(nbVI, rowI)] {
                    weight += getHETan(HETan, rowI, point)
                }
                weight /= (V.row(rowI) - V.row(nbVI)).norm()
                
                triplets.append(.init(i: rowI, j: rowI, value: -weight))
                triplets.append(.init(i: rowI, j: nbVI, value: weight))
            }
        }
        MVCMtr.setFromTriplets(triplets)
    }
    
     static func fixedBoundaryParam_MVC<V: Vector, M1: Matrix, M2: Matrix>
    (_ A: SparseMatrix<Double>, _ bnd: V, bnd_uv: M1, _ UV_Tutte: inout M2!)
    where V.Element == Int, M1.Element == Double, M2.Element == Double, M1.ColType: Vector {
        
        var A = A
        assert(bnd.count == bnd_uv.rows)
        assert(bnd.maxCoeff() < A.rows)
        assert(A.rows == A.cols)
        
        let vN: Int = A.rows
        A.conservativeResize(vN + bnd.count, vN + bnd.count)
        A.reserve(A.nonZeros + bnd.count * 2)
        for pcI in 0..<bnd.count {
            A.insert(row: vN + pcI, col: bnd[pcI], value: 1.0)
            A.insert(row: bnd[pcI], col: vN + pcI, value: 1.0)
        }
        
        let spLUSolver = SparseLLT(A, transpose: false)
        UV_Tutte.resize(A.rows, 2)
        var rhs = Vec<Double>(A.rows)
        
        for dimI in 0..<2 {
            rhs <<== [Vec<Double>.Zero(vN), bnd_uv.col(dimI)]
            UV_Tutte.col(dimI) <<== spLUSolver.solve(b: rhs)
        }
        
        UV_Tutte.conservativeResize(vN, 2)
    }
    
     static func mapTriangleTo2D(_ v: [Vec3d], _ u: inout [Vec2d]) {
        let e: [Vec3d] = [
            v[1] - v[0],
            v[2] - v[0]
        ]
        
        u[0] <<== [0.0, 0.0]
        u[1] <<== [e[0].norm(), 0.0]
        u[2] <<== [e[0].dot(e[1]) / u[1][0], e[0].cross(e[1]).norm() / u[1][0]]
    }
    
     static func computeDeformationGradient(_ v: [Vec3d], _ u: [Vec2d], _ F: inout Mat2d) {
        var x: [Vec2d] = .init(repeating: Vec2d(), count: 3)
        IglUtils.mapTriangleTo2D(v, &x)
        
        let u01: Vec2d = u[1] - u[0]
        let u02: Vec2d = u[2] - u[0]
        let u01Len: Double = u01.norm()
        let U = Mat2d([u01Len, u01.dot(u02) / u01Len, 0.0, (u01[0] * u02[1] - u01[1] * u02[0]) / u01Len], [2, 2])
        var V = Mat2d()
        V <<== [x[1], x[2]]
        F = V * U.inverse()
    }
    
    // to a circle with the perimeter equal to the length of the boundary on the mesh
     static func map_vertices_to_circle<MUV: Matrix>(_ V: Matd,
                                                    _ bnd: Veci,
                                                    _ UV: inout MUV!) where MUV.Element == Double {
        // Get sorted list of boundary vertices
        var interior = [Int]()
        var map_ij = [Int].init(repeating: -1, count: V.rows)
        
        var isOnBnd: [Bool] = .init(repeating: false, count: V.rows)
        
        for i in 0..<bnd.count {
            isOnBnd[bnd[i]] = true
            map_ij[bnd[i]] = i
        }
        
        for i in 0..<isOnBnd.count {
            if (!isOnBnd[i]) {
                map_ij[i] = interior.count
                interior.append(i)
            }
        }
        
        // Map boundary to circle
        var len = [Double].init(repeating: 0.0, count: bnd.count)
        
        for i in 1..<bnd.count {
            len[i] = len[i - 1] + (V.row(bnd[i - 1]) - V.row(bnd[i])).norm()
        }
        let total_len: Double = len[len.count - 1] + (V.row(bnd[0]) - V.row(bnd[bnd.count - 1])).norm()
        
        UV = .init(bnd.count, 2)
        let radius: Double = total_len / 2.0 / .pi
        for i in 0..<bnd.count {
            let frac: Double = len[i] * 2.0 * .pi / total_len
            UV.row(map_ij[bnd[i]]) <<== [radius * cos(frac), radius * sin(frac)]
        }
    }
    
     static func mapScalarToColor_bin<V: Vector, M: Matrix>(_ scalar: V, _ color: inout M, _ thres: Double)
    where V.Element == Double, M.Element == Double {
        fatalError("To be implemented")
    }
    
     static func mapScalarToColor<V: Vector, M: Matrix>(_ scalar: V, _ color: M, _ lowerBound: Double, upperBound: Double, _ opt: Int = 0) where V.Element == Double, M.Element == Double {
        fatalError("To be implemented")
    }
    
     static func addBlockToMatrix<M: Matrix, V: Vector, S: MatrixElement>(_ mtr: inout SparseMatrix<S>, _ block: M, _ index: V, _ dim: Int) where M.Element == S, V.Element == Int {
        
        assert(block.rows == block.cols)
        assert(index.count * dim == block.rows)
        assert(mtr.rows == mtr.cols)
        assert(index.maxCoeff() * dim + dim - 1 < mtr.rows)
        
        for indI in 0..<index.count {
            if (index[indI] < 0) {
                continue
            }
            
            let startIndI: Int = index[indI] * dim
            let startIndI_block: Int = indI * dim
            
            for indJ in 0..<index.count {
                if (index[indJ] < 0) {
                    continue
                }
                let startIndJ = index[indJ] * dim
                let startIndJ_block = indJ * dim
                
                for dimI in 0..<dim {
                    for dimJ in 0..<dim {
                        mtr.coeffRef(row: startIndI + dim, col: startIndJ + dim) { value in
                            value + block[startIndI_block + dimI, startIndJ_block + dimJ]
                        }
                    }
                }
            }
        }
    }
    
     static func addBlockToMatrix<M: Matrix, V1: Vector, V2: Vector, V3: Vector, V4: Vector>
    (_ block: M, _ index: V1, _ dim: Int, _ V: inout V2, _ I: inout V3, _ J: inout V4)
    where M.Element == Double, V1.Element == Int, V2.Element == Double, V3.Element == Int, V4.Element == Int {
        
        var num_free: Int = 0
        for indI in 0..<index.count {
            if (index[indI] >= 0) {
                num_free += 1
            }
        }
        if (num_free == 0) {
            return
        }
        
        assert(block.rows == block.cols)
        assert(index.count * dim == block.rows)
        
        var tripletInd: Int = V.count
        let entryAmt: Int = dim * dim * num_free * num_free
        V.conservativeResize(tripletInd + entryAmt)
        I.conservativeResize(tripletInd + entryAmt)
        J.conservativeResize(tripletInd + entryAmt)
        
        for indI in 0..<index.count {
            if (index[indI] < 0) {
                continue
            }
            let startIndI: Int = index[indI] * dim
            let startIndI_block: Int = indI * dim
            
            for indJ in 0..<index.count {
                if (index[indJ] < 0) {
                    continue
                }
                let startIndJ: Int = index[indJ] * dim
                let startIndJ_block: Int = indJ * dim
                
                for dimI in 0..<dim {
                    for dimJ in 0..<dim {
                        V[tripletInd] = block[startIndI_block + dimI, startIndJ_block + dimJ]
                        I[tripletInd] = startIndI + dimI
                        J[tripletInd] = startIndJ + dimJ
                        tripletInd += 1
                    }
                }
            }
        }
        assert(tripletInd == V.count)
    }
    
     static func addDiagonalToMatrix<V1: Vector, V2: Vector, S: MatrixElement, V3: Vector, V4: Vector, V5: Vector>
    (_ diagonal: V1, _ index: V2, _ dim: Int, _ V: inout V3, _ I: inout V4, _ J: inout V5)
    where V1.Element == S, V2.Element == Int, V3.Element == S, V4.Element == Int, V5.Element == Int {
        
        assert(index.count * dim == diagonal.count)
        
        var tripletInd: Int = V.count
        let entryAmt: Int = diagonal.count
        V.conservativeResize(tripletInd + entryAmt)
        
        I.conservativeResize(tripletInd + entryAmt)
        J.conservativeResize(tripletInd + entryAmt)
        
        for indI in 0..<index.count {
            if (index[indI] < 0) {
                assert(false, "currently doensn't support fixed vertices here!")
                continue
            }
            let startIndI: Int = index[indI] * dim
            let startIndI_diagoanl: Int = indI * dim
            
            for dimI in 0..<dim {
                V[tripletInd] = diagonal[startIndI_diagoanl + dimI]
                J[tripletInd] = startIndI + dimI
                I[tripletInd] = startIndI + dimI
                tripletInd += 1
            }
        }
    }
    
     static func addBlockToMatrix<M1: Matrix, V: Vector, M2: Matrix>
    (_ block: M1, _ index: V, _ dim: Int, _ mtr: inout M2)
    where M1.Element == M2.Element, V.Element == Int, M1.Element: AdditiveArithmetic {
        
        var num_free: Int = 0
        for indI in 0..<index.count {
            if (index[indI] >= 0) {
                num_free += 1
            }
        }
        if (num_free == 0) {
            return
        }
        
        assert(block.rows == block.cols)
        assert(index.count * dim == block.rows)
        assert(mtr.rows == mtr.cols)
        assert(mtr.rows % dim == 0)
        
        for indI in 0..<index.count {
            if (index[indI] < 0) {
                continue
            }
            let startIndI: Int = index[indI] * dim
            let startIndI_block: Int = indI * dim
            
            for indJ in 0..<index.count {
                if (index[indJ] < 0) {
                    continue
                }
                let startIndJ: Int = index[indJ] * dim
                let startIndJ_block: Int = indJ * dim
                
                mtr.block(startIndI, startIndJ, dim, dim) +=
                    block.block(startIndI_block, startIndJ_block, dim, dim)
            }
        }
    }
    
     static func addDiagonalToMatrix<V1: Vector, V2: Vector, M1: Matrix>
    (_ diagonal: V1, _ index: V2, _ dim: Int, _ mtr: inout M1)
    where V1.Element == M1.Element, V2.Element == Int {
        
        assert(index.count * dim == diagonal.count)
        assert(mtr.rows == mtr.cols)
        assert(mtr.rows % dim == 0)
        
        for indI in 0..<index.count {
            if (index[indI] < 0) {
                assert(false, "currently doesn't support fixed vertices here!")
                continue
            }
            let startIndI: Int = index[indI] * dim
            let startIndI_diagonal: Int = indI * dim
            
            mtr[startIndI, startIndI] = diagonal[startIndI_diagonal]
            mtr[startIndI + 1, startIndI + 1] = diagonal[startIndI_diagonal + 1]
        }
    }
    
     static func symmetrizeMatrix(_ mtr: inout any Matrix)  {
        fatalError("To be implemented")
    }
    
    // project a symmetric real matrix to the nearest SPD matrix
     static func makePD<M: Matrix>(_ symMtr: inout M) throws where M.Element == Double {
        /*
        let solver = try EigenSolver(squareMatrix: symMtr)
        if (solver.realEigenValues[0] >= 0.0) {
            return
        }
        var D: Mat<Double> = solver.realEigenValues.asDiagonal()
        let rows: Int = symMtr.rows
        
        for i in 0..<rows {
            if (D[i, i] < 0.0) {
                D[i, i] = 0.0
            } else {
                break
            }
        }
        symMtr = solver.eigenVectors * D * solver.eigenVectors.transpose()*/
        
        let solver = try SelfAdjointEigenSolver(squareMatrix: symMtr)
         
        if (solver.eigenValues[0] >= 0.0) {
            return
        }
        var D: Mat<Double> = solver.eigenValues.asDiagonal()
        let rows: Int = symMtr.rows
        
        for i in 0..<rows {
            if (D[i, i] < 0.0) {
                D[i, i] = 0.0
            } else {
                break
            }
        }
        
        symMtr = solver.eigenVectors * D * solver.eigenVectors.transpose()
    }
    
     static func writeSparseMatrixToFile(_ filePath: String, _ mtr: SparseMatrix<Double>, _ MATLAB: Bool = false) {
        fatalError("To be implemented")
    }
    
     static func writeSparseMatrixToFile<V1: Vector, V2: Vector, V3: Vector>
    (_ filePath: String, _ I: V1, _ J: V2, _ V: V3, _ MATLAB: Bool = false)
    where V1.Element == Int, V2.Element == Int, V3.Element == Double {
        fatalError("To be implemented")
    }
    
     static func loadSparseMatrixFromFile<S: MatrixElement>(_ filePath: String, _ mtr: inout SparseMatrix<S>) {
        fatalError("To be implemented")
    }
    
     static func sparseMatrixToTriplet<S: MatrixElement, V1: Vector, V2: Vector, V3: Vector>
    (_ mtr: SparseMatrix<S>, _ I: inout V1, _ J: inout V2, _ V: inout V3)
    where V1.Element == Int, V2.Element == Int, V3.Element == S {
        
        I.resize(mtr.nonZeros)
        J.resize(mtr.nonZeros)
        V.resize(mtr.nonZeros)
        
        var entryI: Int = 0
        for k in 0..<mtr.outerSize {
            for it in mtr.innerIterator(k) {
                I[entryI] = it.row
                J[entryI] = it.col
                V[entryI] = it.value
                entryI += 1
            }
        }
    }
    
     static func sparseMatrixToTriplet<S: MatrixElement, V1: Vector>(_ mtr: SparseMatrix<S>, _ V: inout V1)
    where V1.Element == S {
        
        V.resize(mtr.nonZeros)
        var entryI: Int = 0
        for k in 0..<mtr.outerSize {
            for it in mtr.innerIterator(k) {
                V[entryI] = it.value
                entryI += 1
            }
        }
    }
    
     static func rtos(_ real: Double) -> String {
        return String(real)
    }
    
     static func computeRotAngle(_ from: RVec2d, _ to: RVec2d) -> Double {
        let angle: Double = acos(max(-1.0, min(1.0, from.dot(to) / from.norm() / to.norm())))
        return ((from[0] * to[1] - from[1] * to[0] < 0.0) ? -angle : angle)
    }
    
    // test whether 2D segments ab intersect with cd
     static func Test2DsegmentSegment(_ a: RVec2d, _ b: RVec2d, _ c: RVec2d, _ d: RVec2d, _ eps: Double = 0.0) -> Bool {
        
        var eps: Double = eps
        var eps_quad: Double = 0.0
        var eps_sq: Double = 0.0
        
        if (eps != 0) {
            eps = abs(eps)
            eps_sq = eps * eps * ((a - b).squaredNorm() + (c - d).squaredNorm()) / 2.0
            eps_quad = eps_sq * eps_sq
        }
        
        // signs of areas correspond to which side of ab points c and b are
        let a1: Double = Signed2DTriArea(a, b, d) // Compute winding of abd (+ or -)
        let a2: Double = Signed2DTriArea(a, b, c) // to intersect, must have sign opposite of a1
        
        // If c and d are on different sides of ab, areas have different signes
        if (a1 * a2 <= eps_quad) { // require unsigned x & y values
            let a3 = Signed2DTriArea(c, d, a) // Compute winding of cda(+ or -)
            let a4 = a3 + a2 - a1 // Since area is constant a1 - a2 = a3 - a4, or a4 = a3 + a2 - a1
            
            // Points a and b on different sides of cd if areas have different signs
            if (a3 * a4 <= eps_quad) {
                if ((abs(a1) <= eps_sq) && (abs(a2) <= eps_sq)) {
                    // colinear
                    let ab: RVec2d = b - a
                    let sqnorm_ab: Double = ab.squaredNorm()
                    let ac: RVec2d = c - a
                    let ad: RVec2d = d - a
                    var coef_c: Double = ac.dot(ab) / sqnorm_ab
                    var coef_d: Double = ad.dot(ab) / sqnorm_ab
                    assert(coef_c != coef_d)
                    
                    if (coef_c > coef_d) {
                        swap(&coef_c, &coef_d)
                    }
                    
                    if ((coef_c > 1.0 + eps) || (coef_d < -eps)) {
                        return false
                    } else {
                        return true
                    }
                } else {
                    // Segments intersect
                    return true
                }
            }
        }
        
        // Segments not intersecting
        return false
    }
    /*
     static func addThickEdge<M1: Matrix, M2: Matrix, M3: Matrix, M4: Matrix>
    (_ V: inout M1, _ F: inout M2, _ UV: inout M3, _ seamColor: inout M4, _ color: RVec3d, _ v0: RVec3d, _ v1: RVec3d, _ halfWidth: Double, _ textScale: Double, _ UVorSurface: Bool = false, _ normal: RVec3d = RVec3d())
    where M1.Element == Double, M2.Element == Int, M3.Element == Double, M4.Element == Double {
        
        if (UVorSurface) {
            let e: RVec3d = v1 - v0
            let n: RVec3d = normal.normalized() * halfWidth
            let bn: RVec3d = e.cross(normal).normalized() * halfWidth
            let vAmt_old: Int = V.rows
            V.conservativeResize(V.rows + 8, 3)
            let VVecs: [RVec3d] = [v0 - n - bn, v0 - n + bn, v0 + n + bn, v0 + n - bn, v1 - n - bn, v1 - n + bn, v1 + n + bn, v1 + n - bn]
            V.bottomRows(8) <<== VVecs
            UV.conservativeResize(UV.rows + 8, 2)
            UV.bottomRows(8) <<== Matd.Ones(8, 2) * textScale
            F.conservativeResize(F.rows + 6, 3)
            F.bottomRows(6) <<== [vAmt_old + 2, vAmt_old + 1, vAmt_old + 5,
                vAmt_old + 2, vAmt_old + 5, vAmt_old + 6,
                vAmt_old + 3, vAmt_old + 2, vAmt_old + 6,
                vAmt_old + 3, vAmt_old + 6, vAmt_old + 7,
                vAmt_old, vAmt_old + 3, vAmt_old + 7,
                vAmt_old, vAmt_old + 7, vAmt_old + 4]
            seamColor.conservativeResize(seamColor.rows + 6, 3)
            seamColor.bottomRows(6) <<== [color, color, color, color, color, color]
        } else {
            let e: RVec3d = v1 - v0
            let n: RVec3d = halfWidth * RVec3d([-e[1], e[0], 0.0]).normalized()
            let vAmt_old: Int = V.rows
            V.conservativeResize(V.rows + 4, 3)
            let VVecs: [RVec3d] = [v0 - n, v0 + n, v1 + n, v1 - n]
            V.bottomRows(4) <<== VVecs
            V.bottomRows(4).col(2) <<== Vec<Double>.Ones(4) * halfWidth // for depth test
            F.conservativeResize(F.rows + 2, 3)
            F.bottomRows(2) <<== [vAmt_old, vAmt_old + 1, vAmt_old + 2,
                vAmt_old, vAmt_old + 2, vAmt_old + 3]
            seamColor.conservativeResize(seamColor.rows + 2, 3)
            seamColor.bottomRows(2) <<== [color, color]
        }
    }*/
    
     static func saveMesh_Seamster<MV: Matrix, MF: Matrix>
    (_ filePath: String, _ V: MV, _ F: MF)
    where MV.Element == Double, MF.Element == Int {
        fatalError("To be implemented")
    }
    
     static func smoothVertField<V: Vector>(_ mesh: TriMesh, _ field: inout V) where V.Element == Double {
        
        assert(field.count == mesh.V.rows)
        let field_copy: Vec<Double> = .init(field, field.rows, field.cols)
        for vI in 0..<field.count {
            for nbVI in mesh.vNeighbor[vI] {
                field[vI] += field_copy[nbVI]
            }
            field[vI] /= Double(mesh.vNeighbor[vI].count + 1)
        }
    }
}

// MARK: Helper Methods
private func getHETan(_ HETan: [Pair<Int, Int> : Double],
                      _ v0: Int,
                      _ v1: Int) -> Double {
    
    if let value = HETan[Pair(v0, v1)] {
        return value
    } else {
        return 0.0
    }
}

////////////////////////////////////////////////////////////
// 2D line segments intersection checking code
// based on Real-Time Collision Detection by Christer Ericson
// (Morgan Kaufmaan Publishers, 2005 Elvesier Inc)
private func Signed2DTriArea(_ a: RVec2d, _ b: RVec2d, _ c: RVec2d) -> Double {
    return (a[0] - c[0]) * (b[1] - c[1]) - (a[1] - c[1]) * (b[0] - c[0])
}
