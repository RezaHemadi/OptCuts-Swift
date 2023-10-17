//
//  SymDirichletEnergy.swift
//  OptCuts
//
//  Created by Reza on 6/20/23.
//

import Foundation
import Matrix

 class SymDirichletEnergy: Energy {
    // MARK: - Initialization
     init() {
        super.init(p_needRefactorize: true)
    }
    
    // MARK: - Methods
    override  func getEnergyValPerElem(_ data: TriMesh,
                                      _ energyValPerElem: inout Vec<Double>,
                                      _ uniformWeight: Bool = false) {
        let normalizer_div: Double = data.surfaceArea
        
        energyValPerElem.resize(data.F.rows)
        for triI in 0..<data.F.rows {
            let triVInd: Vec3i = data.F.row(triI)
            
            let U1: RVec2d = data.V.row(triVInd[0])
            let U2: RVec2d = data.V.row(triVInd[1])
            let U3: RVec2d = data.V.row(triVInd[2])
            
            let U2m1: RVec2d = U2 - U1
            let U3m1: RVec2d = U3 - U1
            
            let area_U: Double = 0.5 * (U2m1[0] * U3m1[1] - U2m1[1] * U3m1[0])
            
            let w: Double = (uniformWeight ? 1.0 : (data.triArea[triI] / normalizer_div))
            energyValPerElem[triI] = w * (1.0 + data.triAreaSq[triI] / area_U / area_U) *
            ((U3m1.squaredNorm() * data.e0SqLen[triI] + U2m1.squaredNorm() * data.e1SqLen[triI]) / 4.0 / data.triAreaSq[triI] - U3m1.dot(U2m1) * data.e0dote1[triI] / 2.0 / data.triAreaSq[triI])
        }
    }
    
    override  func getEnergyValByElemID(_ data: TriMesh,
                                              _ elemI: Int,
                                              _ energyVal: inout Double,
                                              _ uniformWeight: Bool = false) {
        
        let normalizer_div: Double = data.surfaceArea
        
        let triI: Int = elemI
        let triVInd: Vec3i = data.F.row(triI)
        
        let U1: Vec2d = data.V.row(triVInd[0])
        let U2: Vec2d = data.V.row(triVInd[1])
        let U3: Vec2d = data.V.row(triVInd[2])
        
        let U2m1: Vec2d = U2 - U1
        let U3m1: Vec2d = U3 - U1
        
        let area_U: Double = 0.5 * (U2m1[0] * U3m1[1] - U2m1[1] * U3m1[0])
        
        let w: Double = (uniformWeight ? 1.0 : (data.triArea[triI] / normalizer_div))
        energyVal = w * (1.0 + data.triAreaSq[triI] / area_U / area_U) * ((U3m1.squaredNorm() * data.e0SqLen[triI] + U2m1.squaredNorm() * data.e1SqLen[triI]) / 4.0 / data.triAreaSq[triI] - U3m1.dot(U2m1) * data.e0dote1[triI] / 2.0 / data.triAreaSq[triI])
    }
    
    override  func computeGradient(_ data: TriMesh,
                                  _ gradient: inout Vec<Double>,
                                  _ uniformWeight: Bool = false) {
        
        let normalizer_div: Double = data.surfaceArea
        
        gradient.resize(data.V.rows * 2)
        gradient.setZero()
        
        for triI in 0..<data.F.rows {
            let triVInd: Vec3i = data.F.row(triI)
            
            let U1: Vec2d = data.V.row(triVInd[0])
            let U2: Vec2d = data.V.row(triVInd[1])
            let U3: Vec2d = data.V.row(triVInd[2])
            
            let U2m1: Vec2d = U2 - U1
            let U3m1: Vec2d = U3 - U1
            
            let area_U: Double = 0.5 * (U2m1[0] * U3m1[1] - U2m1[1] * U3m1[0])
            
            let leftTerm: Double = 1.0 + data.triAreaSq[triI] / area_U / area_U
            let rightTerm: Double = (U3m1.squaredNorm() * data.e0SqLen[triI] + U2m1.squaredNorm() * data.e1SqLen[triI]) / 4.0 / data.triAreaSq[triI] - U3m1.dot(U2m1) * data.e0dote1[triI] / 2.0 / data.triAreaSq[triI]
            
            let areaRatio: Double = data.triAreaSq[triI] / area_U / area_U / area_U
            let w: Double = (uniformWeight ? 1.0 : (data.triArea[triI] / normalizer_div))
            
            let edge_oppo1: Vec2d = U3 - U2
            let dLeft1: Vec2d = areaRatio * Vec2d(edge_oppo1[1], -edge_oppo1[0])
            let dRight1: Vec2d = ((data.e0dote1[triI] - data.e0SqLen[triI]) * U3m1 + (data.e0dote1[triI] - data.e1SqLen[triI]) * U2m1) / 2.0 / data.triAreaSq[triI]
            gradient.block(triVInd[0] * 2, 0, 2, 1) += w * (dLeft1 * rightTerm + dRight1 * leftTerm)
            
            let edge_oppo2: Vec2d = U1 - U3
            let dLeft2: Vec2d = areaRatio * Vec2d(edge_oppo2[1], -edge_oppo2[0])
            let dRight2: Vec2d = (data.e1SqLen[triI] * U2m1 - data.e0dote1[triI] * U3m1) / 2.0 / data.triAreaSq[triI]
            gradient.block(triVInd[1] * 2, 0, 2, 1) += w * (dLeft2 * rightTerm + dRight2 * leftTerm)
            
            let edge_oppo3: Vec2d = U2 - U1
            let dLeft3: Vec2d = areaRatio * Vec2d(edge_oppo3[1], -edge_oppo3[0])
            let dRight3: Vec2d = (data.e0SqLen[triI] * U3m1 - data.e0dote1[triI] * U2m1) / 2.0 / data.triAreaSq[triI]
            gradient.block(triVInd[2] * 2, 0, 2, 1) += w * (dLeft3 * rightTerm + dRight3 * leftTerm)
        }
        
        for fixedVI in data.fixedVert {
            gradient[2 * fixedVI] = 0.0
            gradient[2 * fixedVI + 1] = 0.0
        }
    }
    
    override  func computeHessian(_ data: TriMesh,
                                 _ V: inout Vec<Double>,
                                 _ I: inout Vec<Int>,
                                 _ J: inout Vec<Int>,
                                 _ uniformWeight: Bool = false) {
        
        let normalizer_div = data.surfaceArea
        
        //var triHessians: [Mat6d] = .init(repeating: Mat6d(), count: data.F.rows)
        let triHessians: UnsafeMutablePointer<Mat6d> = .allocate(capacity: data.F.rows)
        //var vInds: [Vec3i] = .init(repeating: Vec3i(), count: data.F.rows)
        let vInds: UnsafeMutablePointer<Vec3i> = .allocate(capacity: data.F.rows)
        defer {
            triHessians.deallocate()
            vInds.deallocate()
        }
        
        DispatchQueue.concurrentPerform(iterations: data.F.rows) { triI in
            (triHessians + triI).initialize(to: Mat6d())
            (vInds + triI).initialize(to: Vec3i())
            let curHessian = triHessians + triI
            let mat2dIdentity = Mat2d.Identity()
            let triAreaSq = data.triAreaSq[triI]
            
            let triVInd: Vec3i = data.F.row(triI)
            let U1: Vec2d = data.V.row(triVInd[0])
            let U2: Vec2d = data.V.row(triVInd[1])
            let U3: Vec2d = data.V.row(triVInd[2])
            
            let U2m1: Vec2d = U2 - U1
            let U3m1: Vec2d = U3 - U1
            
            let area_U = 0.5 * (U2m1[0] * U3m1[1] - U2m1[1] * U3m1[0])
            let areaRatio = triAreaSq / area_U / area_U / area_U
            let dAreaRatio_div_dArea_mult: Double = 3.0 / 2.0 * areaRatio / area_U
            
            let w: Double = (uniformWeight ? 1.0 : (data.triArea[triI] / normalizer_div))
            
            let e0SqLen_div_dbAreaSq: Double = data.e0SqLen_div_dbAreaSq[triI]
            let e1SqLen_div_dbAreaSq: Double = data.e1SqLen_div_dbAreaSq[triI]
            let e0dote1_div_dbAreaSq: Double = data.e0dote1_div_dbAreaSq[triI]
            
            // compute energy terms
            let leftTerm: Double = 1.0 + triAreaSq / area_U / area_U
            let rightTerm: Double = (U3m1.squaredNorm() * e0SqLen_div_dbAreaSq + U2m1.squaredNorm() * e1SqLen_div_dbAreaSq) / 2.0 - U3m1.dot(U2m1) * e0dote1_div_dbAreaSq
            
            let edge_oppo1: Vec2d = U3 - U2
            let edge_oppo2: Vec2d = U1 - U3
            let edge_oppo3: Vec2d = U2 - U1
            let edge_oppo1_Ortho = Vec2d(edge_oppo1[1], -edge_oppo1[0])
            let edge_oppo2_Ortho = Vec2d(edge_oppo2[1], -edge_oppo2[0])
            let edge_oppo3_Ortho = Vec2d(edge_oppo3[1], -edge_oppo3[0])
            let dOrtho_div_dU: Mat2d = .init([0.0, -1.0, 1.0, 0.0], [2, 2])
            
            // compute 1st order derivatives
            let dLeft1: Vec2d = areaRatio * edge_oppo1_Ortho
            let dRight1: Vec2d = ((e0dote1_div_dbAreaSq - e0SqLen_div_dbAreaSq) * U3m1 +
                                  (e0dote1_div_dbAreaSq - e1SqLen_div_dbAreaSq) * U2m1)
            
            let dLeft2: Vec2d = areaRatio * edge_oppo2_Ortho
            let dRight2: Vec2d = (e1SqLen_div_dbAreaSq * U2m1 - e0dote1_div_dbAreaSq * U3m1)
            
            let dLeft3: Vec2d = areaRatio * edge_oppo3_Ortho
            let dRight3 = (e0SqLen_div_dbAreaSq * U3m1 - e0dote1_div_dbAreaSq * U2m1)
            
            // compute second order derivatives for g_U1
            let d2Left11: Mat2d = dAreaRatio_div_dArea_mult * edge_oppo1_Ortho * edge_oppo1_Ortho.transpose()
            let d2Right11: Double = (e0SqLen_div_dbAreaSq + e1SqLen_div_dbAreaSq - 2.0 * e0dote1_div_dbAreaSq)
            let dLeft1dRight1T: Mat2d = dLeft1 * dRight1.transpose()
            let exp1: Mat2d = d2Left11 * rightTerm + dLeft1dRight1T
            let exp2: Mat2d = mat2dIdentity * d2Right11 * leftTerm
            let exp3: Mat2d = dLeft1dRight1T.transpose()
            curHessian.pointee.block(0, 0, 2, 2) <<== w * (exp1 + exp2 + exp3)
            let d2Left12: Mat2d = dAreaRatio_div_dArea_mult * edge_oppo1_Ortho * edge_oppo2_Ortho.transpose() + areaRatio * dOrtho_div_dU
            let d2Right12: Double = (e0dote1_div_dbAreaSq - e1SqLen_div_dbAreaSq)
            let exp13: Mat2d = d2Left12 * rightTerm + dLeft1 * dRight2.transpose()
            let exp14: Mat2d = mat2dIdentity * d2Right12 * leftTerm
            let exp15: Mat2d = dRight1 * dLeft2.transpose()
            let result1: Mat2d = w * (exp13 + exp14 + exp15)
            curHessian.pointee.block(0, 2, 2, 2) <<== result1
            curHessian.pointee.block(2, 0, 2, 2) <<== result1.transpose()
            
            let d2Left13: Mat2d = dAreaRatio_div_dArea_mult * edge_oppo1_Ortho * edge_oppo3_Ortho.transpose() +
                                  areaRatio * (-dOrtho_div_dU)
            let d2Right13: Double = (e0dote1_div_dbAreaSq - e0SqLen_div_dbAreaSq)
            let exp16: Mat2d = d2Left13 * rightTerm + dLeft1 * dRight3.transpose()
            let exp17: Mat2d = mat2dIdentity * d2Right13 * leftTerm
            let exp18: Mat2d = dRight1 * dLeft3.transpose()
            let result2: Mat2d = w * (exp16 + exp17 + exp18)
            curHessian.pointee.block(0, 4, 2, 2) <<== result2
            curHessian.pointee.block(4, 0, 2, 2) <<== result2.transpose()
            
            // compute second order derivatives for g_U2
            let d2Left22: Mat2d = dAreaRatio_div_dArea_mult * edge_oppo2_Ortho * edge_oppo2_Ortho.transpose()
            let d2Right22 = e1SqLen_div_dbAreaSq
            let exp4 = d2Left22 * rightTerm + dLeft2 * dRight2.transpose()
            let exp5 = mat2dIdentity * d2Right22 * leftTerm
            let exp6 = dRight2 * dLeft2.transpose()
            curHessian.pointee.block(2, 2, 2, 2) <<== w * (exp4 + exp5 + exp6)
            
            let d2Left23: Mat2d = dAreaRatio_div_dArea_mult * edge_oppo2_Ortho * edge_oppo3_Ortho.transpose() +
                areaRatio * dOrtho_div_dU
            let d2Right23: Double = -e0dote1_div_dbAreaSq
            let exp7 = d2Left23 * rightTerm + dLeft2 * dRight3.transpose()
            let exp8 = mat2dIdentity * d2Right23 * leftTerm
            let exp9 = dRight2 * dLeft3.transpose()
            let result3 = w * (exp7 + exp8 + exp9)
            curHessian.pointee.block(2, 4, 2, 2) <<== result3
            curHessian.pointee.block(4, 2, 2, 2) <<== result3.transpose()
            
            // compute second order derivatives for g_U3
            let d2Left33: Mat2d = dAreaRatio_div_dArea_mult * edge_oppo3_Ortho * edge_oppo3_Ortho.transpose()
            let d2Right33: Double = e0SqLen_div_dbAreaSq
            let exp10 = d2Left33 * rightTerm + dLeft3 * dRight3.transpose()
            let exp11 = mat2dIdentity * d2Right33 * leftTerm
            let exp12 = dRight3 * dLeft3.transpose()
            curHessian.pointee.block(4, 4, 2, 2) <<== w * (exp10 + exp11 + exp12)
            
            // project to nearest SPD matrix
            do {
                try IglUtils.makePD(&curHessian.pointee)
            } catch {
                print(error.localizedDescription)
                fatalError("fix EigenSolver")
            }
            
            
            vInds[triI] = triVInd
            for vI in 0..<3 {
                if (data.fixedVert.contains(vInds[triI][vI])) {
                    vInds[triI][vI] = -1
                }
            }
        }
        
        /*
        for triI in 0..<data.F.rows {
            let triVInd: Vec3i = data.F.row(triI)
            
            let U1: Vec2d = data.V.row(triVInd[0])
            let U2: Vec2d = data.V.row(triVInd[1])
            let U3: Vec2d = data.V.row(triVInd[2])
            
            let U2m1: Vec2d = U2 - U1
            let U3m1: Vec2d = U3 - U1
            
            let area_U = 0.5 * (U2m1[0] * U3m1[1] - U2m1[1] * U3m1[0])
            let areaRatio = data.triAreaSq[triI] / area_U / area_U / area_U
            let dAreaRatio_div_dArea_mult: Double = 3.0 / 2.0 * areaRatio / area_U
            
            let w: Double = (uniformWeight ? 1.0 : (data.triArea[triI] / normalizer_div))
            
            let e0SqLen_div_dbAreaSq: Double = data.e0SqLen_div_dbAreaSq[triI]
            let e1SqLen_div_dbAreaSq: Double = data.e1SqLen_div_dbAreaSq[triI]
            let e0dote1_div_dbAreaSq: Double = data.e0dote1_div_dbAreaSq[triI]
            
            // compute energy terms
            let leftTerm: Double = 1.0 + data.triAreaSq[triI] / area_U / area_U
            let rightTerm: Double = (U3m1.squaredNorm() * e0SqLen_div_dbAreaSq + U2m1.squaredNorm() * e1SqLen_div_dbAreaSq) / 2.0 - U3m1.dot(U2m1) * e0dote1_div_dbAreaSq
            
            let edge_oppo1: Vec2d = U3 - U2
            let edge_oppo2: Vec2d = U1 - U3
            let edge_oppo3: Vec2d = U2 - U1
            let edge_oppo1_Ortho = Vec2d(edge_oppo1[1], -edge_oppo1[0])
            let edge_oppo2_Ortho = Vec2d(edge_oppo2[1], -edge_oppo2[0])
            let edge_oppo3_Ortho = Vec2d(edge_oppo3[1], -edge_oppo3[0])
            let dOrtho_div_dU: Mat2d = .init([0.0, -1.0, 1.0, 0.0], [2, 2])
            
            // compute 1st order derivatives
            let dLeft1: Vec2d = areaRatio * edge_oppo1_Ortho
            let dRight1: Vec2d = ((e0dote1_div_dbAreaSq - e0SqLen_div_dbAreaSq) * U3m1 +
                                  (e0dote1_div_dbAreaSq - e1SqLen_div_dbAreaSq) * U2m1)
            
            let dLeft2: Vec2d = areaRatio * edge_oppo2_Ortho
            let dRight2: Vec2d = (e1SqLen_div_dbAreaSq * U2m1 - e0dote1_div_dbAreaSq * U3m1)
            
            let dLeft3: Vec2d = areaRatio * edge_oppo3_Ortho
            let dRight3 = (e0SqLen_div_dbAreaSq * U3m1 - e0dote1_div_dbAreaSq * U2m1)
            
            // compute second order derivatives for g_U1
            let d2Left11: Mat2d = dAreaRatio_div_dArea_mult * edge_oppo1_Ortho * edge_oppo1_Ortho.transpose()
            let d2Right11: Double = (e0SqLen_div_dbAreaSq + e1SqLen_div_dbAreaSq - 2.0 * e0dote1_div_dbAreaSq)
            let dLeft1dRight1T: Mat2d = dLeft1 * dRight1.transpose()
            var exp1: Mat2d = d2Left11 * rightTerm + dLeft1dRight1T
            var exp2: Mat2d = Mat2d.Identity() * d2Right11 * leftTerm
            var exp3: Mat2d = dLeft1dRight1T.transpose()
            triHessians[triI].block(0, 0, 2, 2) <<== w * (exp1 + exp2 + exp3)
            let d2Left12: Mat2d = dAreaRatio_div_dArea_mult * edge_oppo1_Ortho * edge_oppo2_Ortho.transpose() + areaRatio * dOrtho_div_dU
            let d2Right12: Double = (e0dote1_div_dbAreaSq - e1SqLen_div_dbAreaSq)
            exp1 = d2Left12 * rightTerm + dLeft1 * dRight2.transpose()
            exp2 = Mat2d.Identity() * d2Right12 * leftTerm
            exp3 = dRight1 * dLeft2.transpose()
            triHessians[triI].block(0, 2, 2, 2) <<== w * (exp1 + exp2 + exp3)
            triHessians[triI].block(2, 0, 2, 2) <<== triHessians[triI].block(0, 2, 2, 2).transpose()
            
            let d2Left13: Mat2d = dAreaRatio_div_dArea_mult * edge_oppo1_Ortho * edge_oppo3_Ortho.transpose() +
                                  areaRatio * (-dOrtho_div_dU)
            let d2Right13: Double = (e0dote1_div_dbAreaSq - e0SqLen_div_dbAreaSq)
            exp1 = d2Left13 * rightTerm + dLeft1 * dRight3.transpose()
            exp2 = Mat2d.Identity() * d2Right13 * leftTerm
            exp3 = dRight1 * dLeft3.transpose()
            triHessians[triI].block(0, 4, 2, 2) <<== w * (exp1 + exp2 + exp3)
            triHessians[triI].block(4, 0, 2, 2) <<== triHessians[triI].block(0, 4, 2, 2).transpose()
            
            // compute second order derivatives for g_U2
            let d2Left22: Mat2d = dAreaRatio_div_dArea_mult * edge_oppo2_Ortho * edge_oppo2_Ortho.transpose()
            let d2Right22 = e1SqLen_div_dbAreaSq
            exp1 = d2Left22 * rightTerm + dLeft2 * dRight2.transpose()
            exp2 = Mat2d.Identity() * d2Right22 * leftTerm
            exp3 = dRight2 * dLeft2.transpose()
            triHessians[triI].block(2, 2, 2, 2) <<== w * (exp1 + exp2 + exp3)
            
            let d2Left23: Mat2d = dAreaRatio_div_dArea_mult * edge_oppo2_Ortho * edge_oppo3_Ortho.transpose() +
                areaRatio * dOrtho_div_dU
            let d2Right23: Double = -e0dote1_div_dbAreaSq
            exp1 = d2Left23 * rightTerm + dLeft2 * dRight3.transpose()
            exp2 = Mat2d.Identity() * d2Right23 * leftTerm
            exp3 = dRight2 * dLeft3.transpose()
            triHessians[triI].block(2, 4, 2, 2) <<== w * (exp1 + exp2 + exp3)
            triHessians[triI].block(4, 2, 2, 2) <<== triHessians[triI].block(2, 4, 2, 2).transpose()
            
            // compute second order derivatives for g_U3
            let d2Left33: Mat2d = dAreaRatio_div_dArea_mult * edge_oppo3_Ortho * edge_oppo3_Ortho.transpose()
            let d2Right33: Double = e0SqLen_div_dbAreaSq
            exp1 = d2Left33 * rightTerm + dLeft3 * dRight3.transpose()
            exp2 = Mat2d.Identity() * d2Right33 * leftTerm
            exp3 = dRight3 * dLeft3.transpose()
            triHessians[triI].block(4, 4, 2, 2) <<== w * (exp1 + exp2 + exp3)
            
            // project to nearest SPD matrix
            do {
                try IglUtils.makePD(&triHessians[triI])
            } catch {
                print(error.localizedDescription)
                fatalError("fix EigenSolver")
            }
            
            
            vInds[triI] = triVInd
            for vI in 0..<3 {
                if (data.fixedVert.contains(vInds[triI][vI])) {
                    vInds[triI][vI] = -1
                }
            }
        }*/
        
        for triI in 0..<data.F.rows {
            IglUtils.addBlockToMatrix(triHessians[triI], vInds[triI], 2, &V, &I, &J)
        }
        
        var fixedVertInd = Vec<Int>(data.fixedVert.count)
        var fVI: Int = 0
        for fixedVI in data.fixedVert {
            fixedVertInd[fVI] = fixedVI
            fVI += 1
        }
        IglUtils.addDiagonalToMatrix(Vec<Double>.Ones(data.fixedVert.count * 2), fixedVertInd, 2, &V, &I, &J)
    }
    
    override  func computeHessian(_ data: TriMesh,
                                        _ Hessian: inout Mat<Double>,
                                        _ uniformWeight: Bool = false) {
        
        let normalizer_div: Double = data.surfaceArea
        
        Hessian.resize(data.V.rows * 2, data.V.rows * 2)
        //Hessian.setZero()
        
        //var triHessians: [Matd] = .init(repeating: Matd(6, 6), count: data.F.rows)
        let triHessians: UnsafeMutablePointer<Matd> = .allocate(capacity: data.F.rows)
        //var vInds: [Vec3i] = .init(repeating: Vec3i(), count: data.F.rows)
        let vInds: UnsafeMutablePointer<Vec3i> = .allocate(capacity: data.F.rows)
        defer {
            triHessians.deallocate()
            vInds.deallocate()
        }
        
        DispatchQueue.concurrentPerform(iterations: data.F.rows) { triI in
            (triHessians + triI).initialize(to: Matd(6, 6))
            (vInds + triI).initialize(to: Vec3i())
            let curHessian = triHessians + triI
            
            let triVInd: Vec3i = data.F.row(triI)
            
            let U1: Vec2d = data.V.row(triVInd[0])
            let U2: Vec2d = data.V.row(triVInd[1])
            let U3: Vec2d = data.V.row(triVInd[2])
            
            let U2m1: Vec2d = U2 - U1
            let U3m1: Vec2d = U3 - U1
            
            let area_U: Double = 0.5 * (U2m1[0] * U3m1[1] - U2m1[1] * U3m1[0])
            let areaRatio: Double = data.triAreaSq[triI] / area_U / area_U / area_U
            let dAreaRatio_div_dArea_mult = 3.0 / 2.0 * areaRatio / area_U
            
            let w: Double = (uniformWeight ? 1.0 : (data.triArea[triI] / normalizer_div))
            
            let e0SqLen_div_dbAreaSq: Double = data.e0SqLen_div_dbAreaSq[triI]
            let e1SqLen_div_dbAreaSq: Double = data.e1SqLen_div_dbAreaSq[triI]
            let e0dote1_div_dbAreaSq: Double = data.e0dote1_div_dbAreaSq[triI]
            
            // compute energy terms
            let leftTerm: Double = 1.0 + data.triAreaSq[triI] / area_U / area_U
            let rightTerm: Double = (U3m1.squaredNorm() * e0SqLen_div_dbAreaSq + U2m1.squaredNorm() * e1SqLen_div_dbAreaSq) / 2.0 - U3m1.dot(U2m1) * e0dote1_div_dbAreaSq
            
            let edge_oppo1: Vec2d = U3 - U2
            let edge_oppo2: Vec2d = U1 - U3
            let edge_oppo3: Vec2d = U2 - U1
            let edge_oppo1_Ortho = Vec2d(edge_oppo1[1], -edge_oppo1[0])
            let edge_oppo2_Ortho = Vec2d(edge_oppo2[1], -edge_oppo2[0])
            let edge_oppo3_Ortho = Vec2d(edge_oppo3[1], -edge_oppo3[0])
            let dOrtho_div_dU = Mat2d([0.0, -1.0, 1.0, 0.0], [2, 2])
            
            // compute 1st order derivatives
            let dLeft1: Vec2d = areaRatio * edge_oppo1_Ortho
            let dRight1: Vec2d = ((e0dote1_div_dbAreaSq - e0SqLen_div_dbAreaSq) * U3m1 +
                                  (e0dote1_div_dbAreaSq - e1SqLen_div_dbAreaSq) * U2m1)
            let dLeft2: Vec2d = areaRatio * edge_oppo2_Ortho
            let dRight2: Vec2d = (e1SqLen_div_dbAreaSq * U2m1 - e0dote1_div_dbAreaSq * U3m1)
            
            let dLeft3: Vec2d = areaRatio * edge_oppo3_Ortho
            let dRight3: Vec2d = (e0SqLen_div_dbAreaSq * U3m1 - e0dote1_div_dbAreaSq * U2m1)
            
            // compute second order derivatives for g_U1
            let d2Left11: Mat2d = dAreaRatio_div_dArea_mult * edge_oppo1_Ortho * edge_oppo1_Ortho.transpose()
            let d2Right11: Double = (e0SqLen_div_dbAreaSq + e1SqLen_div_dbAreaSq - 2.0 * e0dote1_div_dbAreaSq)
            let dLeft1dRight1T: Mat2d = dLeft1 * dRight1.transpose()
            var exp1: Mat2d = d2Left11 * rightTerm + dLeft1dRight1T
            var exp2: Mat2d = Mat2d.Identity() * d2Right11 * leftTerm
            curHessian.pointee.block(0, 0, 2, 2) <<== w * (exp1 + exp2 + dLeft1dRight1T.transpose())
            
            let d2Left12: Mat2d = dAreaRatio_div_dArea_mult * edge_oppo1_Ortho * edge_oppo2_Ortho.transpose() +
                areaRatio * dOrtho_div_dU
            let d2Right12: Double = (e0dote1_div_dbAreaSq - e1SqLen_div_dbAreaSq)
            exp1 = (d2Left12 * rightTerm) + (dLeft1 * dRight2.transpose())
            exp2 = Mat2d.Identity() * d2Right12 * leftTerm
            let result1 = w * (exp1 + exp2 + dRight1 * dLeft2.transpose())
            curHessian.pointee.block(0, 2, 2, 2) <<== result1
            curHessian.pointee.block(2, 0, 2, 2) <<== result1.transpose()
            
            let d2Left13: Mat2d = dAreaRatio_div_dArea_mult * edge_oppo1_Ortho * edge_oppo3_Ortho.transpose() +
                areaRatio * (-dOrtho_div_dU)
            let d2Right13: Double = (e0dote1_div_dbAreaSq - e0SqLen_div_dbAreaSq)
            exp1 = d2Left13 * rightTerm + dLeft1 * dRight3.transpose()
            exp2 = Mat2d.Identity() * d2Right13 * leftTerm
            let result2 = w * (exp1 + exp2 + dRight1 * dLeft3.transpose())
            curHessian.pointee.block(0, 4, 2, 2) <<== result2
            curHessian.pointee.block(4, 0, 2, 2) <<== result2.transpose()
            
            // compute second order derivatives for g_U2
            let d2Left22: Mat2d = dAreaRatio_div_dArea_mult * edge_oppo2_Ortho * edge_oppo2_Ortho.transpose()
            let d2Right22: Double = e1SqLen_div_dbAreaSq
            exp1 = d2Left22 * rightTerm + dLeft2 * dRight2.transpose()
            exp2 = Mat2d.Identity() * d2Right22 * leftTerm
            var exp3: Mat2d = dRight2 * dLeft2.transpose()
            curHessian.pointee.block(2, 2, 2, 2) <<== w * (exp1 + exp2 + exp3)
            
            let d2Left23: Mat2d = dAreaRatio_div_dArea_mult * edge_oppo2_Ortho * edge_oppo3_Ortho.transpose() + areaRatio * dOrtho_div_dU
            let d2Right23: Double = -e0dote1_div_dbAreaSq
            exp1 = d2Left23 * rightTerm + dLeft2 * dRight3.transpose()
            exp2 = Mat2d.Identity() * d2Right23 * leftTerm
            exp3 = dRight2 * dLeft3.transpose()
            let result3 = w * (exp1 + exp2 + exp3)
            curHessian.pointee.block(2, 4, 2, 2) <<== result3
            curHessian.pointee.block(4, 2, 2, 2) <<== result3.transpose()
            
            // compute second order derivatives for g_U3
            let d2Left33: Mat2d = dAreaRatio_div_dArea_mult * edge_oppo3_Ortho * edge_oppo3_Ortho.transpose()
            let d2Right33: Double = e0SqLen_div_dbAreaSq
            exp1 = d2Left33 * rightTerm + dLeft3 * dRight3.transpose()
            exp2 = Mat2d.Identity() * d2Right33 * leftTerm
            exp3 = dRight3 * dLeft3.transpose()
            curHessian.pointee.block(4, 4, 2, 2) <<== w * ( exp1 + exp2 + exp3)
            
            // project to nearest SPD matrix
            do {
                try IglUtils.makePD(&curHessian.pointee)
            } catch {
                print(error.localizedDescription)
                fatalError("fix EigenSolver")
            }
            
            vInds[triI] = triVInd
            for vI in 0..<3 {
                if (data.fixedVert.contains(vInds[triI][vI])) {
                    vInds[triI][vI] = -1
                }
            }
        }
        
        /*
        for triI in 0..<data.F.rows { // TODO: Make this for loop parallel
            
            
            
        }*/
        
        for triI in 0..<data.F.rows {
            IglUtils.addBlockToMatrix(triHessians[triI], vInds[triI], 2, &Hessian)
        }
        
        var fixedVertInd = Veci(data.fixedVert.count)
        var fVI: Int = 0
        for fixedVI in data.fixedVert {
            fixedVertInd[fVI] = fixedVI
            fVI += 1
        }
        IglUtils.addDiagonalToMatrix(Vec<Double>.Ones(data.fixedVert.count * 2),
                                     fixedVertInd, 2, &Hessian)
    }
    
    // To prevent element inversion
    override  func initStepSize(_ data: TriMesh,
                                      _ searchDir: Vec<Double>,
                                      _ stepSize: inout Double) {
        
        assert(stepSize > 0.0)
        
        for triI in 0..<data.F.rows {
            let triVInd: Vec3i = data.F.row(triI)
            
            let U1: Vec2d = data.V.row(triVInd[0])
            let U2: Vec2d = data.V.row(triVInd[1])
            let U3: Vec2d = data.V.row(triVInd[2])
            
            let V1 = Vec2d(searchDir[triVInd[0] * 2], searchDir[triVInd[0] * 2 + 1])
            let V2 = Vec2d(searchDir[triVInd[1] * 2], searchDir[triVInd[1] * 2 + 1])
            let V3 = Vec2d(searchDir[triVInd[2] * 2], searchDir[triVInd[2] * 2 + 1])
            
            let U2m1: Vec2d = U2 - U1
            let U3m1: Vec2d = U3 - U1
            let V2m1: Vec2d = V2 - V1
            let V3m1: Vec2d = V3 - V1
            
            let a: Double = V2m1[0] * V3m1[1] - V2m1[1] * V3m1[0]
            let b: Double = U2m1[0] * V3m1[1] - U2m1[1] * V3m1[0] + V2m1[0] * U3m1[1] - V2m1[1] * U3m1[0]
            let c: Double = U2m1[0] * U3m1[1] - U2m1[1] * U3m1[0]
            assert(c > 0.0)
            let delta: Double = b * b - 4.0 * a * c
            
            var bound: Double = stepSize
            if (a > 0.0) {
                if (b < 0.0) && (delta >= 0.0) {
                    bound = 2.0 * c / (-b + sqrt(delta))
                    // (same in math as (-b - sqrt(delta)) / 2.0 / a
                    //  but smaller numerical error when b < 0.0)
                    assert(bound > 0.0)
                }
            } else if (a < 0.0) {
                assert(delta > 0.0)
                if (b < 0.0) {
                    bound = 2.0 * c / (-b + sqrt(delta))
                    // (same in math as (-b - sqrt(delta)) / 2.0 / a
                    //  but smaller numerical error when b < 0.0)
                } else {
                    bound = (-b - sqrt(delta)) / 2.0 / a
                }
                assert(bound > 0.0)
            } else {
                if (b < 0.0) {
                    bound = -c / b
                    assert(bound > 0.0)
                }
            }
            
            if (bound < stepSize) {
                stepSize = bound
            }
        }
    }
    
    override  func checkEnergyVal(_ data: TriMesh) {
        
        let normalizer_div: Double = data.surfaceArea
        
        var energyValPerTri: Vec<Double> = .init(data.F.rows)
        var err: Double = 0.0
        for triI in 0..<data.F.rows {
            let triVInd: Vec3i = data.F.row(triI)
            
            let P1: Vec3d = data.V_rest.row(triVInd[0])
            let P2: Vec3d = data.V_rest.row(triVInd[1])
            let P3: Vec3d = data.V_rest.row(triVInd[2])
            
            // fake isometric UV coordinates
            let P: [Vec3d] = [P1, P2, P3]
            var U: [Vec2d] = .init(repeating: Vec2d(), count: 3)
            IglUtils.mapTriangleTo2D(P, &U)
            let U2m1: Vec2d = U[1]
            let U3m1: Vec2d = U[2]
            
            let area_U: Double = 0.5 * (U2m1[0] * U3m1[1] - U2m1[1] * U3m1[0])
            
            let w: Double = data.triArea[triI] / normalizer_div
            energyValPerTri[triI] = w * (1.0 + data.triAreaSq[triI] / area_U / area_U) *
            ((U3m1.squaredNorm() * data.e0SqLen[triI] * U2m1.squaredNorm() * data.e1SqLen[triI]) / 4.0 / data.triAreaSq[triI] -
             U3m1.dot(U2m1) * data.e0dote1[triI] / 2.0 / data.triAreaSq[triI])
            err += energyValPerTri[triI] - w * 4.0
        }
    }
    
     func getEnergyValPerVert(_ data: TriMesh,
                             _ energyValPerVert: inout Vec<Double>) {
        var energyValPerElem = Vec<Double>()
        getEnergyValPerElem(data, &energyValPerElem)
        
        var totalWeight = Vec<Double>()
        totalWeight.resize(data.V_rest.rows)
        totalWeight.setZero()
        energyValPerVert.resize(data.V_rest.rows)
        energyValPerVert.setZero()
        for triI in 0..<data.F.rows {
            for i in 0..<3 {
                energyValPerVert[data.F[triI, i]] += energyValPerElem[triI]
                totalWeight[data.F[triI, i]] += data.triArea[triI]
                // TODO: verify the scale if the value will be used rather than rank!
            }
        }
        for vI in 0..<data.V_rest.rows {
            energyValPerVert[vI] /= totalWeight[vI] //!!! is normalization needed?
        }
    }
    
     func getMaxUnweightedEnergyValPerVert(_ data: TriMesh,
                                          _ MaxUnweightedEnergyValPerVert: inout Vec<Double>) {
        var energyValPerElem = Vec<Double>()
        getEnergyValPerElem(data, &energyValPerElem, true)
        
        MaxUnweightedEnergyValPerVert.resize(data.V_rest.rows)
        MaxUnweightedEnergyValPerVert.setZero()
        for triI in 0..<data.F.rows {
            for i in 0..<3 {
                if (MaxUnweightedEnergyValPerVert[data.F[triI, i]] < energyValPerElem[triI]) {
                    MaxUnweightedEnergyValPerVert[data.F[triI, i]] = energyValPerElem[triI]
                }
            }
        }
    }
    
     func computeLocalGradient(_ data: TriMesh,
                                     _ localGradients: inout Mat<Double>) {
        
        let normalizer_div: Double = data.surfaceArea
        
        localGradients.resize(data.F.rows * 3, 2)
        for triI in 0..<data.F.rows {
            let triVInd: Vec3i = data.F.row(triI)
            
            let U1: Vec2d = data.V.row(triVInd[0])
            let U2: Vec2d = data.V.row(triVInd[1])
            let U3: Vec2d = data.V.row(triVInd[2])
            
            let U2m1: Vec2d = U2 - U1
            let U3m1: Vec2d = U3 - U1
            
            let area_U: Double = 0.5 * (U2m1[0] * U3m1[1] - U2m1[1] * U3m1[0])
            
            let leftTerm: Double = 1.0 + data.triAreaSq[triI] / area_U / area_U
            let rightTerm: Double = (U3m1.squaredNorm() * data.e0SqLen[triI] + U2m1.squaredNorm() * data.e1SqLen[triI]) / 4 / data.triAreaSq[triI] - U3m1.dot(U2m1) * data.e0dote1[triI] / 2.0 / data.triAreaSq[triI]
            let areaRatio: Double = data.triAreaSq[triI] / area_U / area_U / area_U
            let w: Double = data.triArea[triI] / normalizer_div
            let startRowI: Int = triI * 3
            
            let edge_oppo1: Vec2d = U3 - U2
            let dLeft1: Vec2d = areaRatio * Vec2d(edge_oppo1[1], -edge_oppo1[0])
            let dRight1: Vec2d = ((data.e0dote1[triI] - data.e0SqLen[triI]) * U3m1 + (data.e0dote1[triI] - data.e1SqLen[triI]) * U2m1) / 2.0 / data.triAreaSq[triI]
            localGradients.row(startRowI) <<== w * (dLeft1 * rightTerm + dRight1 * leftTerm)
            
            let edge_oppo2: Vec2d = U1 - U3
            let dLeft2: Vec2d = areaRatio * Vec2d(edge_oppo2[1], -edge_oppo2[0])
            let dRight2: Vec2d = (data.e1SqLen[triI] * U2m1 - data.e0dote1[triI] * U3m1) / 2.0 / data.triAreaSq[triI]
            localGradients.row(startRowI + 1) <<== w * (dLeft2 * rightTerm + dRight2 * leftTerm)
            
            let edge_oppo3: Vec2d = U2 - U1
            let dLeft3: Vec2d = areaRatio * Vec2d(edge_oppo3[1], -edge_oppo3[0])
            let dRight3: Vec2d = (data.e0SqLen[triI] * U3m1 - data.e0dote1[triI] * U2m1) / 2.0 / data.triAreaSq[triI]
            localGradients.row(startRowI + 2) <<== w * (dLeft3 * rightTerm + dRight3 * leftTerm)
        }
    }
    
     func getDivGradPerElem(_ data: TriMesh,
                           _  divGradPerElem: inout Vec<Double>) {
        
        var divGrad_vert = Vec<Double>()
        computeDivGradPerVert(data, &divGrad_vert)
        
        divGradPerElem.resize(data.F.rows)
        for triI in 0..<data.F.rows {
            let triVInd: RVec3i = data.F.row(triI)
            divGradPerElem[triI] = (divGrad_vert[triVInd[0]] + divGrad_vert[triVInd[1]] + divGrad_vert[triVInd[2]]) / 3.0
        }
    }
    
     func computeDivGradPerVert(_ data: TriMesh,
                                      _ divGradPerVert: inout Vec<Double>) {
        
        var localGradients = Matd()
        computeLocalGradient(data, &localGradients)
        
        // NOTE: no need to weight by area since gradient already contains area information
        // TODO: don't need to compute mean because it's always zero at staionary, go back to the simpler version?
        let mean: Matd = .Zero(data.V_rest.rows, 2)
        var incTriAmt: Veci = .Zero(data.V_rest.rows)
        for triI in 0..<data.F.rows {
            let triVInd: RVec3i = data.F.row(triI)
            let locGradStartInd: Int = triI * 3
            for i in 0..<3 {
                mean.row(triVInd[i]) += localGradients.row(locGradStartInd + i)
                incTriAmt[triVInd[i]] += 1
            }
        }
        for vI in 0..<data.V_rest.rows {
            mean.row(vI) /= Double(incTriAmt[vI])
        }
        
        var standardDeviation: Vec<Double> = .Zero(data.V_rest.rows)
        for triI in 0..<data.F.rows {
            let triVInd: RVec3i = data.F.row(triI)
            let locGradStartInd = triI * 3
            for i in 0..<3 {
                standardDeviation[triVInd[i]] += (localGradients.row(locGradStartInd + i) - mean.row(triVInd[i])).squaredNorm()
            }
        }
        
        divGradPerVert = .Zero(data.V_rest.rows)
        for vI in 0..<data.V_rest.rows {
            if (incTriAmt[vI] == 1) {
                // impossible to be splitted
                divGradPerVert[vI] = 0.0
            } else {
                divGradPerVert[vI] = sqrt(standardDeviation[vI] / (Double(incTriAmt[vI]) - 1.0))
            }
        }
    }
}
