//
//  Shared.swift
//  OptCuts
//
//  Created by Reza on 6/23/23.
//

import Foundation
import Matrix



 let __DBL_MAX__: Double = Double.greatestFiniteMagnitude
 var methodType: MethodType = .MT_OPTCUTS
 var filterExp_in: Double = 0.6

 var energyChanges_bSplit: [Pair<Double, Double>] = []
 var energyChanges_iSplit: [Pair<Double, Double>] = []
 var energyChanges_merge: [Pair<Double, Double>] = []

 var paths_bSplit: [[Int]] = []
 var paths_iSplit: [[Int]] = []
 var paths_merge: [[Int]] = []

 var newVertPoses_bSplit: [Mat<Double>] = []
 var newVertPoses_iSplit: [Mat<Double>] = []
 var newVertPoses_merge: [Mat<Double>] = []

 var inSplitTotalAmt: Int = 0

 var triSoup: [TriMesh] = []
 var channel_result: Int = 1

 var energyParams: [UnsafeMutablePointer<Double>] = []

 var V = Matd()
 var UV = Matd()
 var N = Matd()
 var F = Mati()
 var FUV = Mati()
 var FN = Mati()

// optimization
 var vertAmt_input: Int = 0
 var triSoup_backup: TriMesh = .init()
 var optimizer: Optimizer!
 var energyTerms: [Energy] = []

 var bijectiveParam: Bool = true
 var rand1PinitCut: Bool = false
 var lambda_init: Double = .zero
 var optimization_on: Bool = false
 var iterNum: Int = 0
 var converged: Int = 0
 var fractureMode: Bool = false
 var fracThres: Double = 0.0
 var topoLineSearch: Bool = true
 var initCutOption: Int = 0
 var outerLoopFinished: Bool = false
 var upperBound: Double = 4.1
 let convTol_upperBound: Double = 1.0e-3

 var opType_queried: Int = -1
 var path_queried: [Int] = []
 var newVertPos_queried: Matd = .init()
 var reQuery: Bool = false

// visualization
 var headlessMode: Bool = false
 var channel_initial: Int = 0
 var channedl_findExtrema: Int = 2
 var viewChannel: Int = 1
 var viewUV: Bool = true
 var texScale: Double = 1.0
 var showSeam: Bool = true
 var seamColor = Matd()
 var showDistortion: Int = 1 // 0: don't show; 1: SD energy value; 2: L2 stretch value
 var showTexture: Bool = true // show checkerboard
 var isLighting: Bool = false
 var showFracTail: Bool = true
 var fracTailSize: Float = 15.0
 var offlineMode: Bool = false
 var saveInfo_postDraw: Bool = false
 var infoName: String = ""
 var isCapture3D: Bool = false
 var capture3DI: Int = 0

 var secPast: Double = 0.0

// updateLambda_stationaryV static variables
 var iterNum_bestFeasible: Int = -1
 var triSoup_bestFeasible: TriMesh = .init()
 var E_se_bestFeasible: Double = __DBL_MAX__
 var lastStationaryIterNum: Int = 0 // still necessary becauseboundary and interior query are with same iterNum
 var configs_stationaryV: [Double : [Pair<Double, Double>]] = [:]
// save info at first feasible stationary VT for comparison
 var saved: Bool = false
 var start: DispatchTime!

// predraw func statics
// save info once bound is reached for comparison
 var preDrawSaved: Bool = false
 var outputFolderPath: URL!
