import FWCore.ParameterSet.Config as cms

process = cms.Process("ChFluctMC")

# -------------------- Standard Setup --------------------
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(3000))

# -------------------- Input files --------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/HINPbPbSpring23MiniAOD/MinBias_Drum5F_5p36TeV_hydjet/MINIAODSIM/NoPU_132X_mcRun3_2023_realistic_HI_v9-v2/2820000/02288831-b588-4380-bac8-9910f24691bd.root'
    ),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    skipBadFiles = cms.untracked.bool(True)
)

# -------------------- Centrality --------------------
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")

# -------------------- Global Tag --------------------
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_mcRun3_2023_realistic_HI_v5', '')
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(
        record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_Nominal"),
        connect = cms.string("sqlite_file:CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_Nominal.db"),
        label = cms.untracked.string("HFtowers")
    ),
])

# -------------------- Event Selection --------------------
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.load('HeavyIonsAnalysis.EventAnalysis.hffilter_cfi')

process.eventFilter = cms.Sequence(
    process.phfCoincFilter2Th4 *
    process.primaryVertexFilter *
    process.clusterCompatibilityFilter
)

# -------------------- Output --------------------
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('Test_gen_reco_mc.root')
)

# -------------------- Analyzer --------------------
process.load("Analyzers.ChFluctuations.ChFluctuations_cfi")

process.defaultCPDC.tracks = cms.InputTag("packedPFCandidates")
process.defaultCPDC.tracksgen = cms.InputTag("packedGenParticles")
process.defaultCPDC.trackschi2 = cms.InputTag("packedPFCandidateTrackChi2")
process.defaultCPDC.vertex = cms.InputTag("offlineSlimmedPrimaryVertices")
process.defaultCPDC.centralitybin = cms.InputTag("centralityBin", "HFtowers")
process.defaultCPDC.IsMC = cms.untracked.bool(True)

# pt, eta, vertex cuts
process.defaultCPDC.ptMin = cms.untracked.double(0.3)
process.defaultCPDC.ptMax = cms.untracked.double(3.0)
process.defaultCPDC.etaMax = cms.untracked.double(2.4)
process.defaultCPDC.zminVtx = cms.untracked.double(-15.0)
process.defaultCPDC.zmaxVtx = cms.untracked.double(15.0)

# Efficiency files (correct InputTag for MiniAOD)
process.defaultCPDC.fGeneral = cms.untracked.InputTag(
    "/afs/cern.ch/user/a/anarzary/CMSSW_13_2_5_patch1/src/Analyzers/ChFluctuations/data/EFF/general_tracks/General_TrackEff_3D_nhits0_noChi2Cut_ptGt10ptResoCut_Nominal_NewSample.root"
)
process.defaultCPDC.fGeneral2 = cms.untracked.InputTag(
    "/afs/cern.ch/user/a/anarzary/CMSSW_13_2_5_patch1/src/Analyzers/ChFluctuations/data/EFF/general_tracks/General_TrackEff_3D_nhits0_noChi2Cut_ptGt10ptResoCut_Nominal_NewSampleQCD.root"
)

# Binning
process.defaultCPDC.pTBins = cms.untracked.vdouble(
    0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7,
    0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
    1.8, 1.9, 2.0, 2.5, 3.0
)
process.defaultCPDC.etaBins = cms.untracked.vdouble(
    -3.0, -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0
)
process.defaultCPDC.centbins = cms.untracked.vdouble(
    0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0
)
process.defaultCPDC.algoParameters = cms.vint32(*range(47))

# -------------------- Path --------------------
process.p = cms.Path(
    #process.eventFilter *   # optional: enable if you want event selection
    process.centralityBin *
    process.defaultCPDC
)

