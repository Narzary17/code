import FWCore.ParameterSet.Config as cms

chFluctuations = cms.EDAnalyzer('ChFluctuations',
    # Input Tags
    tracks    = cms.InputTag("packedPFCandidates"),
    tracksgen = cms.InputTag("packedGenParticles"),
    srcTower  = cms.InputTag("packedPFCandidates"),
    vertex    = cms.InputTag("offlineSlimmedPrimaryVertices"),
    centralitySrc = cms.InputTag("hiCentrality"),

    IsMC = cms.untracked.bool(False), 

    # Selection defaults
    noffmin       = cms.untracked.int32(0),
    noffmax       = cms.untracked.int32(10000),
    ptnoffmin     = cms.untracked.double(0.4),
    ptnoffmax     = cms.untracked.double(10000.0),
    dzdzerrornoff = cms.untracked.double(3.0),
    d0d0errornoff = cms.untracked.double(3.0),
    pterrorptnoff = cms.untracked.double(0.1),

    # Analysis Cuts
    etamin    = cms.untracked.double(-2.4),
    etamax    = cms.untracked.double(2.4),
    ptmin     = cms.untracked.double(0.3),
    ptmax     = cms.untracked.double(3.0),

    # Efficiency
    cweight   = cms.untracked.bool(True),
    fname     = cms.untracked.InputTag('Eff_OO_2025_Hijing_MB_Centrality_NoPU_3D_Nominal_Official.root'),

    # Binning
    etaBins = cms.untracked.vdouble(-3.0, -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0),
    centbins = cms.untracked.vdouble(0, 10, 20, 40, 60, 80, 100, 120, 140, 160)
)
