import FWCore.ParameterSet.Config as cms

process = cms.Process("RaghuEffTest")

process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/hidata/OORun2025/IonPhysics1/MINIAOD/PromptReco-v1/000/394/217/00000/fdc2d47a-9fc0-4ef3-ace4-400e62e785cb.root'
    ),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    skipBadFiles=cms.untracked.bool(True)
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Express_v4', '')

# Define Output
process.TFileService = cms.Service("TFileService", fileName = cms.string('OO_data_output.root'))

# Trigger Selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltMB = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltMB.HLTPaths = ["HLT_MinimumBiasHF_OR_BptxAND_v1"]
process.hltMB.andOr = cms.bool(True)
process.hltMB.throw = cms.bool(False)

# Analyzer
process.load("Analyzers.ChFluctuations.ChFluctuations_cfi")

process.p = cms.Path(
    process.hltMB *
    process.defaultCPDC
)
