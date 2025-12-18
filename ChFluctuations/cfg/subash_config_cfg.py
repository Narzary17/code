import FWCore.ParameterSet.Config as cms

process = cms.Process("RaghuEffTest")

process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')


# __________________ General _________________

# Configure the logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool( True ),
)

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1) 
)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                
                                '/store/mc/HINOOSpring25MiniAOD/MinBias_OO_5p36TeV_hijing/MINIAODSIM/NoPU_150X_mcRun3_2025_forOO_realistic_v7-v2/90000/7f2d71d7-987c-46d1-95b3-5c0101060cd9.root'
                                ),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            skipBadFiles=cms.untracked.bool(True)
)


# Global Tag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '150X_mcRun3_2025_forOO_realistic_v7', '')

# Output File
process.TFileService = cms.Service("TFileService", fileName = cms.string('OO_MC_output.root'))

process.load("Analyzers.ChFluctuations.ChFluctuations_cfi")

process.p = cms.Path(
    process.defaultCPDC
)
