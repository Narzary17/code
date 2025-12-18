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

# Configure the number of maximum event the analyser run on in interactive mode
# -1 == ALL
process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(3000) 
)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                
                                #'/store/mc/HINPbPbSpring23MiniAOD/MinBias_Drum5F_5p36TeV_hydjet/MINIAODSIM/NoPU_132X_mcRun3_2023_realistic_HI_v9-v2/2820000/02288831-b588-4380-bac8-9910f24691bd.root',
                                #'/store/hidata/HIRun2024A/HIPhysicsRawPrime0/MINIAOD/PromptReco-v1/000/387/748/00000/c574f914-1d52-4ae8-9bcc-12affad8404b.root'
                                '/store/hidata/HIRun2023A/HIPhysicsRawPrime0/MINIAOD/PromptReco-v2/000/374/668/00000/06179488-b7e6-44f6-bec9-eb242a290ffd.root'
                                #'file:/afs/cern.ch/work/s/subehera/SayanDptDpt/nechstudy/CMSSW_13_2_5_patch1/src/Analyzers/ChFluctuations/cfg/f9e38f93-d48b-4efc-aaa3-033556b66514.root'
                            ),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            skipBadFiles=cms.untracked.bool(True)
)


process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

### centrality ###
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")

# Set the global tag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_HI_v2', '')
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
             tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_Nominal"),
             connect = cms.string("sqlite_file:CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_Nominal.db"),
             label = cms.untracked.string("HFtowers")
    ),
])

