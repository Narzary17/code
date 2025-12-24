import FWCore.ParameterSet.Config as cms

process = cms.Process("OOFlowAnalysisMC")

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

# mc dataset
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # Using Global Redirector to ensure access
        'root://cms-xrd-global.cern.ch//store/mc/HINOOSpring25MiniAOD/MinBias-Angantyr_5p36TeV_pythia8/MINIAODSIM/NoPU_CustomTrack_NoPU_CustomTrack_150X_mcRun3_2025_forOO_realistic_v8-v3/2520000/01e4c607-8a4a-4a4c-8c25-7577364df774.root'
    ),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    skipBadFiles=cms.untracked.bool(True)
)

# Global Tag 
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '150X_mcRun3_2025_forOO_realistic_v8', '')

# Trigger
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltMB = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltMB.HLTPaths = ["HLT_MinimumBiasHF_OR_BptxAND_v1"] 
process.hltMB.andOr = cms.bool(True)
process.hltMB.throw = cms.bool(False)
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.load('HeavyIonsAnalysis.EventAnalysis.hffilterPF_cfi')

from HeavyIonsAnalysis.TrackAnalysis.unpackedTracksAndVertices_cfi import unpackedTracksAndVertices
process.unpackedTracksAndVertices = unpackedTracksAndVertices

process.load('HeavyIonsAnalysis.VertexAnalysis.pileupvertexfilter_cfi')
process.pileupvertexfilter.doOO = True

process.eventFilter_HM = cms.Sequence(
    process.phfCoincFilterPF2Th4 *
    process.primaryVertexFilter *
    process.unpackedTracksAndVertices *
    process.pileupvertexfilter *
    process.clusterCompatibilityFilter
)

process.TFileService = cms.Service("TFileService", fileName = cms.string('OO_MC_Output_Angantyr.root'))

process.load('Analyzers.Cumulants.ChFluctuations_cfi')

process.defaultCPDC = process.chFluctuations.clone(
    IsMC = cms.untracked.bool(False),  
    noffmin = cms.untracked.int32(10)  
)

process.p = cms.Path(
    # process.hltMB *           
    process.eventFilter_HM *    
    process.defaultCPDC         
)
