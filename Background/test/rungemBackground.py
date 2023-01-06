import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras


process = cms.Process("gemBackground", eras.phase2_muon)

process.load('gemsw.Geometry.GeometryTestBeam2022_cff')

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )

process.source = cms.Source(
        "PoolSource",
        fileNames = cms.untracked.vstring('file:/home/jheo/gem-background/data/357696/03c7e084-84a1-4658-852c-3b87f8049886.root')
        )

process.TFileService = cms.Service("TFileService", fileName=cms.string("hits.root"))

process.gembackground = cms.EDAnalyzer('gemBackground',
        gemRecHits = cms.InputTag("gemRecHits"),
        onlineMetaDataDigis = cms.InputTag("onlineMetaDataDigis"),
        )

process.p = cms.Path(process.gembackground)
