import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
HFLambdasDump = cms.EDAnalyzer(
    "HFLb2JpsiL0",
    verbose            = cms.untracked.int32(0), 
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    muonQualityString  = cms.untracked.string("AllGlobalMuons"), 
    muonPt             = cms.untracked.double(4.0),
    psiMuons           = cms.untracked.int32(2),
    psiWindow          = cms.untracked.double(0.2),
    L0Window           = cms.untracked.double(0.3),
    LbWindow           = cms.untracked.double(0.7),
    trackPt            = cms.untracked.double(0.5),
    deltaR             = cms.untracked.double(99.0),
    maxDoca            = cms.untracked.double(0.5),
    maxD0              = cms.untracked.double(99.0),
    maxDz              = cms.untracked.double(99.0),
    pvWeight           = cms.untracked.double(0.6),
    pAngle	       = cms.untracked.double(0.02),
    maxVtxChi2         = cms.untracked.double(3.85), # corresponds to prob > 0.05
    type               = cms.untracked.int32(5122)
)

# ######################################################################
# Sequences
# ######################################################################
HFLambdasSequence     = cms.Sequence(HFLambdasDump)

