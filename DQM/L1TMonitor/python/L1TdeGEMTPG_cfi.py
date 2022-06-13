import FWCore.ParameterSet.Config as cms

l1tdeGEMTPGCommon = cms.PSet(
    monitorDir = cms.string("L1TEMU/L1TdeGEMTPG"),
    verbose = cms.bool(False),
    ## when multiple chambers are enabled, order them by station number!
    chambers = cms.vstring("GE11", "GE21"),
    dataEmul = cms.vstring("data","emul"),
    clusterVars = cms.vstring("size", "pad", "bx"),
    clusterNBin = cms.vuint32(20,192,20),
    clusterMinBin = cms.vdouble(0,0,-10),
    clusterMaxBin = cms.vdouble(20,192,10),
    ## GEM reads out only one BX (BX=0) in DAQ path
    ## For debugging purposes, show all the trigger clusters
    useDataClustersOnlyInBX0 = cms.bool(True),
    B904Setup = cms.bool(False),
)

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
l1tdeGEMTPG = DQMEDAnalyzer(
    "L1TdeGEMTPG",
    l1tdeGEMTPGCommon,
    data = cms.InputTag("emtfStage2Digis"),
    emul = cms.InputTag("simMuonGEMPadDigiClusters"),
)
