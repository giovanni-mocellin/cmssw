import FWCore.ParameterSet.Config as cms

process = cms.Process('DQMTEST')

# options
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('analysis')

options.register('runNum',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Run number")
options.parseArguments()




process.MessageLogger = cms.Service("MessageLogger",
  statistics = cms.untracked.vstring(),
  destinations = cms.untracked.vstring('cerr'),
  cerr = cms.untracked.PSet(
#      threshold = cms.untracked.string('WARNING')
      threshold = cms.untracked.string('DEBUG')
  )
)

process.load('Configuration.Geometry.GeometryExtended2017Reco_cff')
process.load("DQM.Integration.config.FrontierCondition_GT_cfi")

process.load('Geometry.GEMGeometry.GeometryGEMCosmicStandDB_cff')

# For debug purposes - use what is already in GEM DB
process.GEMQC8ConfESSource.WriteDummy = cms.untracked.int32(-2) #-1 -- P5 chambers, -2 -- special case
process.GEMQC8ConfESSource.runNumber = cms.int32( options.runNum )
process.GEMQC8ConfESSource.printValues = cms.untracked.bool( False )
process.myPrefer = cms.ESPrefer('GEMQC8ConfESSource')



process.load("DQM.Integration.config.environment_cfi")
process.dqmEnv.subSystemFolder = "GEM"
process.dqmEnv.eventInfoFolder = "EventInfo"
process.dqmSaver.path = ""
process.dqmSaver.tag = "GEM"

process.source = cms.Source("PoolSource",
  labelRawDataLikeMC = cms.untracked.bool(False),
  fileNames = cms.untracked.vstring(
#    'file:/eos/cms/store/express/Commissioning2018/ExpressCosmics/FEVT/Express-v1/000/310/292/00000/6C23251D-4F18-E811-AEC5-02163E01A41D.root'
#        'file:gem_EDM.root'
        'file:../../../EventFilter/GEMRawToDigi/test/gem_EDM-qc8spec.root'
  ),
  #inputCommands = cms.untracked.vstring(
  #  'drop *',
  #  'keep FEDRawDataCollection_*_*_*'
  #)
)

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)

#process.load("EventFilter.GEMRawToDigi.muonGEMDigis_cfi")
#process.load('RecoLocalMuon.GEMRecHit.gemRecHits_cfi')
process.load("DQM.GEM.GEMDQM_cff")

#process.muonGEMDigis.useDBEMap = True
#process.muonGEMDigis.unPackStatusDigis = True

process.path = cms.Path(
  #process.muonGEMDigis *
  #process.gemRecHits *
  process.GEMDQM
)

process.end_path = cms.EndPath(
  process.dqmEnv +
  process.dqmSaver
)

process.schedule = cms.Schedule(
  process.path,
  process.end_path
)
