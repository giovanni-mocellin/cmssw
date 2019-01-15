# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: SingleElectronPt10_cfi.py -s GEN,SIM,DIGI,L1 --pileup=NoPileUp --geometry DB --conditions=auto:startup -n 1 --no_exec
import FWCore.ParameterSet.Config as cms


# options
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('analysis')
options.register('skipEvents',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Number of events to skip")
options.register('streamer',
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Read input from streamer file")
options.register('localMode',
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Read input from file produced in 'local-mode'")
options.register('debug',
                 True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Enable debug data")
options.register('dumpRaw',
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Print RAW data")
options.register('dumpDigis',
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Print digis")
options.register('histos',
                 True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Produce standard histograms")
options.register('edm',
                 True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Produce EDM file")
options.register('valEvents',
                 True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Filter on validation events")
options.register('process',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Rename process if used")
options.register('mps',
                 '',
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.int,
                 "List of MPs to process")
options.register('json',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "JSON file with list of good lumi sections")
options.register('evtDisp',
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 'Produce histos for individual events')
options.register('runNum',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Run number")

options.parseArguments()


pname="Raw2Digi"
if (options.process!=""):
    pname=options.process
#process = cms.Process(pname)

#from Configuration.StandardSequences.Eras import eras
#process = cms.Process(pname, eras.Run2_2017, eras.run2_GEM_2017)
from Configuration.StandardSequences.Eras import eras

process = cms.Process('RECO',eras.Run2_2017)
#process = cms.Process('RECO',eras.Run2_2017,eras.run2_GEM_2017)

process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
#process.load('Configuration.StandardSequences.AlCaRecoStreams_cff')


# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.Geometry.GeometryExtended2017Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('Geometry.GEMGeometry.GeometryGEMCosmicStandDB_cff')

# For debug purposes - use what is already in GEM DB
process.GEMQC8ConfESSource.WriteDummy = cms.untracked.int32(-2) # -1 -- P5 chambers, -2 -- special case
process.GEMQC8ConfESSource.runNumber = cms.int32( options.runNum )
process.GEMQC8ConfESSource.printValues = cms.untracked.bool( False )
process.GEMQC8ConfESSource.OnlyConfDef = cms.untracked.int32( 0 )
#process.myPrefer = cms.ESPrefer('DDCompactView','GEMQC8ConfESSource')
process.myPrefer = cms.ESPrefer('GEMQC8ConfESSource')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
if (options.localMode) :
    process.source = cms.Source(
        "GEMLocalModeDataSource",
        fileNames = cms.untracked.vstring (options.inputFiles),
        skipEvents=cms.untracked.uint32(options.skipEvents),
        fedId = cms.untracked.int32( 1472 ),  # which fedID to assign
        runNumber = cms.untracked.int32( options.runNum ),
        #processEvents =cms.untracked.vuint32( 118 , 153 , 250 , 282 , 533 , 534 , 595 , 603 , 630 , 794 , 797 , 885 , 915 , 928 , 1128 , 1151 , 1269 , 1302 , 1377 , 1630 , 1883 , 1946 , 1988 , 2197 , 2384 , 2387 , 2405 , 2424 , 2458 , 2461 , 2813 , 2827 , 2987 , 3067 , 3120 , 3208 , 3209 , 3215 , 3265 , 3305 , 3422 , 3628 , 3649 , 3706 , 3831 , 3979 , 4137 , 4290 , 4293 , 4495 , 4533 , 4585 , 4633 , 4656 , 4713 , 4738 , 4817 , 4850 , 4860 , 4887 , 4895 , 4922 , 5013 , 5030 , 5068 , 5085 , 5120 , 5167 , 5300 , 5366 , 5898 , 5953 , 6295 , 6403 , 6494 , 6574 , 6763 , 6792 , 6844 , 6852 , 6867 , 6879 , 6900 , 7004 , 7024 , 7109 , 7115 , 7123 , 7236 , 7353 , 7408 , 7495 , 7500 , 7673 )
    )
elif (options.streamer) :
    process.source = cms.Source(
        "NewEventStreamFileReader",
        fileNames = cms.untracked.vstring (options.inputFiles),
        skipEvents=cms.untracked.uint32(options.skipEvents)
    )
else :
    process.source = cms.Source (
        "PoolSource",
        fileNames = cms.untracked.vstring (options.inputFiles),
        skipEvents=cms.untracked.uint32(options.skipEvents)
    )

if (options.json):
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = options.json).getVLuminosityBlockRange()

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)


# Additional output definition
# TTree output file
process.load("CommonTools.UtilAlgos.TFileService_cfi")

# enable debug message logging for our modules
#process.MessageLogger.categories.append('L1TCaloEvents')

if (options.dumpRaw):
    process.MessageLogger.infos.placeholder = cms.untracked.bool(False)
    process.MessageLogger.infos.INFO = cms.untracked.PSet(limit = cms.untracked.int32(0))

if (options.debug):
#    process.MessageLogger.debugModules = cms.untracked.vstring('L1TRawToDigi:caloStage2Digis', 'MP7BufferDumpToRaw:stage2MPRaw', 'MP7BufferDumpToRaw:stage2DemuxRaw')
    process.MessageLogger.debugModules = cms.untracked.vstring('*')
    process.MessageLogger.cerr.threshold = cms.untracked.string('DEBUG')


# Other statements
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')

# validation event filter
process.load('EventFilter.L1TRawToDigi.validationEventFilter_cfi')

# MP selectah
process.load('EventFilter.L1TRawToDigi.tmtFilter_cfi')
process.tmtFilter.mpList = cms.untracked.vint32(options.mps)

# dump raw data
process.dumpRaw = cms.EDAnalyzer(
    "DumpFEDRawDataProduct",
    #label = cms.untracked.string("rawDataCollector"),
    inputTag = cms.untracked.InputTag("source","gemLocalModeDataSource"),
    feds = cms.untracked.vint32 ( 1472 ),
    dumpPayload = cms.untracked.bool ( options.dumpRaw )
)

# raw to digi
process.load('EventFilter.GEMRawToDigi.muonGEMDigis_cfi')
#process.load('EventFilter.GEMRawToDigi.GEMSQLiteCabling_cfi')
process.muonGEMDigis.InputLabel = cms.InputTag("source","gemLocalModeDataSource")
process.muonGEMDigis.useDBEMap = True
process.muonGEMDigis.unPackStatusDigis = True

#process.load('Geometry.GEMGeometryBuilder.gemGeometry_cfi')
process.load('RecoLocalMuon.GEMRecHit.gemRecHits_cfi')

process.gemRecHits = cms.EDProducer("GEMRecHitProducer",
    recAlgoConfig = cms.PSet(),
    recAlgo = cms.string('GEMRecHitStandardAlgo'),
    gemDigiLabel = cms.InputTag("muonGEMDigis"),
    # maskSource = cms.string('File'),
    # maskvecfile = cms.FileInPath('RecoLocalMuon/GEMRecHit/data/GEMMaskVec.dat'),
    # deadSource = cms.string('File'),
    # deadvecfile = cms.FileInPath('RecoLocalMuon/GEMRecHit/data/GEMDeadVec.dat')
)

process.reader_qc8conf = cms.EDAnalyzer( "GEMQC8ConfRcdReader",
  dumpFileName = cms.untracked.string( "dumpQC8conf-qc8spec.out" )
  #dumpFileName = cms.untracked.string( "" ) # no dump
)

process.reader_elmap = cms.EDAnalyzer( "GEMELMapRcdReader",
  dumpFileName = cms.untracked.string( "dumpELMap-from-qc8spec.out" )
  #dumpFileName = cms.untracked.string( "" ) # no dump
)

process.load("DQM.GEM.GEMDQM_cff")
process.load("DQM.Integration.config.environment_cfi")
process.dqmEnv.subSystemFolder = "GEM"
process.dqmEnv.eventInfoFolder = "EventInfo"
process.dqmSaver.path = ""
process.dqmSaver.tag = "GEM"

# Path and EndPath definitions
process.path = cms.Path(
    #process.validationEventFilter
    process.dumpRaw
    +process.muonGEMDigis
    #+process.reader_elmap
    #+process.reader_qc8conf
    +process.gemRecHits
    +process.GEMDQM
)

process.end_path = cms.EndPath(
  process.dqmEnv +
  process.dqmSaver
)

process.schedule = cms.Schedule(
  process.path,
  process.end_path
)


# enable validation event filtering
if (not options.valEvents):
    process.path.remove(process.validationEventFilter)

# enable validation event filtering
if (len(options.mps)==0):
    process.path.remove(process.tmtFilter)

# enable RAW printout
if (not options.dumpRaw):
    process.path.remove(process.dumpRaw)

# optional EDM file
if (options.edm):
    process.output = cms.OutputModule(
        "PoolOutputModule",
        outputCommands = cms.untracked.vstring("keep *"),
        SelectEvents = cms.untracked.PSet(
            SelectEvents = cms.vstring('path')
        ),
        fileName = cms.untracked.string('gem_EDM-qc8spec.root')
    )

    process.out = cms.EndPath(
        process.output
    )


# print ESModules
for name, module in process.es_sources_().iteritems():
    if 'GEM' in name:
        print "ESModules> provider:%s '%s'" % ( name, module.type_() )
for name, module in process.es_producers_().iteritems():
    if 'GEM' in name:
        print "ESModules> provider:%s '%s'" % ( name, module.type_() )
