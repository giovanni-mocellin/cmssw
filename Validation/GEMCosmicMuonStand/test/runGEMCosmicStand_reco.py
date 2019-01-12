# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: SingleMuPt100_cfi -s GEN,SIM,DIGI,L1,DIGI2RAW,RAW2DIGI,L1Reco,RECO --conditions auto:run2_mc --magField 38T_PostLS1 --datatier GEN-SIM --geometry GEMCosmicStand --eventcontent FEVTDEBUGHLT --era phase2_muon -n 100 --fileout out_reco.root
import datetime
print datetime.datetime.now()
import FWCore.ParameterSet.Config as cms
import configureRun_cfi as runConfig
from Configuration.StandardSequences.Eras import eras

process = cms.Process('RECO',eras.phase2_muon)

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

# The superchambers in the 15 slots
SuperChType = runConfig.StandConfiguration

print(SuperChType)

# Calculation of SuperChSeedingLayers from SuperChType
SuperChSeedingLayers = []

for i in range (0,30):
	SuperChSeedingLayers.append(0)

for j in range (0,3):
	for i in range (5*j,5*(j+1)):
		if (SuperChType[i]!='0'):
			SuperChSeedingLayers[i*2]=1
			SuperChSeedingLayers[i*2+1]=3
			break
	for i in range (5*(j+1)-1,5*j-1,-1):
		if (SuperChType[i]!='0'):
			SuperChSeedingLayers[i*2]=4
			SuperChSeedingLayers[i*2+1]=2
			break
			
print(SuperChSeedingLayers)

# import configurations
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Geometry.GEMGeometry.GeometryGEMCosmicStand_cff')
process.load('Configuration.Geometry.GeometryExtended2017Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('Geometry.GEMGeometry.GeometryGEMCosmicStandDB_cff')

# For debug purposes - use what is already in GEM DB
process.GEMQC8ConfESSource.WriteDummy = cms.untracked.int32(-2) # -1 -- P5 chambers, -2 -- special case
process.GEMQC8ConfESSource.runNumber = cms.int32( options.runNum )
process.GEMQC8ConfESSource.printValues = cms.untracked.bool( False )
process.GEMQC8ConfESSource.OnlyConfDef = cms.untracked.int32( 0 )
process.myPrefer = cms.ESPrefer('GEMQC8ConfESSource')

# DEFINITION OF THE SUPERCHAMBERS INSIDE THE STAND
#for i in xrange(len(SuperChType)):
#    column_row = '_c%d_r%d' % ((i/5)+1, i%5+1)
#    if SuperChType[i]=='L' : size = 'L'
#    if SuperChType[i]=='S' : size = 'S'
#    if SuperChType[i]!='0' :
#    	geomFile = 'Validation/GEMCosmicMuonStand/data/gem11'+size+column_row+'.xml'
#    	print(geomFile)
#    if SuperChType[i]!='0' :
#    	process.XMLIdealGeometryESSource.geomXMLFiles.append(geomFile)
#    	print('-> Appended')

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
strOutput = runConfig.OutputFileName
process.load('RecoMuon.TrackingTools.MuonServiceProxy_cff')

# enable debug message logging for our modules
if (options.dumpRaw):
    process.MessageLogger.infos.placeholder = cms.untracked.bool(False)
    process.MessageLogger.infos.INFO = cms.untracked.PSet(limit = cms.untracked.int32(0))

if (options.debug):
    process.MessageLogger.debugModules = cms.untracked.vstring('*')
    process.MessageLogger.cerr.threshold = cms.untracked.string('DEBUG')

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
    
# Reconstruction of muon track

process.MuonServiceProxy.ServiceParameters.Propagators.append('StraightLinePropagator')

process.GEMCosmicMuonForQC8 = cms.EDProducer("GEMCosmicMuonForQC8",
                                       process.MuonServiceProxy,
                                       gemRecHitLabel = cms.InputTag("gemRecHits"),
                                       maxClusterSize = cms.double(runConfig.maxClusterSize),
                                       minClusterSize = cms.double(runConfig.minClusterSize),
                                       trackChi2 = cms.double(runConfig.trackChi2),
                                       trackResX = cms.double(runConfig.trackResX),
                                       trackResY = cms.double(runConfig.trackResY),
                                       MulSigmaOnWindow = cms.double(runConfig.MulSigmaOnWindow),
                                       SuperChamberType = cms.vstring(SuperChType),
                                       SuperChamberSeedingLayers = cms.vdouble(SuperChSeedingLayers),
                                       MuonSmootherParameters = cms.PSet(
                                           PropagatorAlong = cms.string('SteppingHelixPropagatorAny'),
                                           PropagatorOpposite = cms.string('SteppingHelixPropagatorAny'),
                                           RescalingFactor = cms.double(5.0)
                                           ),
                                       )
process.GEMCosmicMuonForQC8.ServiceParameters.GEMLayers = cms.untracked.bool(True)
process.GEMCosmicMuonForQC8.ServiceParameters.CSCLayers = cms.untracked.bool(False)
process.GEMCosmicMuonForQC8.ServiceParameters.RPCLayers = cms.bool(False)

fScale = 1.0

process.gemcrValidation = cms.EDProducer('gemcrValidation',
    process.MuonServiceProxy,
    verboseSimHit = cms.untracked.int32(1),
    simInputLabel = cms.InputTag('g4SimHits',"MuonGEMHits"),
    genVtx = cms.InputTag("generator","unsmeared", "RECO"),
    recHitsInputLabel = cms.InputTag('gemRecHits'),
    tracksInputLabel = cms.InputTag('GEMCosmicMuonForQC8','','RECO'),
    seedInputLabel = cms.InputTag('GEMCosmicMuonForQC8','','RECO'),
    trajInputLabel = cms.InputTag('GEMCosmicMuonForQC8','','RECO'),
    chNoInputLabel = cms.InputTag('GEMCosmicMuonForQC8','','RECO'),
    seedTypeInputLabel = cms.InputTag('GEMCosmicMuonForQC8','','RECO'),
    genParticleLabel = cms.InputTag('genParticles','','RECO'),
    gemDigiLabel = cms.InputTag("muonGEMDigis","","RECO"),
    nBinGlobalZR = cms.untracked.vdouble(200,200,200,150,180,250),
    RangeGlobalZR = cms.untracked.vdouble(564,572,786,794,786,802,110,260,170,350,100,350),
    maxClusterSize = cms.double(runConfig.maxClusterSize),
    minClusterSize = cms.double(runConfig.minClusterSize),
    maxResidual = cms.double(runConfig.maxResidual),
    isMC = cms.bool(True),
    SuperChamberType = cms.vstring(SuperChType),
    SuperChamberSeedingLayers = cms.vdouble(SuperChSeedingLayers),
    MuonSmootherParameters = cms.PSet(
                      PropagatorAlong = cms.string('SteppingHelixPropagatorAny'),
                      PropagatorOpposite = cms.string('SteppingHelixPropagatorAny'),
                      RescalingFactor = cms.double(5.0)
                      )
)

process.TFileService = cms.Service("TFileService",
    							   fileName = cms.string('temp_'+strOutput)
								   )

# Path and EndPath definitions
process.rawTOhits_step = cms.Path(process.dumpRaw+process.muonGEMDigis+process.gemRecHits)
process.reconstruction_step = cms.Path(process.GEMCosmicMuonForQC8)
process.validation_step = cms.Path(process.gemcrValidation)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(process.rawTOhits_step,
								process.reconstruction_step,
                                process.validation_step,
                                process.endjob_step
                                )
                                
# enable validation event filtering
if (not options.valEvents):
    process.rawTOhits_step.remove(process.validationEventFilter)

# enable validation event filtering
if (len(options.mps)==0):
    process.rawTOhits_step.remove(process.tmtFilter)

# enable RAW printout
if (not options.dumpRaw):
    process.rawTOhits_step.remove(process.dumpRaw)

process.gemSegments.maxRecHitsInCluster = cms.int32(10)
process.gemSegments.minHitsPerSegment = cms.uint32(3)
process.gemSegments.clusterOnlySameBXRecHits = cms.bool(True)
process.gemSegments.dEtaChainBoxMax = cms.double(1.05)
process.gemSegments.dPhiChainBoxMax = cms.double(1.12)
process.gemSegments.dXclusBoxMax = cms.double(10.0)
process.gemSegments.dYclusBoxMax = cms.double(50.0)
process.gemSegments.preClustering = cms.bool(False)
process.gemSegments.preClusteringUseChaining = cms.bool(False)
