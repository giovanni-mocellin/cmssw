import time
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from CondCore.CondDB.CondDB_cfi import *

options = VarParsing.VarParsing()
options.register('connectionString',
                 'sqlite_file:GEM-QC8Conf.db', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Connection string")
options.register('tag',
                 'GEMQC8Conf_v1', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Tag")
options.register('runNumber',
                 1, #default value, int
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Run number; default gives latest IOV")
options.register('numberOfRuns',
                 1, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "number of runs in the job")
options.register('messageLevel',
                 0, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Message level; default to 0")
options.parseArguments()

CondDBReference = CondDB.clone( connect = cms.string( options.connectionString ) )
CondDBReference.DBParameters.messageLevel = cms.untracked.int32( options.messageLevel )

process = cms.Process( "DBTest" )

process.MessageLogger = cms.Service( "MessageLogger",
                                     destinations = cms.untracked.vstring( 'cout' ),
                                     cout = cms.untracked.PSet( threshold = cms.untracked.string( 'INFO' ) ),
                                     )

process.source = cms.Source("EmptyIOVSource",
                            timetype = cms.string('runnumber'),
                            firstValue = cms.uint64(options.runNumber),
                            lastValue = cms.uint64(options.runNumber + options.numberOfRuns),
                            interval = cms.uint64(1)
                        )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( options.numberOfRuns ) ) #options.numberOfRuns runs per job

process.GEMQC8ConfESSource = cms.ESSource( "PoolDBESSource",
                                   CondDBReference,
                                   toGet = cms.VPSet( cms.PSet( record = cms.string('GEMQC8ConfRcd'),
                                                                tag = cms.string('GEMQC8Conf_v1')
                                                                )
                                                      ),
                                   )

process.reader = cms.EDAnalyzer( "GEMQC8ConfRcdReader",
                                 dumpFileName = cms.untracked.string( "dumpQC8conf.out" )
)

process.recordDataGetter = cms.EDAnalyzer( "EventSetupRecordDataGetter",
                                           toGet =  cms.VPSet(),
                                           verbose = cms.untracked.bool( True )
                                           )
process.escontent = cms.EDAnalyzer( "PrintEventSetupContent",
                                    compact = cms.untracked.bool( True ),
                                    printProviders = cms.untracked.bool( True )
                                    )
process.esretrieval = cms.EDAnalyzer( "PrintEventSetupDataRetrieval",
                                      printProviders = cms.untracked.bool( True )
                                      )

#Path definition
process.GEMQC8ConfReaderSourcePath = cms.Path( process.reader + process.recordDataGetter )
process.esout = cms.EndPath( process.escontent + process.esretrieval )

#Schedule definition
process.schedule = cms.Schedule( process.GEMQC8ConfReaderSourcePath,
                                 process.esout
                                 )

for name, module in process.es_sources_().iteritems():
    print "ESModules> provider:%s '%s'" % ( name, module.type_() )
for name, module in process.es_producers_().iteritems():
    print "ESModules> provider:%s '%s'" % ( name, module.type_() )
