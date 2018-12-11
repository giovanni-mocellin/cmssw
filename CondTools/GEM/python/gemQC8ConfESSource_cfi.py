import FWCore.ParameterSet.Config as cms

from CondCore.CondDB.CondDB_cfi import *
SourceDBConnection = CondDB.clone( connect = cms.string('oracle://INT2R/CMS_GEM_MUON_VIEW') )
SourceDBConnection.DBParameters.authenticationPath = cms.untracked.string('.')
SourceDBConnection.DBParameters.authenticationSystem = cms.untracked.int32(2)

GEMQC8ConfESSource = cms.ESSource("GEMQC8ConfESSource",
    SourceDBConnection,
    runNumber = cms.int32(1),
    DebugMode = cms.untracked.int32(1),
    WriteDummy = cms.untracked.int32(0),#fakeData for testing
    NoDBOutput = cms.untracked.int32(1), # whether PoolDBOutputService is missing
    OnlyConfDef = cms.untracked.int32(1), # whether EMAP should be skipped
    printValues = cms.untracked.bool( True ), # whether to print obtained values

   # geometry definitions
   geomXMLFiles = cms.vstring('Geometry/CMSCommonData/data/materials.xml',
        'Geometry/CMSCommonData/data/rotations.xml',
        'Geometry/MuonCommonData/data/GEMQC8/cms.xml',
        #'Geometry/MuonCommonData/data/GEMQC8/muonBase.xml', # Phase-2 Muon
        #'Geometry/MuonCommonData/data/GEMQC8/cmsMuon.xml',
        #'Geometry/MuonCommonData/data/GEMQC8/mf.xml',
        'Geometry/MuonCommonData/data/gemf/TDR_BaseLine/gemf.xml',
        'Geometry/MuonCommonData/data/GEMQC8/gem11.xml',
        'Geometry/MuonCommonData/data/GEMQC8/muonNumbering.xml',
        'Geometry/MuonCommonData/data/GEMQC8/muonSens.xml',
        'Geometry/MuonCommonData/data/GEMQC8/muonProdCuts.xml',
        'Geometry/MuonCommonData/data/GEMQC8/GEMSpecsFilter.xml',   # Phase-2 Muon
        'Geometry/MuonCommonData/data/GEMQC8/GEMSpecs.xml',
        'Geometry/MuonCommonData/data/GEMQC8/gembox.xml',
        ## these are added on the fly
        ## 'Geometry/MuonCommonData/data/GEMQC8/gem11L_c1_r1.xml',
        ## 'Geometry/MuonCommonData/data/GEMQC8/gem11L_c1_r2.xml',
        ## 'Geometry/MuonCommonData/data/GEMQC8/gem11L_c1_r3.xml',
        ## 'Geometry/MuonCommonData/data/GEMQC8/gem11L_c1_r4.xml',
        ## 'Geometry/MuonCommonData/data/GEMQC8/gem11L_c1_r5.xml',
        ## 'Geometry/MuonCommonData/data/GEMQC8/gem11L_c2_r1.xml',
        ## 'Geometry/MuonCommonData/data/GEMQC8/gem11L_c2_r2.xml',
        ## 'Geometry/MuonCommonData/data/GEMQC8/gem11L_c2_r3.xml',
        ## 'Geometry/MuonCommonData/data/GEMQC8/gem11L_c2_r4.xml',
        ## 'Geometry/MuonCommonData/data/GEMQC8/gem11L_c2_r5.xml',
        ## 'Geometry/MuonCommonData/data/GEMQC8/gem11L_c3_r1.xml',
        ## 'Geometry/MuonCommonData/data/GEMQC8/gem11L_c3_r2.xml',
        ## 'Geometry/MuonCommonData/data/GEMQC8/gem11L_c3_r3.xml',
        ## 'Geometry/MuonCommonData/data/GEMQC8/gem11L_c3_r4.xml',
        ## 'Geometry/MuonCommonData/data/GEMQC8/gem11L_c3_r5.xml',

        ),
    rootNodeName = cms.string('cms:OCMS')
)
