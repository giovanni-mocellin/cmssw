set -x
nEvents=$1
if [ ${#nEvents} -eq 0 ] ; then nEvents=3; fi
inpFiles=$2
if [ ${#inpFiles} -eq 0 ] ; then inpFiles="file:run000000_Cosmic_CERNQC8_2018-11-09_chunk_321.root"; fi
runNum=$3
if [ ${#runNum} -eq 0 ] ; then runNum=1; fi

cmsRun unpackData-QC8try.py inputFiles=${inpFiles} localMode=True maxEvents=${nEvents} edm=True runNum=${runNum} #dumpRaw=True

cp gem_EDM-qc8spec.root gem_EDM-qc8spec-runNum${runNum}.root

cmsRun test_qc8_cfg.py runNum=${runNum}

cp DQM_V0001_GEM_R000000001.root DQM_V0001_GEM_R000000001-runNum${runNum}.root
