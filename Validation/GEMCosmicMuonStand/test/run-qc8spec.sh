set -x
#cmsRun unpackData-QC8try.py inputFiles=file:/afs/cern.ch/user/d/dorney/public/v3Hack/run000000_Cosmic_CERNQC8_2018-11-09_chunk_321.dat streamer=True maxEvents=1 edm=True
#cmsRun unpackData-QC8try.py inputFiles="file:/tmp/andriusj/run000000_Cosmic_CERNQC8_2018-11-09_chunk_321.root,file:/tmp/andriusj/tmp.root" localMode=True maxEvents=100 edm=True #dumpRaw=True
#cmsRun unpackData-QC8try.py inputFiles="file:run000000_Cosmic_CERNQC8_2018-11-09_chunk_321.root" localMode=True maxEvents=100 edm=True #dumpRaw=True
nEvents=$1
if [ ${#nEvents} -eq 0 ] ; then nEvents=3; fi
inpFiles=$2
if [ ${#inpFiles} -eq 0 ] ; then inpFiles="file:run000000_Cosmic_CERNQC8_2018-11-09_chunk_321.root"; fi
runNum=$3
if [ ${#runNum} -eq 0 ] ; then runNum=1; fi
cmsRun unpackData-QC8try.py inputFiles=${inpFiles} localMode=True maxEvents=${nEvents} edm=True runNum=${runNum} #dumpRaw=True
#cmsRun unpackData-QC8try.py inputFiles="inpFiles.lst" localMode=True maxEvents=${nEvents} edm=True #dumpRaw=True

cp gem_EDM-qc8spec.root gem_EDM-qc8spec-runNum${runNum}.root
