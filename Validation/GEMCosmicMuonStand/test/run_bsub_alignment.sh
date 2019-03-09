#! /bin/sh

NeventsPerJob=$1
RUN=$2
JobId=$3
#echo geomType${geomType} column${column} row${row} trueDx${tDx} trueDy${tDy} trueDz${tDz} trueRz${tRz} NeventsPerJob${NeventsPerJob} run$RUN jobid$JobId

srcFile=runGEMCosmicStand_sim_align.py

srcDir=/afs/cern.ch/user/j/jslee/CMSSW_10_1_0/src
workDir=$srcDir/Validation/GEMCosmicMuonStand/test
tempDir=/tmp/jslee
resultDir=/afs/cern.ch/work/j/jslee/QC8

#newIndex=${RUN}_${JobId}_c${column}_r${row}_dx${dx}_dy${dy}_dz${dz}_rx${rx}_ry${ry}_rz${rz}
newIndex=${RUN}_${JobId}
runFile=runGEMCosmicStand_sim_temp_${newIndex}.py

#geomDir=$CMSSW_BASE/src/Geometry/MuonCommonData/data/GEMQC8
#geomFile=gem11L_c${column}_r${row}.xml

#newGeomDir=$newIndex
#mkdir $geomDir/$newGeomDir
#cp $geomDir/*.xml $geomDir/$newGeomDir

cd $srcDir
eval `scramv1 runtime -sh`
#eval `scramv1 runtime -csh`
cd $workDir

sed -e 's/"runNum",1,/"runNum",'${RUN}',/' \
    -e 's/"eventsPerJob",10000,/"eventsPerJob",'${NeventsPerJob}',/' \
    -e 's/"idxJob","-1",/"idxJob","'${JobId}'",/' \
    -e 's#file:#file:'${tempDir}'/#' \
    -e 's#temp_#'${tempDir}'/temp_#' ${srcFile} > ${runFile}
#    -e 's/out_reco_MC.root/out_reco_MC_'${RUN}'_'${JobId}'.root/' \
#    -e 's/AddTrueDx/'${trueDx}'/' \
#    -e 's/AddTrueDy/'${trueDy}'/' \
#    -e 's/AddTrueDz/'${trueDz}'/' \
#    -e 's/AddTrueRz/'${trueRz}'/' ${srcFile} > ${runFile}

chmod 777 ${runFile}
cmsRun ${runFile}
rm ${tempDir}/out_reco_MC_${RUN}_${JobId}.root
mv ${tempDir}/temp_out_reco_MC_${RUN}_${JobId}.root ${resultDir}
rm ${runFile}
#rm -rf $geomDir/$newGeomDir

#sleep 5
fileName=`find . -type f -name 'core.*'`
for file in $fileName
do
  rm $file
done
LSFjobName=`find . -type d -name 'LSFJOB_*'`
for file in $LSFjobName
do
  rm -rf ${file}
done


