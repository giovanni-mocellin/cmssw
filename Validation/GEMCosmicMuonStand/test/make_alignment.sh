#! /bin/sh

RUN=332
NeventsPerJob=10000
LastJob=100
StartJob=0

JobId=$StartJob

while [ $JobId -lt $LastJob ]
do
  echo NeventsPerJob${NeventsPerJob} run$RUN jobid$JobId
  #cp run_bsub.sh run_bsub_${RUN}_${JobId}.sh
  #cp run_bsub_alignment.sh run_bsub_alignment_${RUN}_${JobId}.sh

  #For batch job
  bsub -q 1nd run_bsub_alignment.sh ${NeventsPerJob} ${RUN} ${JobId}

  #For local job
  #./run_bsub_alignment.sh ${NeventsPerJob} ${RUN} ${JobId}

  sleep 3
  #rm run_bsub_alignment_${RUN}_${JobId}.sh
  JobId=`expr $JobId + 1`
done

