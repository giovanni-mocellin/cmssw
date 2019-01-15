#!/bin/bash
#cmsenv
if [ "$1" == "all" ] ; then
    cd ../../../
else
    cd ../../
fi
pwd
echo
scram b -j 4
