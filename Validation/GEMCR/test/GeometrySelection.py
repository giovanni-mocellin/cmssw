#!/usr/bin/env python
# Copyright [2017] <III Phys. Inst A, RWTH Aachen University>

import sys
import os
import random
import datetime

from datetime import datetime
import optparse,os,time,cPickle,subprocess,shutil,sys
import logging
log = logging.getLogger( 'remote' )
import yaml

# Possible geometries: A , B , C , P5
GeomSelected = "A"

def getSamplesFromLocalDisk():
    filestring = "ls /home/home1/institut_3a/mocellin/tests/src/Geometry/MuonCommonData/data/cosmic1/GeomFiles/geom" + GeomSelected + "/gem11*"
    proc = subprocess.Popen([filestring], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    listOfSamples = out.split("\n")
    newList = []
    for item in listOfSamples:
        if not item:
            continue
        newList.append(item)
    return newList
    
def main():
	directory = "/home/home1/institut_3a/mocellin/tests/src/Geometry/MuonCommonData/data/cosmic1/GeomFiles/geom" + GeomSelected
	destination = "/home/home1/institut_3a/mocellin/tests/src/Geometry/MuonCommonData/data/cosmic1"
	print getSamplesFromLocalDisk()
	for item in getSamplesFromLocalDisk():
		name = item.split("/")[-1]
		print ("copy from /Geometry/MuonCommonData/data/cosmic1/GeomFiles/geom" + GeomSelected + "/" + name + " to " + "/Geometry/MuonCommonData/data/cosmic1/" + name)
		os.system("cp %s/%s %s/%s"%(directory,name,destination,name))
		
if __name__ == '__main__':
    main()
