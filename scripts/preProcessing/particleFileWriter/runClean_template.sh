#!/bin/bash

pfwLocation=$(pwd)
# Set full input file path where the input file is pfw_input_${fileName}.py 
inputFile="examples/bar/pfw_input_bar.py"
# Set full path to GEOS executable, which for example on quartz may be
pathToGEOS="${pfwLocation}/../../../build-quartz-clang@14-release/bin/geosx"
# Set user bank for this parallel job allocation
userBank="uco"

# This is where you want to run the simulation, which should be on a large 
#     parallel file system, for example lustre1 on LC. 
#     A sub-directory /fileName will be created within runLocation.
runLocation="/p/lustre1/$(whoami)/testRuns" 

fileName=$(echo "$inputFile" | grep -oP "pfw_input_\K[^.py]+")

# Create runLocation directory if it does not exist 
if [ ! -d "$runLocation" ]; then
	mkdir -p "$runLocation"
fi

# If fileName directory already exists delete contents, else mkdir
if [ -d "${runLocation}/${fileName}" ]; then
	rm -rf "${runLocation}/${fileName}/"
fi
# Create the /fileName subdirectory
mkdir $runLocation/$fileName/										

# Copy: input file, pfw, restart, geometry objects, and userDefs
cp $inputFile $runLocation/$fileName
cp particleFileWriter.py $runLocation/$fileName
cp pfw_check.py $runLocation/$fileName
cp pfw_geometryObjects.py $runLocation/$fileName

# Move to runLocation
cd $runLocation/$fileName

# Create userDefs file
echo "# -*- coding: utf-8 -*-" > userDefs_$(whoami).py
echo "geosPath=\"${pathToGEOS}\"" >> userDefs_$(whoami).py
echo "testRunDirectory=\"${runLocation}\"" >> userDefs_$(whoami).py
echo "defaultRunDirectory=\"${runLocation}\"" >> userDefs_$(whoami).py
echo "defaultBank=\"${userBank}\"" >> userDefs_$(whoami).py

# Run particleFileWriter.py:
python3 particleFileWriter.py pfw_input_$fileName