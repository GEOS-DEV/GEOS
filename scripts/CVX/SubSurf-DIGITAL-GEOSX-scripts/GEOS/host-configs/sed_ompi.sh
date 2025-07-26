#!/bin/bash 

for F in $* ; do 
    echo "## Working with $F "
    $ECHO sed -i \
	-e 's/gcc/$ENV{GEOSX_CC}/g' \
	-e 's/g++/$ENV{GEOSX_CXX}/g' \
	-e 's/gfortran/$ENV{GEOSX_FORT}/g' \
	-e 's/mpicc/$ENV{GEOSX_MPICC}/g' \
	-e 's/mpicxx/$ENV{GEOSX_MPICXX}/g' \
	-e 's/mpif90/$ENV{GEOSX_MPIFORT}/g' \
	-e 's/mpirun/$ENV{GEOSX_MPIRUN}/g' $F
done
