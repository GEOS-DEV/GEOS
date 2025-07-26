#!/bin/bash 

for F in $* ; do 
    echo "## Working with $F "
    $ECHO sed -i \
	-e 's/gcc/$ENV{GEOSX_CC}/g' \
	-e 's/g++/$ENV{GEOSX_CXX}/g' \
	-e 's/gfortran/$ENV{GEOSX_FORT}/g' \
	-e 's/mpiicc/$ENV{GEOSX_MPICC}/g' \
	-e 's/mpiicpc/$ENV{GEOSX_MPICXX}/g' \
	-e 's/mpiifort/$ENV{GEOSX_MPIFORT}/g' \
	-e 's/mpirun/$ENV{GEOSX_MPIRUN}/g' $F
done

 
