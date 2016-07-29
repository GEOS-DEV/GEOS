#!/bin/bash

svn update

export ME=`whoami`
export SCRIPT_PATH=../../scripts/

export N=16
export MYCOMPILER=icc
if [ "$OSTYPE" = "darwin" ]; then
    export MYCOMPILER=mac
    export N=4
else
. /usr/local/tools/dotkit/init.sh
use ic-13.1.163
use mvapich2-intel-1.7
#  - or - 
use gcc-4.6.1
#mvapich2-gnu-1.7
fi

if [ "1" -eq "1" ]; then
make clean
make -j$N COMPILER=$MYCOMPILER DEBUG=0 2>&1 | tee c_opt
mv apgen apgen.x
fi

if [ "1" -eq "1" ]; then
make clean
make -j$N COMPILER=$MYCOMPILER DEBUG=1 2>&1 | tee c_debug
mv apgen apgen.dbg
fi

#g++ -O3 -I../ StatisticalDistributionBaseT.cpp -c
#g++ -O3 -I../ Aperturator.cpp -c
#g++ -O3 -I../ main.cpp -c
#g++ -O3 -I../ aperturatorlib.cpp -c
#g++ -O3 StatisticalDistributionBaseT.o Aperturator.o main.o -o apgen;
#ar -cvq libfractalsurface.a StatisticalDistributionBaseT.o Aperturator.o aperturatorlib.o
cp apgen.x $SCRIPT_PATH/apgen
