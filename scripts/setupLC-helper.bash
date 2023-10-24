#!/bin/bash

## Builds GEOS for a specific system and host config.
## Usage ./setupLC-helper.bash machine compiler commandToGetANode [extra arguments to config-build ]
MACHINE=$1
COMPILER=$2
GET_A_NODE=$3

## Eat up the command line arguments so the rest can be forwarded to config-build.
shift
shift
shift

CONFIG=$MACHINE-$COMPILER
LOG_FILE=$CONFIG.log
HOST_CONFIG=host-configs/LLNL/$CONFIG.cmake

echo "Building GEOS on $MACHINE for $HOST_CONFIG. Progress will be written to $LOG_FILE."

ssh $MACHINE -t "
. /etc/profile  &&
cd $PWD &&
module load cmake/3.23.1 &&
python3 scripts/config-build.py -hc $HOST_CONFIG -bt Release $@ &&
cd build-$CONFIG-release &&
$GET_A_NODE make -j 36 &&
exit" > $LOG_FILE 2>&1
