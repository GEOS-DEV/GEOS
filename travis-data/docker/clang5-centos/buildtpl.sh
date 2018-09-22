#!/bin/sh
set -e
set -x

sudo scl enable llvm-toolset-7 bash
git clone https://github.com/GEOSX/thirdPartyLibs.git
cd thirdPartyLibs
git checkout bugfix/addLumberjackAcom
git submodule update --init --recursive
python scripts/config-build.py -hc host-configs/default.cmake -bt Release -DNUM_PROC:STRING=2
cd build-default-release
make
cd ..
git submodule deinit .
rm -rf build-default-release
