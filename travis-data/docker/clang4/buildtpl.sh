#!/bin/sh
set -e
set -x

git clone https://github.com/GEOSX/thirdPartyLibs.git
cd thirdPartyLibs
git submodule update --init --recursive
git checkout bugfix/siloBuild
python scripts/config-build.py -hc host-configs/default.cmake -bt Release
cd build-default-release
make
cd ..
git submodule deinit .
rm -rf build-default-release
