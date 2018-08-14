#!/bin/sh
set -e
set -x


wget -q --no-check-certificate http://releases.llvm.org/5.0.2/${llvmtar}${tarext} \
    && tar xf ${llvmtar}${tarext} \
    && sudo cp -fR ${llvmtar}/* /usr \
    && rm -rf ${llvmtar} \
    && rm ${llvmtar}${tarext}

sudo apt-get --assume-yes install openmpi-bin libopenmpi-dev
git clone https://github.com/GEOSX/thirdPartyLibs.git
cd thirdPartyLibs
git submodule update --init --recursive
python scripts/config-build.py -hc host-configs/default.cmake -bt Release
cd build-default-release
make
cd ..
git submodule deinit .
rm -rf build-default-release
