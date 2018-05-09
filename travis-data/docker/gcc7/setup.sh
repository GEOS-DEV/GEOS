#!/bin/sh
set -e
set -x


sudo apt-get -qq install -y --no-install-recommends \
     g++-${gccver} \
    && sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-${gccver} 100 \
    && sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-${gccver} 100 \
    && sudo update-alternatives --install /usr/bin/cc cc /usr/bin/gcc-${gccver} 100 \
    && sudo update-alternatives --install /usr/bin/c++ c++ /usr/bin/g++-${gccver} 100
#sudo apt-get --assume-yes install -f  mpich
sudo apt-get --assume-yes install openmpi-bin libopenmpi-dev
git clone https://github.com/GEOSX/thirdPartyLibs.git
cd thirdPartyLibs
git submodule init
git submodule update
python scripts/config-build.py -hc host-configs/default.cmake -bt Release
cd build-default-release
make hdf5
make
cd ..
git submodule deinit .
rm -rf build-default-release
