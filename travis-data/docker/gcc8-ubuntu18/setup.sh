#!/bin/sh
set -e
set -x


sudo apt-get -qq install -y --no-install-recommends \
     g++-${gccver} \
    && sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-${gccver} 100 \
    && sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-${gccver} 100 \
    && sudo update-alternatives --install /usr/bin/cc cc /usr/bin/gcc-${gccver} 100 \
    && sudo update-alternatives --install /usr/bin/c++ c++ /usr/bin/g++-${gccver} 100
sudo apt-get --assume-yes install gfortran-8
cd /usr/bin
sudo ln -s gfortran-8 gfortran
sudo ln -s g77
sudo ln -s g90
cd ~

#sudo apt-get --assume-yes install -f  mpich
sudo apt-get --assume-yes install openmpi-bin libopenmpi-dev
