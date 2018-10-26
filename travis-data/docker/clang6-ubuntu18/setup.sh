#!/bin/sh
set -e
set -x

sudo apt-get --assume-yes install apt-utils
sudo apt-get --assume-yes install clang-6.0

sudo apt-get --assume-yes install gfortran-7
cd /usr/bin
sudo ln -s gfortran-7 gfortran
sudo ln -s g77
sudo ln -s g90
cd ~

sudo apt-get --assume-yes install openmpi-bin libopenmpi-dev
mpicxx --version