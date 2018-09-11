#!/bin/sh
set -e
set -x

sudo apt-get --assume-yes install apt-utils
sudo apt-get --assume-yes install clang-5.0
sudo apt-get --assume-yes install openmpi-bin libopenmpi-dev
mpicxx --version