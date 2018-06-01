#!/bin/sh
set -e
set -x


apt-get -qq update
apt-get -qq install -y --no-install-recommends wget
wget "http://keyserver.ubuntu.com/pks/lookup?op=get&search=0x60C317803A41BA51845E371A1E9377A2BA9EF27F" -O out && apt-key add out && rm out
echo deb http://ppa.launchpad.net/ubuntu-toolchain-r/test/ubuntu xenial main >> /etc/apt/sources.list
echo deb-src http://ppa.launchpad.net/ubuntu-toolchain-r/test/ubuntu xenial main >> /etc/apt/sources.list
apt-get -qq update
apt-get -qq install -y --no-install-recommends python-dev build-essential sudo git dh-autoreconf ninja-build ca-certificates libtbb-dev zlib1g-dev libblas-dev liblapack-dev python-sphinx uncrustify

#wget -q --no-check-certificate https://cmake.org/files/v3.11/cmake-3.11.1-Linux-x86_64.tar.gz
#tar -xzf cmake-3.11.1-Linux-x86_64.tar.gz
#cp -fR cmake-3.11.1-Linux-x86_64/* /usr
#rm -rf cmake-3.11.1-Linux-x86_64
#rm cmake-3.11.1-Linux-x86_64.tar.gz

wget -q --no-check-certificate https://cmake.org/files/v3.8/cmake-3.8.2-Linux-x86_64.tar.gz
tar -xzf cmake-3.8.2-Linux-x86_64.tar.gz
cp -fR cmake-3.8.2-Linux-x86_64/* /usr
rm -rf cmake-3.8.2-Linux-x86_64
rm cmake-3.8.2-Linux-x86_64.tar.gz

wget -q --no-check-certificate https://bootstrap.pypa.io/get-pip.py
python get-pip.py
rm get-pip.py
pip install -q -U pip
useradd -ms /bin/bash geosx
printf "geosx:geosx" | chpasswd
adduser geosx sudo
printf "geosx ALL= NOPASSWD: ALL\\n" >> /etc/sudoers

sudo mkdir /home/geosx/geosx_repo
