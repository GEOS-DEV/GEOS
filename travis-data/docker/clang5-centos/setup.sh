#!/bin/sh
set -e
set -x

sudo yum -y install centos-release-scl
sudo yum -y install llvm-toolset-7
sudo scl enable llvm-toolset-7 bash
sudo yum -y install openmpi-devel

export PATH=/opt/rh/llvm-toolset-7/root/usr/bin:/opt/rh/llvm-toolset-7/root/usr/sbin:/usr/lib64/openmpi/bin:$PATH
mpicxx --version