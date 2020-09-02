#!/bin/bash

mkdir ./tmp_anaconda
cd ./tmp_anaconda

wget --no-check-certificate http://repo.continuum.io/archive/Anaconda2-4.1.1-Linux-x86_64.sh
bash ./Anaconda2-4.1.1-Linux-x86_64.sh
cd ..
rm -r ./tmp_anaconda