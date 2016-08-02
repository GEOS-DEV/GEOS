#!/bin/bash
. /usr/local/tools/dotkit/init.sh
use gcc-4.9.3p
use clang-3.8.1

g++      -std=c++14 -O3  -o gcc49.x    main.cpp
clang++  -std=c++14 -O3  -o clang38.x  main.cpp

echo gcc49
for i in `seq 1 10`;
do
    srun -n1 -p pdebug ./gcc49.x   200 400 200 100 2
done

echo clang38
for i in `seq 1 10`;
do
    srun -n1 -ppdebug ./clang38.x  200 400 200 100 2
done
