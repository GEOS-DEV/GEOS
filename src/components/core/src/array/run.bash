#!/bin/bash

clang++        -std=c++14 -O3 -march=native   -I.. -o clang.x   main.cpp
clang++-mp-3.7 -std=c++14 -O3 -march=native   -I.. -o clang37.x main.cpp
g++-mp-5       -std=c++14 -O3 -march=westmere -I.. -o gcc5.x    main.cpp
g++-mp-4.9     -std=c++14 -O3 -march=westmere -I.. -o gcc49.x   main.cpp

echo clang-apple
for i in `seq 1 10`;
do
    ./clang.x   500 200 400 100
done

echo clang-37
for i in `seq 1 10`;
do
    ./clang37.x 500 100 400 100
done

echo gcc49
for i in `seq 1 10`;
do
    ./gcc49.x   500 100 400 100
done

echo gcc5
for i in `seq 1 10`;
do
    ./gcc5.x    500 100 400 100
done
