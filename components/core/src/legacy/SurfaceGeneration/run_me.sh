#!/usr/bin/bash

#ARGS = 2/3D point file, mean, stdev

export SCRIPT_PATH=../../scripts

perl $SCRIPT_PATH/points2minmax.pl $1 > mm

#create file1 input for apgen
cat mm | awk '/^0/{print $2}' > file1
cat mm | awk '/^1/{print $2}' >> file1
cat mm | awk '/^0/{print $3}' >> file1
cat mm | awk '/^1/{print $3}' >> file1
echo "5 5 1.2 0 5 $2 $3 1.3" >> file1

cat mm | awk '/^0/{print $2}' > file1b
cat mm | awk '/^1/{print $2}' >> file1b
cat mm | awk '/^2/{print $2}' >> file1b
cat mm | awk '/^0/{print $3}' >> file1b
cat mm | awk '/^1/{print $3}' >> file1b
cat mm | awk '/^2/{print $3}' >> file1b
echo "5 5 5 1.2 0 3 $2 $3 1.3" >> file1b

#create file2 input for apgen
cat $1 | wc | awk '{print $1}' > file2
cat $1 | awk '{print $1" "$2}' >> file2

cat $1 | wc | awk '{print $1}' > file2b
cat $1 | awk '{print $1" "$2" "$3}' >> file2b


#make aperture distribution
echo "X Y Z f" > aps.3D
echo "X Y Z f" > aps3.3D
./apgen.x file1 file2 > aps
./volgen file1b file2b > aps3
cat aps >> aps.3D
cat aps3 >> aps3.3D

rm mm
