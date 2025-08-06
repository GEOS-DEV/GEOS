#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --partition=parallel
#SBATCH -A rhurley6
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0
#SBATCH --export=ALL

echo "Starting GEOS compile by..."
start=`date +%s`
cd build-rockfish-clang-release/
#geos-env
make -j48 geosx
end=`date +%s`
runtime=$((end-start))
echo "Script execution time: $runtime seconds"
echo "GEOS compiled successfully, YAYAY!"
