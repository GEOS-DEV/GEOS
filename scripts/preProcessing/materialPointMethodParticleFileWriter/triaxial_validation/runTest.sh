rm CDTriaxial_results_updated.png

cd ..
echo "y" | ./runClean_sghosh29.sh

cd /data1/sghosh29/Working_MPM_LLNL/testGEOS
cd CDtriaxial

mpirun -np 27 /data1/sghosh29/Working_MPM_LLNL/GEOS/build-system76-pc-clang-release/bin/geosx -i mpm_CDtriaxial.xml -x 3 -y 3 -z 3

cd /data1/sghosh29/Working_MPM_LLNL/GEOS/scripts/preProcessing/materialPointMethodParticleFileWriter/triaxial_validation

python plotReactions_CDTriaxial.py
