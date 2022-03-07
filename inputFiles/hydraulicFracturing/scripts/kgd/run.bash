rm serial.log
rm parallel.log
#srun -n1 /g/g90/cusini1/geosx/GEOSX/build-toss_3_x86_64_ib-clang@6.0.0-debug/bin/geosx -i /g/g90/cusini1/geosx/GEOSX/examples/hydraulicFracturing/KGD/KGD_rotated.xml -x 1 -y 1 -z 1 > serial.log 2>&1
#srun -n2 /g/g90/cusini1/geosx/GEOSX/build-toss_3_x86_64_ib-clang@6.0.0-debug/bin/geosx -i /g/g90/cusini1/geosx/GEOSX/examples/hydraulicFracturing/KGD/KGD_rotated.xml -x 2 -y 1 -z 1 > parallel.log 2>&1

OMP_NUM_THREADS=1 ../../../build-darwin-clang-debug/bin/geosx -i KGD_rotated.xml > serial.log 2>&1
OMP_NUM_THREADS=1 mpirun -np 2 ../../../build-darwin-clang-debug/bin/geosx -i KGD_rotated.xml -x 2 -y 1 -z 1 > parallel.log 2>&1
