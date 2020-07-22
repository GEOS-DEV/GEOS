# Delete me

# wget https://www.python.org/ftp/python/3.8.5/Python-3.8.5.tgz
# https://bugs.python.org/issue39697

mkdir -p /usr/WS1/GEOS/GEOSX/python/quartz-clang@10.0.0/lib
./configure --prefix=/usr/WS1/GEOS/GEOSX/python/quartz-clang@10.0.0 CXX=/usr/tce/packages/clang/clang-10.0.0/bin/clang++ CC=/usr/tce/packages/clang/clang-10.0.0/bin/clang --enable-shared LDFLAGS="-Wl,-rpath /usr/WS1/GEOS/GEOSX/python/quartz-clang@10.0.0/lib"
make -j && make install
/usr/WS1/GEOS/GEOSX/python/quartz-clang@10.0.0/bin/python3 -m pip install scipy
MPICC=/usr/tce/packages/mvapich2/mvapich2-2.3-clang-10.0.0/bin/mpicc /usr/WS1/GEOS/GEOSX/python/quartz-clang@10.0.0/bin/python3 -m pip install mpi4py


mkdir -p /usr/WS1/GEOS/GEOSX/python/quartz-gcc@8.1.0/lib
./configure --prefix=/usr/WS1/GEOS/GEOSX/python/quartz-gcc@8.1.0 CXX=/usr/tce/packages/gcc/gcc-8.1.0/bin/g++ CC=/usr/tce/packages/gcc/gcc-8.1.0/bin/gcc --enable-shared LDFLAGS="-Wl,-rpath /usr/WS1/GEOS/GEOSX/python/quartz-gcc@8.1.0/lib"
make -j && make install
/usr/WS1/GEOS/GEOSX/python/quartz-gcc@8.1.0/bin/python3 -m pip install scipy
MPICC=/usr/tce/packages/mvapich2/mvapich2-2.3-gcc-8.1.0/bin/mpicc /usr/WS1/GEOS/GEOSX/python/quartz-gcc@8.1.0/bin/python3 -m pip install mpi4py

