#!/bin/bash
set -e

# See https://bugs.python.org/issue39697

# # For gcc8
# _CC=/usr/tce/packages/gcc/gcc-8.1.0/bin/gcc
# _CXX=/usr/tce/packages/gcc/gcc-8.1.0/bin/g++
# _MPICC=/usr/tce/packages/mvapich2/mvapich2-2.3-gcc-8.1.0/bin/mpicc
# _OPENMPCFLAGS="-fopenmp"
# INSTALL_DIR=/usr/gapps/GEOSX/python/quartz-gcc@8.1.0

# For clang10
_CC=/usr/tce/packages/clang/clang-10.0.0/bin/clang
_CXX=/usr/tce/packages/clang/clang-10.0.0/bin/clang++
_MPICC=/usr/tce/packages/mvapich2/mvapich2-2.3-clang-10.0.0/bin/mpicc
_OPENMPCFLAGS="-fopenmp=libomp"
INSTALL_DIR=/usr/gapps/GEOSX/python/quartz-clang@10.0.0


_CFLAGS="-O3 -DNDEBUG -march=native -mtune=native"


# # Create the install directory. The lib subdirectory must exist for the compilation to succeed.
# if [ -d $INSTALL_DIR ]; then
#   rm -rf $INSTALL_DIR
# fi

# mkdir -p $INSTALL_DIR/lib

# # Build python
# if [ ! -f Python-3.8.5.tgz ]; then
#   wget https://www.python.org/ftp/python/3.8.5/Python-3.8.5.tgz
# fi

# if [ -d Python-3.8.5 ]; then
#   rm -rf Python-3.8.5
# fi

# tar -xzvf Python-3.8.5.tgz
# cd Python-3.8.5

# ./configure --prefix=$INSTALL_DIR CXX=$_CXX CC=$_CC CFLAGS="$_CFLAGS" --enable-shared LDFLAGS="-Wl,-rpath $INSTALL_DIR/lib -lpthread"
# make -j && make install

# cd ..

# # Pip install cython
# $INSTALL_DIR/bin/python3 -m pip install cython

# # Build numpy
# if [ ! -d "numpy" ]; then
#   mkdir numpy
# fi

# cd numpy

# if [ ! -f "v1.19.1.tar.gz" ]; then
#   wget https://github.com/numpy/numpy/archive/v1.19.1.tar.gz
# fi
# if [ -d "numpy-1.19.1" ]; then
#   rm -rf numpy-1.19.1
# fi

# tar -xzvf v1.19.1.tar.gz
# cd numpy-1.19.1

# echo "
# [mkl]
# library_dirs = /usr/tce/packages/mkl/mkl-2019.0/lib
# include_dirs = /usr/tce/packages/mkl/mkl-2019.0/include
# mkl_libs = mkl_intel_lp64,mkl_gnu_thread,mkl_core
# lapack_libs =
# " > site.cfg

# CFLAGS="$_CFLAGS $_OPENMPCFLAGS" LDFLAGS="/usr/tce/packages/mkl/mkl-2019.0/lib/libmkl_intel_lp64.so /usr/tce/packages/mkl/mkl-2019.0/lib/libmkl_gnu_thread.so /usr/tce/packages/mkl/mkl-2019.0/lib/libmkl_core.so /usr/tce/packages/mkl/mkl-2019.0/lib/libmkl_avx2.so /usr/tce/packages/mkl/mkl-2019.0/lib/libmkl_def.so -Wl,-rpath /usr/tce/packages/mkl/mkl-2019.0/lib" $INSTALL_DIR/bin/python3 setup.py build -j 36 install

# cd ../..

# Build scipy (We might be able to just pip-install scipy and have everything work but it doesn't link to MKL)
if [ ! -d "scipy" ]; then
  mkdir scipy
fi

cd scipy

if [ ! -f v1.5.2.tar.gz ]; then
  wget https://github.com/scipy/scipy/archive/v1.5.2.tar.gz
fi
if [ -d scipy-1.5.2 ]; then
  rm -rf scipy-1.5.2
fi

if [ -f scipy-1.5.2.patch ]; then
  rm scipy-1.5.2.patch
fi

# The following patch is needed for clang, it seems to be a bug in scipy.
echo "
diff --git a/scipy/special/_faddeeva.cxx b/scipy/special/_faddeeva.cxx
index 9134516..159122c 100644
--- a/scipy/special/_faddeeva.cxx
+++ b/scipy/special/_faddeeva.cxx
@@ -130,7 +130,7 @@ double faddeeva_voigt_profile(double x, double sigma, double gamma)
 
     if(sigma == 0){
         if (gamma == 0){
-            if (isnan(x))
+            if (std::isnan(x))
                 return x;
             if (x == 0)
                 return NPY_INFINITY;
" > scipy-1.5.2.patch

tar -xzvf v1.5.2.tar.gz
cd scipy-1.5.2
git apply --reject ../scipy-1.5.2.patch
$INSTALL_DIR/bin/python3 setup.py build -j 36 install
cd ../..

# Now pip install other packages
MPICC=$_MPICC $INSTALL_DIR/bin/python3 -m pip install mpi4py
