##################################
# uberenv host-config
##################################
# chaos_5_x86_64_ib-gcc@4.9.3
##################################

# cmake from uberenv
# cmake executable path: /usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/cmake-3.3.1-3gc4unffj5rqcq35gecg4wv3roecpldt/bin/cmake
#######
# using gcc@4.9.3 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/apps/gnu/4.9.3/bin/gcc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/apps/gnu/4.9.3/bin/g++" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN OFF CACHE PATH "")

set(ATK_DIR "/usr/gapps/GEOS/asctoolkit/install-cab-chaos_5_x86_64_ib-gcc@4.9.3-debug" CACHE PATH "")

#set(RAJA_DIR "/usr/gapps/GEOS/RAJA/install-gcc-4.9.3-release" CACHE PATH "")
#set(RAJA_DIR "/g/g15/settgast/workspace/Codes/RAJA/install-gcc-4.9.3-release" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "/usr/gapps/GEOS/asctoolkit/uberenv_libs/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/conduit-github-2016-05-18-xqbkgfstnxnbt43ptpb6d26iv5pvytyk" CACHE PATH "")

# Python path
set(PYTHON_DIR "/usr/gapps/GEOS/python/anaconda2" CACHE PATH "")

set(UNCRUSTIFY_EXECUTABLE "/usr/gapps/GEOS/asctoolkit/uberenv_libs/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/uncrustify-0.61-px2meiscmkbwcnmmom3qnlzdzmf2yx7x/bin/uncrustify" CACHE PATH "")



#######
# MPI - manually added these for now.
#######
set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/usr/local/tools/mvapich2-gnu-2.0/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/local/tools/mvapich2-gnu-2.0/bin/mpicxx" CACHE PATH "")
set(MPI_Fortran_COMPILER "/usr/local/tools/mvapich2-gnu-2.0/bin/mpif90" CACHE PATH "")

#####
# GCOV - manually added for now
#####
set(GCOV_PATH "/usr/apps/gnu/4.9.3/bin/gcov" CACHE PATH "")

include("${CMAKE_CURRENT_LIST_DIR}/hc-defaults.cmake")
set(GEOSX_LINK_PREPEND_FLAG  "-Wl,--whole-archive"    CACHE PATH "" FORCE)
set(GEOSX_LINK_POSTPEND_FLAG "-Wl,--no-whole-archive" CACHE PATH "" FORCE)
