
set(ATK_ROOT "/Users/rrsettgast/Codes/geosx/asctoolkit" CACHE PATH "")
set(CONFIG_NAME "rockhopper.local-darwin-x86_64-gcc\@mp5" CACHE PATH "") 
set(ATK_DIR "${ATK_ROOT}/install-${CONFIG_NAME}-debug" CACHE PATH "")
set(RAJA_DIR "/Users/rrsettgast/Codes/geosx/RAJA/install-clang-release" CACHE PATH "")

message("ATK_DIR=${ATK_DIR}")
include("${CMAKE_CURRENT_LIST_DIR}/hc-defaults.cmake")
include("${ATK_ROOT}/uberenv_libs.gcc5/${CONFIG_NAME}.cmake")


set(GEOSX_LINK_PREPEND_FLAG "-Wl,-force_load" CACHE PATH "" FORCE)
set(GEOSX_LINK_POSTPEND_FLAG "" CACHE PATH "" FORCE)


#set(ENABLE_MPI ON CACHE PATH "")
#set(MPI_C_COMPILER "/opt/local/bin/mpicc-openmpi-gcc5" CACHE PATH "")
#set(MPI_CXX_COMPILER "/opt/local/bin/mpicxx-openmpi-gcc5" CACHE PATH "")
#set(MPI_Fortran_COMPILER "/opt/local/bin/mpifort-openmpi-gcc5" CACHE PATH "")
