set(CONFIG_NAME "sherlock-gcc10-cuda11-ampere" CACHE PATH "") 

set(GCC_ROOT "/share/software/user/open/gcc/10.1.0" CACHE PATH "")
set(MPI_ROOT "/share/software/user/open/openmpi/4.1.0" CACHE PATH "")

set(ENABLE_CUDA ON CACHE BOOL "" FORCE)
set(ENABLE_HYPRE_CUDA ON CACHE BOOL "" FORCE)
set(CUDA_ARCH "sm_80" CACHE STRING "") 
set(CMAKE_CUDA_ARCHITECTURES "80" CACHE STRING "")
set(CUDA_TOOLKIT_ROOT_DIR "/share/software/user/open/cuda/11.2.0" CACHE STRING "") 

include(/home/groups/tchelepi/geosx/host-configs/sherlock-base.cmake)

