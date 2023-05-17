set(CONFIG_NAME "sherlock-gcc8-cuda10-volta" CACHE PATH "") 

set(GCC_ROOT "/share/software/user/open/gcc/8.1.0" CACHE PATH "")
set(MPI_ROOT "/share/software/user/open/openmpi/4.0.3" CACHE PATH "")

set(ENABLE_CUDA ON CACHE BOOL "" FORCE)
set(ENABLE_HYPRE_CUDA ON CACHE BOOL "" FORCE)
set(CUDA_ARCH "sm_70" CACHE STRING "") 
set(CMAKE_CUDA_ARCHITECTURES "70" CACHE STRING "")
set(CUDA_TOOLKIT_ROOT_DIR "/share/software/user/open/cuda/10.2.89" CACHE STRING "") 

include(/home/groups/tchelepi/geosx/host-configs/sherlock-base.cmake)

