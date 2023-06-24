site_name(HOST_NAME)
set(CONFIG_NAME "cvx-wsl" CACHE PATH "" FORCE)
message( "CONFIG_NAME=${CONFIG_NAME}" )

set(GCC_ROOT "/usr" CACHE PATH "")
include(${CMAKE_CURRENT_LIST_DIR}/cvx-gcc.cmake)

set(MPI_ROOT "/usr" CACHE PATH "" FORCE)
include(${CMAKE_CURRENT_LIST_DIR}/cvx-mpi.cmake)

include(${CMAKE_CURRENT_LIST_DIR}/cvx-base.cmake)
