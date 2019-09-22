include(${CMAKE_CURRENT_LIST_DIR}/../host-configs/lassen-clang@upstream-NoMPI.cmake)

set(ENABLE_MPI OFF CACHE BOOL "" FORCE)
set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC OFF CACHE BOOL "" FORCE )
