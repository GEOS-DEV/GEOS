set(ENABLE_FORTRAN OFF CACHE BOOL  "Enables C++11 language support")

set(BLT_CXX_STD "c++14" CACHE STRING "Version of C++ standard")

set( OpenMP_CXX_FLAGS "-fopenmp" CACHE STRING "" FORCE)
set( OpenMP_C_FLAGS   "-fopenmp" CACHE STRING "" FORCE)
message("OpenMP_C_FLAGS=${OpenMP_C_FLAGS}")
message("OpenMP_CXX_FLAGS=${OpenMP_CXX_FLAGS}")
option(ENABLE_OPENMP     "Enables OpenMP compiler support" ON)
