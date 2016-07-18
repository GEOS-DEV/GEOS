message("Processing hc-defaults.cmake")
set(ENABLE_FORTRAN OFF CACHE BOOL  "Disables Fortran support")

set(BLT_CXX_STD "c++14" CACHE STRING "Version of C++ standard")

#option(ENABLE_OPENMP     "Enables OpenMP compiler support" ON)

message("Leaving hc-defaults.cmake")
