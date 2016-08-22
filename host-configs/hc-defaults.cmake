message("Processing hc-defaults.cmake")
set(ENABLE_FORTRAN OFF CACHE BOOL  "Disables Fortran support")

set(BLT_CXX_STD "c++14" CACHE STRING "Version of C++ standard")

#message("CMAKE_HOST_SYSTEM_NAME = ${CMAKE_HOST_SYSTEM_NAME}")
#message("CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
#message("CMAKE_HOST_APPLE = ${CMAKE_HOST_APPLE}")
if( CMAKE_HOST_APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang" )
  option(ENABLE_OPENMP     "Enables OpenMP compiler support" OFF FORCE)  
else()
  option(ENABLE_OPENMP     "Enables OpenMP compiler support" ON FORCE)
endif()

message("Leaving hc-defaults.cmake")
