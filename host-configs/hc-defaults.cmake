message("Processing hc-defaults.cmake/n")
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


option(ENABLE_CONTAINERARRAY_RETURN_PTR     "Enables ViewWrapper to return pointers instead of references" ON )
message( "ENABLE_CONTAINERARRAY_RETURN_PTR = ${ENABLE_CONTAINERARRAY_RETURN_PTR}" )




option( ENABLE_CHAI "Enables CHAI" ON )
option( ENABLE_RAJA "Enables RAJA" ON )
#option( ENABLE_CALIPER "Enables CALIPER" OFF )
option( ENABLE_FPARSER "Enables FPARSER" OFF )



message("Leaving hc-defaults.cmake\n")
