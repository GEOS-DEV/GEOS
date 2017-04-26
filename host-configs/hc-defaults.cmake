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


option(ENABLE_CONTAINERARRAY_RETURN_PTR     "Enables ViewWrapper to return pointers instead of references" OFF FORCE)




set( thirdPartyLibs "")
option( ENABLE_CHAI ON CACHE BOOL "Enables CHAI" FORCE )
if( ENABLE_CHAI )
  set( thirdPartyLibs ${thirdPartyLibs} chai )
endif()

option( ENABLE_RAJA ON CACHE BOOL "Enables RAJA" FORCE )
if( ENABLE_RAJA )
  set( thirdPartyLibs ${thirdPartyLibs} raja )
endif()

option( ENABLE_CALIPER OFF CACHE BOOL "Enables CALIPER" FORCE )
if( ENABLE_CALIPER )
  set( thirdPartyLibs ${thirdPartyLibs} caliper )
endif()

option( ENABLE_FPARSER ON CACHE BOOL "Enables FPARSER" FORCE )
if( ENABLE_FPARSER )
  set( thirdPartyLibs ${thirdPartyLibs} fparser )
endif()


message("Leaving hc-defaults.cmake")
