
#set(ENABLE_FORTRAN OFF CACHE BOOL "Enables Fortran" FORCE)
#set(ENABLE_CXX11 OFF CACHE BOOL "Enables C++11 language support" FORCE)
#set(ENABLE_CXX14 ON CACHE BOOL "Enables C++14 language support" FORCE)

if(CXX_STD STREQUAL c++11)
    MESSAGE(FATAL_ERROR "c++11 is enabled. GEOSX requires c++14")
endif(CXX_STD STREQUAL c++11)

if(NOT CXX_STD STREQUAL c++14)
    MESSAGE(FATAL_ERROR "c++14 is NOT enabled. GEOSX requires c++14")
endif(NOT CXX_STD STREQUAL c++14)

blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -pedantic)
blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -Wabi)
blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -Wno-variadic-macros)