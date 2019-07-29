message( "Entering /src/coreComponents/linearAlgebraInterface/CMakeLists.txt")

add_subdirectory(src)
if( ENABLE_TESTS )
  add_subdirectory(tests)
endif()

message( "Leaving /src/coreComponents/linearAlgebraInterface/CMakeLists.txt")
