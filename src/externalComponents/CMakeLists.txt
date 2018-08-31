message( "Entering /src/externalComponents/CMakeLists.txt")


add_subdirectory(cxx-utilities/src)

if( ENABLE_PVT_PACKAGE )
    add_subdirectory(PVTPackage/PVTPackage/source)
endif()


message( "Leaving /src/externalComponents/CMakeLists.txt")
