add_subdirectory(core)

add_subdirectory(newComponentTemplate)

blt_add_executable(NAME geosx
                   SOURCES core/src/main/main.cpp
                   DEPENDS_ON slic sidre conduit hdf5 geosx_core componentName)