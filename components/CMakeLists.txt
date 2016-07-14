blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -I${CMAKE_SOURCE_DIR}/components/core/src)

blt_add_component(COMPONENT_NAME core DEFAULT_STATE ON)
blt_add_component(COMPONENT_NAME newComponentTemplate DEFAULT_STATE ON)

#add_subdirectory(core)
#add_subdirectory(newComponentTemplate)


blt_add_executable(NAME geosx
                   SOURCES core/src/main/main.cpp
                   DEPENDS_ON slic sidre conduit hdf5 geosx_core raja componentName)