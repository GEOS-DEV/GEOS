blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -I${CMAKE_SOURCE_DIR}/components/core/src)

blt_add_component(COMPONENT_NAME core DEFAULT_STATE ON)
blt_add_component(COMPONENT_NAME newComponentTemplate DEFAULT_STATE ON)


blt_add_executable(NAME geosx
                   SOURCES core/src/main/main.cpp
                   DEPENDS_ON ${GEOSX_LINK_PREPEND_FLAG} geosx_core ${GEOSX_LINK_POSTPEND_FLAG} 
                              ${GEOSX_LINK_PREPEND_FLAG} componentName ${GEOSX_LINK_POSTPEND_FLAG}
                   SYSTEM_DEPENDS_ON slic sidre conduit hdf5 raja )
                   
if (ENABLE_DOCS)
  add_subdirectory(docs)
endif()
                   