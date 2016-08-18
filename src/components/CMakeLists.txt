blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -I${CMAKE_SOURCE_DIR}/components/core/src)
blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -L${PYTHON_DIR}/lib)

blt_add_component(COMPONENT_NAME core DEFAULT_STATE ON)
blt_add_component(COMPONENT_NAME cxx-utilities DEFAULT_STATE ON)
#blt_add_component(COMPONENT_NAME newComponentTemplate DEFAULT_STATE ON)

file(GLOB children RELATIVE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/*)
message("children = ${children}")


set(componentlist geosx_core cxx-utilities)
foreach(child ${children})
  if(IS_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/${child})
    if( NOT ( ${child} STREQUAL "core" OR ${child} STREQUAL "docs"  OR ${child} STREQUAL "cxx-utilities" ) )
#      list(APPEND dirlist ${child})
      blt_add_component(COMPONENT_NAME ${child} DEFAULT_STATE ON)    
      set(componentlist ${componentlist} ${child} )  
    endif()
  endif()
endforeach()

message("componentlist = ${componentlist}")
set(componentlistlink "")

foreach(component ${componentlist})
  set( componentlistlink ${componentlistlink} ${GEOSX_LINK_PREPEND_FLAG} ${component} ${GEOSX_LINK_POSTPEND_FLAG} ) 
endforeach()
message("componentlistlink = ${componentlistlink}")


blt_add_executable(NAME geosx
                   SOURCES core/src/main/main.cpp
                   DEPENDS_ON ${componentlistlink}
                   SYSTEM_DEPENDS_ON slic sidre conduit hdf5 raja python_interp)

#blt_add_executable(NAME geosx
#                   SOURCES core/src/main/main.cpp
#                   DEPENDS_ON ${GEOSX_LINK_PREPEND_FLAG} geosx_core ${GEOSX_LINK_POSTPEND_FLAG} 
#                              ${GEOSX_LINK_PREPEND_FLAG} newComponentTemplate ${GEOSX_LINK_POSTPEND_FLAG}
#                              ${GEOSX_LINK_PREPEND_FLAG} cxx-utilities ${GEOSX_LINK_POSTPEND_FLAG}
#                   SYSTEM_DEPENDS_ON slic sidre conduit hdf5 raja )
                   
if (ENABLE_DOCS)
  add_subdirectory(docs)
endif()
