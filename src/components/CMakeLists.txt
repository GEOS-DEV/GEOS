blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -I${CMAKE_SOURCE_DIR}/components/core/src)
blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -isystem${CMAKE_SOURCE_DIR}/thirdparty)


# Python config
set(ENABLE_PYTHON OFF CACHE BOOL  "Disables Python Plugin")
set( extraComponentsListLink "")
if( ENABLE_PYTHON )
  if( PYTHON_DIR )
    set( extraComponentsListLink ${extraComponentsListLink} python_interp )
  	blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -DUSE_PYTHON=1)
		blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -L${PYTHON_DIR}/lib)
	else()
		MESSAGE(FATAL_ERROR "ENABLE_PYTHON is true, but PYTHON_DIR is not defined.")
	endif()
else()
	blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -DUSE_PYTHON=0)
endif()


blt_add_component(COMPONENT_NAME core DEFAULT_STATE ON)
add_subdirectory(cxx-utilities/src)

file(GLOB children RELATIVE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/*)
#message("children = ${children}")


set(componentlist geosx_core cxx-utilities)
foreach(child ${children})
  if(IS_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/${child})
    if( NOT ( ${child} STREQUAL "core" OR 
              ${child} STREQUAL "docs" OR 
              ${child} STREQUAL "cxx-utilities" OR 
              ${child} STREQUAL "ConstitutivePackage" OR 
              ${child} STREQUAL "MeshRefinementPackage" OR 
              ${child} STREQUAL "PhysicsSolverPackage2" OR 
              ${child} STREQUAL "ReservoirSimulationPackage" ) )
#      list(APPEND dirlist ${child})
      blt_add_component(COMPONENT_NAME ${child} DEFAULT_STATE ON)    
      set(componentlist ${componentlist} ${child} )  
    endif()
  endif()
endforeach()

#message("componentlist = ${componentlist}")
set(componentlistlink "")

foreach(component ${componentlist})
  set( componentlistlink ${componentlistlink} ${GEOSX_LINK_PREPEND_FLAG} ${component} ${GEOSX_LINK_POSTPEND_FLAG} ) 
endforeach()
#message("componentlistlink = ${componentlistlink}")


blt_add_executable(NAME geosx
                   SOURCES core/src/main/main.cpp
                   DEPENDS_ON ${componentlistlink} "-lz"
                   SYSTEM_DEPENDS_ON slic sidre conduit hdf5 raja chai pugixml ${extraComponentsListLink}) 

#set_property ( TARGET geosx APPEND_STRING PROPERTY LINK_FLAGS " -lz" )


                   
if (ENABLE_DOCS)
  add_subdirectory(docs)
endif()
