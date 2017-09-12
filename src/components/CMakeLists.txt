blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -I${CMAKE_SOURCE_DIR}/components/core/src)


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


add_subdirectory(core/src)
add_subdirectory(cxx-utilities/src)

file(GLOB children RELATIVE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/*)
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
      add_subdirectory(${child})
      set(componentlist ${componentlist} ${child} )  
    endif()
  endif()
endforeach()

set(componentlistlink "")

foreach(component ${componentlist})
  set( componentlistlink ${componentlistlink} ${GEOSX_LINK_PREPEND_FLAG} ${component} ${GEOSX_LINK_POSTPEND_FLAG} ) 
endforeach()
message("componentlistlink = ${componentlistlink}")

if( ENABLE_MATHPRESSO )
  set( extraComponentsListLink ${extraComponentsListLink} mathpresso )
endif()

blt_add_executable(NAME geosx
                   SOURCES core/src/main/main.cpp
                   DEPENDS_ON ${componentlistlink} ${extraComponentsListLink}) 

                   
if (ENABLE_DOCS)
  add_subdirectory(docs)
endif()
