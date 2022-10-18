

#
# Specify all headers
#
set( schema_headers
     schemaUtilities.hpp )

#
# Specify all sources
#
set( schema_sources
     schemaUtilities.cpp )

set( dependencyList ${parallelDeps} dataRepository )

if( ENABLE_CUDA_NVTOOLSEXT )
  set( dependencyList ${dependencyList} CUDA::nvToolsExt )
endif()

blt_add_library( NAME       schema
                 SOURCES    ${schema_sources}
                 HEADERS    ${schema_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOSX_BUILD_OBJ_LIBS}
               )


target_include_directories( schema PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX schema )
