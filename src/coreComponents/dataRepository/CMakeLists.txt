#
# Specify all headers
#
set( dataRepository_headers
     BufferOps.hpp
     BufferOpsDevice.hpp
     BufferOps_inline.hpp
     ConduitRestart.hpp
     DefaultValue.hpp
     ExecutableGroup.hpp
     Group.hpp
     HistoryDataSpec.hpp
     InputFlags.hpp
     KeyIndexT.hpp
     KeyNames.hpp
     MappedVector.hpp
     ObjectCatalog.hpp
     ReferenceWrapper.hpp
     RestartFlags.hpp
     Wrapper.hpp
     WrapperBase.hpp
     wrapperHelpers.hpp
     xmlWrapper.hpp
   )

#
# Specify all sources
#
set( dataRepository_sources
     BufferOpsDevice.cpp
     ConduitRestart.cpp
     ExecutableGroup.cpp
     Group.cpp
     WrapperBase.cpp
     xmlWrapper.cpp
   )

set( dependencyList codingUtilities )

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

if( ENABLE_PYGEOSX )
  set( dependencyList ${dependencyList} pylvarray )
  add_subdirectory( python )
endif()
 

blt_add_library( NAME       dataRepository
                 SOURCES    ${dataRepository_sources}
                 HEADERS    ${dataRepository_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOSX_BUILD_OBJ_LIBS}
               )

target_include_directories( dataRepository PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

geosx_add_code_checks(PREFIX dataRepository )

add_subdirectory( unitTests )
