
#
# Specify all headers
#
set(common_headers
    ${CMAKE_BINARY_DIR}/include/common/GeosxConfig.hpp
    DataTypes.hpp
    GeosxMacros.hpp
    TimingMacros.hpp
    Logger.hpp
   )

#
# Specify all sources
#
set(common_sources
    DataTypes.cpp
    Logger.cpp
   )

set( dependencyList cxx-utilities pugixml )

if( ENABLE_MPI )
  set( dependencyList ${dependencyList} mpi )
endif()

if( ATK_FOUND )
  set( dependencyList ${dependencyList} conduit )
endif()

if( ENABLE_CALIPER )
  set( dependencyList ${dependencyList} caliper )
endif()

#set( dependencyList2 mpi )
#foreach( lib ${dependencyList} )
#    set( dependencyList2 ${dependencyList2} ${GEOSX_LINK_PREPEND_FLAG} ${lib} ${GEOSX_LINK_POSTPEND_FLAG} )
#endforeach()
#message( "dependencyList2 = ${dependencyList2}" )


blt_add_library( NAME                  common
                 SOURCES               ${common_sources}
                 HEADERS               ${common_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
#                 SHARED TRUE
               )
               
               
               
target_include_directories( common PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)
               
geosx_add_code_checks(PREFIX common )

add_subdirectory( unitTests )
