
#
# Specify all headers
#
set(common_headers
    integer_conversion.hpp
    DataTypes.hpp
    Logger.hpp
    SortedArray.hpp )

#
# Specify all sources
#
set(common_sources
    Logger.cpp
    DataTypes.cpp )


set( dependencyList cxx-utilities mpi )

if( ATK_FOUND )
  set( dependencyList ${dependencyList} slic conduit lumberjack )
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
#                 SHARED TRUE
               )
               
               
               
target_include_directories( common PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)
               
geosx_add_code_checks(PREFIX common )

add_subdirectory( unitTests )
