
#
# Specify all headers
#
set(meshUtilities_headers
    MeshManager.hpp
    MeshGeneratorBase.hpp
    InternalMeshGenerator.hpp
    MeshUtilities.hpp
    SimpleGeometricObjects/GeometricObjectManager.hpp
    SimpleGeometricObjects/SimpleGeometricObjectBase.hpp
    SimpleGeometricObjects/Box.hpp
    StructuredGridUtilities.hpp
   )

#
# Specify all sources
#
set(meshUtilities_sources
    MeshManager.cpp
    MeshGeneratorBase.cpp
    InternalMeshGenerator.cpp
    MeshUtilities.cpp
    SimpleGeometricObjects/GeometricObjectManager.cpp
    SimpleGeometricObjects/SimpleGeometricObjectBase.cpp
    SimpleGeometricObjects/Box.cpp
   )
    
    
#set( dependencyList fileIO common  )
#set( dependencyList2 trilinos )
#foreach( lib ${dependencyList} )
#    set( dependencyList2 ${dependencyList2} ${GEOSX_LINK_PREPEND_FLAG} ${lib} ${GEOSX_LINK_POSTPEND_FLAG} )
#endforeach()
#message( "dependencyList2 = ${dependencyList2}" )
#
#blt_add_library( NAME                  meshUtilities
#                 SOURCES               ${meshUtilities_sources}
#                 HEADERS               ${meshUtilities_headers}
#                 DEPENDS_ON            ${dependencyList2}
#                 SHARED FALSE
#               )

blt_add_library( NAME                  meshUtilities
                 SOURCES               ${meshUtilities_sources}
                 HEADERS               ${meshUtilities_headers}
                 DEPENDS_ON            fileIO common trilinos
                 SHARED FALSE
                 OBJECT TRUE
               )
                              
target_include_directories( meshUtilities PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src)
geosx_add_code_checks(PREFIX meshUtilities )