
#
# Specify all headers
#
set(mesh_headers
    CellBlockManager.hpp
    CellBlock.hpp
    CellBlockSubRegion.hpp
    FaceManager.hpp
    NodeManager.hpp
    ElementRegion.hpp
    ElementRegionManager.hpp
    MeshBody.hpp
    MeshLevel.hpp    
  )


#
# Specify all sources
#
set(mesh_sources
    CellBlockManager.cpp
    CellBlock.cpp
    CellBlockSubRegion.cpp
    FaceManager.cpp
    NodeManager.cpp
    ElementRegion.cpp
    ElementRegionManager.cpp
    MeshBody.cpp
    MeshLevel.cpp    
   )
    

#set( dependencyList managers  )
#set( dependencyList2 "" )
#foreach( lib ${dependencyList} )
#    set( dependencyList2 ${dependencyList2} ${GEOSX_LINK_PREPEND_FLAG} ${lib} ${GEOSX_LINK_POSTPEND_FLAG} )
#endforeach()
#message( "dependencyList2 = ${dependencyList2}" )

blt_add_library( NAME                  mesh
                 SOURCES               ${mesh_sources}
                 HEADERS               ${mesh_headers}
                 DEPENDS_ON            managers
                 SHARED FALSE
               )
               
target_include_directories( mesh PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src)
geosx_add_code_checks(PREFIX mesh )