message( "Entering /src/externalComponents/CMakeLists.txt")



set(externalComponentsList "" )

if( ENABLE_PVT_PACKAGE )
  add_subdirectory( PVTPackage/PVTPackage/source )
  set( externalComponentsList ${externalComponentsList} PVTPackage )
endif()

file(GLOB children RELATIVE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/*)
foreach(child ${children})
  if(IS_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/${child})
    if( NOT ( ${child} STREQUAL "cxx-utilities" OR
              ${child} STREQUAL "PVTPackage" ) )
      add_subdirectory(${child})
      set(externalComponentsList ${externalComponentsList} ${child} )  
    endif()
  endif()
endforeach()

if( ${BUILD_SHARED_LIBS} )
  set( externalComponentsLinkList ${externalComponentsList} PARENT_SCOPE )
else()
  set(externalComponentsLinkList "" PARENT_SCOPE )
  foreach(component ${externalComponentsList})
    set( externalComponentsLinkList ${externalComponentsLinkList} ${GEOSX_LINK_PREPEND_FLAG} ${component} ${GEOSX_LINK_POSTPEND_FLAG} PARENT_SCOPE ) 
  endforeach()
endif()



message( "Leaving /src/externalComponents/CMakeLists.txt" )