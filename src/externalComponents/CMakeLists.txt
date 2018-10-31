message( "Entering src/externalComponents/CMakeLists.txt")



set(externalComponents "")

file(GLOB children RELATIVE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/*)
foreach(child ${children})
  if(IS_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/${child})
    if( EXISTS "${CMAKE_CURRENT_LIST_DIR}/${child}/CMakeLists.txt" AND ENABLE_${child} )
      add_subdirectory(${child})
      set( ENABLE_${child} ON CACHE BOOL "" FORCE )
      list( APPEND externalComponents ${child} )
    else()
      set( ENABLE_${child} OFF CACHE BOOL "" FORCE )
    endif()
  endif()
endforeach()
set(externalComponentsList "${externalComponents}" CACHE STRING "" FORCE)


set( externalComponentsLinkList "" )

if( ${BUILD_SHARED_LIBS} )
  set( externalComponentsLinkList ${externalComponentsList} )
else()
  foreach(component ${externalComponentsList})
    set( externalComponentsLinkList ${externalComponentsLinkList} ${GEOSX_LINK_PREPEND_FLAG} ${component} ${GEOSX_LINK_POSTPEND_FLAG} ) 
  endforeach()
endif()

set( externalComponentsLinkList ${externalComponentsLinkList} PARENT_SCOPE )



message( "Leaving /src/externalComponents/CMakeLists.txt" )
