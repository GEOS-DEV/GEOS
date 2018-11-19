message( "Entering /src/externalComponents/CMakeLists.txt")



set(externalComponentsList "" CACHE STRING "" FORCE)

file(GLOB children RELATIVE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/*)
foreach(child ${children})
  if(IS_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/${child})
    if( EXISTS "${CMAKE_CURRENT_LIST_DIR}/${child}/CMakeLists.txt")
      add_subdirectory(${child})
      set( ENABLE_${child} ON CACHE BOOL "" FORCE )
      list( APPEND externalComponentsList ${child} )
    else()
      set( ENABLE_${child} OFF CACHE BOOL "" FORCE )
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
