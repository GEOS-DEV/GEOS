message( "Entering src/externalComponents/CMakeLists.txt")



set(externalComponents "")

file(GLOB children RELATIVE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/*)
foreach(child ${children})
  if(IS_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/${child})
    message("Found directory at src/externalComponents/${child}")
    if( ${ENABLE_${child}} )
      if( (EXISTS "${CMAKE_CURRENT_LIST_DIR}/${child}/CMakeLists.txt")  )
        message("Found CMakeLists.txt at src/externalComponents/${child}....adding component!")
        add_subdirectory(${child})
        list( APPEND externalComponents ${child} )
      else()
        message( FATAL_ERROR "ENABLE_${child} is ON, but cannot find CMakeLists.txt" )
      endif()
    endif()
  endif()
endforeach()
set(externalComponentsList "${externalComponents}" CACHE STRING "" FORCE)

set( externalComponentsLinkList "" )

set( externalComponentsLinkList ${externalComponentsList} )

set( externalComponentsLinkList ${externalComponentsLinkList} PARENT_SCOPE )



message( "Leaving /src/externalComponents/CMakeLists.txt" )
