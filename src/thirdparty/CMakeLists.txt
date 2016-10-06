
#blt_add_component(COMPONENT_NAME optionparser DEFAULT_STATE ON)


# find RAJA, or build if RAJA_DIR is not defined
#if(RAJA_DIR)
#   include(../cmake/thirdparty/FindRAJA.cmake)
#   if (NOT RAJA_FOUND)
#	MESSAGE(FATAL_ERROR "RAJA not found. Do not define RAJA_DIR to attempt internal build of RAJA.")
#   endif()
#elseif(EXISTS ${CMAKE_CURRENT_LIST_DIR}/raja/CMakeLists.txt)
#  MESSAGE("including raja CMakeLists.txt") 
#  set(RAJA_DIR ${CMAKE_CURRENT_LIST_DIR}/raja-install/share/raja/cmake)
#  find_package(RAJA REQUIRED)
#  include(../cmake/thirdparty/FindRAJA.cmake) 
#else()
#   MESSAGE(FATAL_ERROR "RAJA_DIR not defined and no local raja found in ${CMAKE_CURRENT_LIST_DIR}/raja/CMakeLists.txt.")
#endif(RAJA_DIR)

add_subdirectory(pugixml)
