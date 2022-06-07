
set( python_headers pygeosx.hpp
		    )

set( python_sources pygeosx.cpp
		    )
	    	  
set( dependencyList ${extraComponentsLinkList} ${externalComponentsLinkList} Python3::Python ) 

blt_add_library( NAME pygeosx
                 SOURCES ${python_sources}
                 HEADERS ${python_headers}
                 DEPENDS_ON ${dependencyList}
                 SHARED TRUE
                 CLEAR_PREFIX TRUE
               )

add_custom_command( TARGET pygeosx
                    PRE_BUILD
                    COMMAND ${Python3_EXECUTABLE} -m venv ${CMAKE_BINARY_DIR}/lib/PYGEOSX --system-site-packages )

add_custom_command( TARGET pygeosx
                    POST_BUILD
                    COMMAND ln -sf ${CMAKE_BINARY_DIR}/lib/pygeosx.so ${CMAKE_BINARY_DIR}/lib/PYGEOSX/lib/python3.8/site-packages/ )

add_custom_command( TARGET pygeosx
                    POST_BUILD
                    COMMAND ln -sf ${CMAKE_BINARY_DIR}/lib/pylvarray.so ${CMAKE_BINARY_DIR}/lib/PYGEOSX/lib/python3.8/site-packages/ )

target_include_directories( pygeosx PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

geosx_add_code_checks( PREFIX pygeosx )

install(TARGETS pygeosx LIBRARY DESTINATION lib)

macro(parse_test_file filename)
  file(READ "${filename}" contents)
  string(REPLACE "\n" ";" rows ${contents})
  foreach(line ${rows})
    string(REGEX MATCHALL "^ *def *test_([A-Za-z_0-9])+ *" found ${line})
    if (found)
      string(REGEX REPLACE "^ *def *test_(([A-Za-z_0-9])+) *" "\\1" name "${found}")
      add_test(NAME ${name}COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/${filename}::test_${name})
    endif()
  endforeach()
endmacro(parse_test_file)

#parse_test_file(test_WaveEquationAdjoint.py)
