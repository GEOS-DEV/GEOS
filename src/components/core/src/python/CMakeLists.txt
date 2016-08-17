set(anaconda
     /g/g17/sherman/anaconda2)

set(python_interp_headers
    ${anaconda}/include/python2.7/Python.h
    ${anaconda}/lib/python2.7/site-packages/numpy/core/include/numpy/arrayobject.h )


message("Linking the python interpreter")
blt_add_library( NAME             python_interp
                 HEADERS          ${python_interp_headers}
                 HEADERS_OUTPUT_SUBDIR python_interp
                 DEPENDS_ON        
                 SYSTEM_DEPENDS_ON
                 )

set_target_properties(python_interp PROPERTIES LINKER_LANGUAGE C)
# LINK_FLAGS "-lpython2.7"
# INCLUDES         ${anaconda}/lib/libpython2.7.so