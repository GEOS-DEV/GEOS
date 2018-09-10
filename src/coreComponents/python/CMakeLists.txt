set(ENV{LD_LIBRARY_PATH}, ENV{LD_LIBRARY_PATH}:${PYTHON_DIR}/lib/ )

set(PYTHON_INTERP_INCLUDE_DIRS
    ${PYTHON_DIR}/include/python2.7/
    ${PYTHON_DIR}/lib/python2.7/site-packages/numpy/core/include/
    )

blt_register_library( NAME python_interp
                      INCLUDES    ${PYTHON_INTERP_INCLUDE_DIRS}
                      HEADERS     ${PYTHON_INTERP_INCLUDE_DIRS}
                      HEADERS_OUTPUT_SUBDIR python_interp
                      LIBRARIES python2.7 
                      )
