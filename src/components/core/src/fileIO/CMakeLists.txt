
#
# Specify all headers
#
set(fileIO_headers
  xmlWrapper.hpp
  silo/SiloFile.hpp
  blueprint/Blueprint.hpp
  utils/utils.hpp )

#
# Specify all sources
#
set(fileIO_sources
  xmlWrapper.cpp
  silo/SiloFile.cpp
  blueprint/Blueprint.cpp )

    
blt_add_library( NAME                  fileIO
                 SOURCES               ${fileIO_sources}
                 HEADERS               ${fileIO_headers}
                 DEPENDS_ON            common silo hdf5 mpi pugixml RAJA
                 SHARED FALSE
               )