####################################
# 3rd Party Dependencies
####################################

################################
# Conduit
################################
#if (CONDUIT_DIR)
#  include(cmake/thirdparty/FindConduit.cmake)
#  blt_register_library( NAME conduit
#                        INCLUDES ${CONDUIT_INCLUDE_DIRS} 
#                        LIBRARIES  conduit)
#  blt_register_library( NAME conduit_io
#                        INCLUDES ${CONDUIT_INCLUDE_DIRS}
#                        LIBRARIES  conduit_io)
#endif()

################################
# ATK
################################
#if (ATK_DIR)
#  include(cmake/thirdparty/FindATK.cmake)
#  blt_register_library( NAME atk
#                        INCLUDES ${ATK_INCLUDE_DIRS} 
#                        LIBRARIES  sidre)
#endif()

if (ATK_DIR)
  include(cmake/thirdparty/FindATK.cmake)
  blt_register_library( NAME sidre
                        INCLUDES ${ATK_INCLUDE_DIRS} 
                        LIBRARIES  sidre)

  blt_register_library( NAME slic
                        INCLUDES ${ATK_INCLUDE_DIRS} 
                        LIBRARIES  slic)
endif()
