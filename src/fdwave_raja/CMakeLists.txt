###############################################################################
# Copyright (c) 2016-21, Lawrence Livermore National Security, LLC
# and RAJA project contributors. See the RAJA/COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
###############################################################################
set (fdwave_sources
     grid.cpp
     pml.cpp
     data_setup.cpp
     target_3d.cpp)

set( fdwave_headers 
     memoryManager.hpp
     constants.hpp
     grid.hpp
     pml.hpp
     data_setup.hpp
     target_3d.hpp)

set( dependencyList lvarray  RAJA chai common)

blt_add_library( NAME          fdwavelib_raja
                 SOURCES       ${fdwave_sources}
                 HEADERS       ${fdwave_headers}
                 DEPENDS_ON    ${dependencyList}
                 OBJECT        ${buildAsObj}
	       )

set( extraComponentsLinkList ${extraComponentsLinkList} fdwavelib_raja)
blt_add_executable(NAME fdwave-eqn_raja
                   SOURCES fdwave-eqn.cpp
                   DEPENDS_ON ${extraComponentsLinkList}
                              ${externalComponentsLinkList}
             )
