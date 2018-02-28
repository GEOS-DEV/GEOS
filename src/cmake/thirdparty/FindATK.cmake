###############################################################################
# Copyright (c) 2014-2015, Lawrence Livermore National Security, LLC.
# 
# Produced at the Lawrence Livermore National Laboratory
# 
# LLNL-CODE-666778
# 
# All rights reserved.
# 
# This file is part of Conduit. 
# 
# For details, see https://lc.llnl.gov/conduit/.
# 
# Please also read conduit/LICENSE
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# 
# * Redistributions of source code must retain the above copyright notice, 
#   this list of conditions and the disclaimer below.
# 
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the disclaimer (as noted below) in the
#   documentation and/or other materials provided with the distribution.
# 
# * Neither the name of the LLNS/LLNL nor the names of its contributors may
#   be used to endorse or promote products derived from this software without
#   specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
# LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
# DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
# POSSIBILITY OF SUCH DAMAGE.
# 
###############################################################################
#
# Setup ATK
# This file defines:
#  ATK_FOUND - If ATK was found
#  ATK_INCLUDE_DIRS - The ATK include directories
#  
#  If found, the conduit CMake targets will also be imported

# first Check for ATK_DIR

if(ATK_DIR)
    if(NOT EXISTS ${ATK_DIR}/lib/cmake/axom_utils-targets.cmake)
        MESSAGE(FATAL_ERROR "Could not find ATK cmake include file (${ATK_DIR}/lib/cmake/axom_utils-targets.cmake)")
    endif()
    include(${ATK_DIR}/lib/cmake/axom_utils-targets.cmake)


    #if(ENABLE_MPI)
      if(NOT EXISTS ${ATK_DIR}/lib/cmake/lumberjack-targets.cmake)
          MESSAGE(FATAL_ERROR "Could not find ATK cmake include file (${ATK_DIR}/lib/cmake/lumberjack-targets.cmake)")
      endif()
    include(${ATK_DIR}/lib/cmake/lumberjack-targets.cmake)
    #endif(ENABLE_MPI)

    if(NOT EXISTS ${ATK_DIR}/lib/cmake/slic-targets.cmake)
        MESSAGE(FATAL_ERROR "Could not find ATK cmake include file (${ATK_DIR}/lib/cmake/slic-targets.cmake)")
    endif()
    include(${ATK_DIR}/lib/cmake/slic-targets.cmake)


    if(NOT EXISTS ${ATK_DIR}/lib/cmake/sidre-targets.cmake)
        MESSAGE(FATAL_ERROR "Could not find ATK cmake include file (${ATK_DIR}/lib/cmake/sidre-targets.cmake)")
    endif()
    include(${ATK_DIR}/lib/cmake/sidre-targets.cmake)


    #include("${ATK_CMAKE}/lumberjack-targets.cmake")
    #include("${ATK_CMAKE}/sidre-targets.cmake")
    #include("${ATK_CMAKE}/slic-targets.cmake")
    include("${ATK_CMAKE}/mint-targets.cmake")
    include("${ATK_CMAKE}/fmt-targets.cmake")
    include("${ATK_CMAKE}/primal-targets.cmake")
    include("${ATK_CMAKE}/slam-targets.cmake")
    include("${ATK_CMAKE}/quest-targets.cmake")
    include("${ATK_CMAKE}/slam-targets.cmake")
    include("${ATK_CMAKE}/spio-targets.cmake")

    set(ATK_FOUND TRUE)    
    
    set(ATK_INCLUDE_DIRS ${ATK_DIR}/include)
else()
    set(ATK_FOUND FALSE)
endif()
