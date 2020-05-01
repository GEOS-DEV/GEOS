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

if (NOT EXISTS ${ATK_DIR})
    message(STATUS "Using axom from thirdPartyLibs")
    set(ATK_DIR "${GEOSX_TPL_DIR}/axom" CACHE PATH "")
endif()

if (EXISTS ${ATK_DIR})
    message(STATUS "Using axom found at ${ATK_DIR}")

    set(ATK_FOUND TRUE)
    set(ATK_INCLUDE_DIRS ${ATK_DIR}/include)
    set(ATK_CMAKE ${ATK_DIR}/lib/cmake)

    include(${ATK_CMAKE}/fmt-targets.cmake)
    include(${ATK_CMAKE}/sparsehash-targets.cmake)
    include(${ATK_CMAKE}/axom-targets.cmake)

    # wrap imported compile flags properly for nvcc
    if(ENABLE_CUDA)
        get_property(_axom_compile_flags TARGET axom PROPERTY INTERFACE_COMPILE_OPTIONS SET)
        if(_axom_compile_flags)
            get_property(_axom_compile_flags TARGET axom PROPERTY INTERFACE_COMPILE_OPTIONS)
            set(_axom_compile_flags_wrapped)
            foreach(_flag ${_axom_compile_flags})
                list(APPEND _axom_compile_flags_wrapped
                     $<$<NOT:$<COMPILE_LANGUAGE:CUDA>>:${_flag}>
                     $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${_flag}>)
            endforeach()
            set_target_properties(axom PROPERTIES INTERFACE_COMPILE_OPTIONS "${_axom_compile_flags_wrapped}")
        endif()
    endif()

    blt_register_library (NAME axom
                          INCLUDES ${ATK_INCLUDE_DIRS}
                          LIBRARIES axom
                          TREAT_INCLUDES_AS_SYSTEM ON)
                          
    set(thirdPartyLibs ${thirdPartyLibs} axom )

else()
    set(ATK_FOUND FALSE)
    message("Not using axom")
endif()
