###############################################################################
# Copyright (c) 2014, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-666778
#
# All rights reserved.
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

################################
# Setup build options and their default values
################################
include(core/src/cmake/geosxOptions.cmake)

################################
# Setup toolkit generate targets
################################
#include(cmake/SetupShroud.cmake)

################################
# ATK's Third party library setup
################################
include(core/src/cmake/thirdparty/SetupGeosxThirdParty.cmake)

#
# We don't try to use this approach for CMake generators that support
# multiple configurations. See: CZ JIRA: ATK-45
#
if(NOT CMAKE_CONFIGURATION_TYPES)
    ######################################################
    # Add define we can use when debug builds are enabled
    ######################################################
    if( (CMAKE_BUILD_TYPE MATCHES Debug)
        OR (CMAKE_BUILD_TYPE MATCHES RelWithDebInfo )
      )
        add_definitions(-DATK_DEBUG)
    endif()
endif()

################################
# Fortran Configuration
################################
if(ENABLE_FORTRAN)
    # Create macros for Fortran name mangling
    include(FortranCInterface)
    FortranCInterface_HEADER(${HEADER_INCLUDES_DIRECTORY}/common/FC.h MACRO_NAMESPACE "FC_")

    if (ENABLE_MPI)
        # Determine if we should use fortran mpif.h header or fortran mpi module
        find_path(mpif_path
            NAMES "mpif.h"
            PATHS ${MPI_Fortran_INCLUDE_PATH}
            NO_DEFAULT_PATH
            )
        
        if(mpif_path)
            set(MPI_Fortran_USE_MPIF ON CACHE PATH "")
            message(STATUS "Using MPI Fortran header: mpif.h")
        else()
            set(MPI_Fortran_USE_MPIF OFF CACHE PATH "")
            message(STATUS "Using MPI Fortran module: mpi.mod")
        endif()
    endif()
endif()


##############################################################################
# Setup some additional compiler options that can be useful in various targets
# These are stored in their own variables.
# Usage: To add one of these sets of flags to some source files:
#   get_source_file_property(_origflags <src_file> COMPILE_FLAGS)
#   set_source_files_properties(<list_of_src_files> 
#        PROPERTIES COMPILE_FLAGS "${_origFlags} ${<flags_variable}" )
##############################################################################

set(custom_compiler_flags_list) # Tracks custom compiler flags for logging

# Flag for disabling warnings about omp pragmas in the code
blt_append_custom_compiler_flag(FLAGS_VAR ATK_DISABLE_OMP_PRAGMA_WARNINGS
                  DEFAULT "-Wno-unknown-pragmas"
                  XL      "-qignprag=omp"
                  INTEL   "-diag-disable 3180"
                  )
list(APPEND custom_compiler_flags_list ATK_DISABLE_OMP_PRAGMA_WARNINGS)

# Flag for disabling warnings about unused parameters.
# Useful when we include external code.
blt_append_custom_compiler_flag(FLAGS_VAR ATK_DISABLE_UNUSED_PARAMETER_WARNINGS
                  DEFAULT "-Wno-unused-parameter"
                  XL      "-qinfo=nopar"
                  )
list(APPEND custom_compiler_flags_list ATK_DISABLE_UNUSED_PARAMETER_WARNINGS)

# Flag for disabling warnings about unused variables
# Useful when we include external code.
blt_append_custom_compiler_flag(FLAGS_VAR ATK_DISABLE_UNUSED_VARIABLE_WARNINGS
                  DEFAULT "-Wno-unused-variable"
                  XL      "-qinfo=nouse"
                  )
list(APPEND custom_compiler_flags_list ATK_DISABLE_UNUSED_VARIABLE_WARNINGS)

# Flag for disabling warnings about variables that may be uninitialized.
# Useful when we are using compiler generated interface code (e.g. in shroud)
blt_append_custom_compiler_flag(FLAGS_VAR ATK_DISABLE_UNINITIALIZED_WARNINGS
                  DEFAULT "-Wno-uninitialized"
                  XL      "-qsuppress=1540-1102"
                  )
list(APPEND custom_compiler_flags_list ATK_DISABLE_UNINITIALIZED_WARNINGS)

# Flag for disabling warnings about strict aliasing.
# Useful when we are using compiler generated interface code (e.g. in shroud)
blt_append_custom_compiler_flag(FLAGS_VAR ATK_DISABLE_ALIASING_WARNINGS
                  DEFAULT "-Wno-strict-aliasing"
                  XL      ""
                  )
list(APPEND custom_compiler_flags_list ATK_DISABLE_ALIASING_WARNINGS)

# Flag for enabling the C preprocessor in fortran.
# (Note KW 5/2016) The XL flag only applies to *.f90 files -- I could not find a more general solution.   
#     xlf only allows one file remapping at a time. If you have *.f files, '-qsuffix=cpp=f' should work.
#     Alternatively, you can rename the file's extension to automatically invoke the preprocessor (e.g. *.f90 ->  *.F90)
blt_append_custom_compiler_flag(FLAGS_VAR ATK_PREPROCESS_FORTRAN
                  DEFAULT "-cpp"
                  XL      "-qsuffix=cpp=f90"  # Note: only invokes the preprocessor on files with extension *.f90
                  )
list(APPEND custom_compiler_flags_list ATK_PREPROCESS_FORTRAN)

   
# message(STATUS "Custom compiler flags:")
# foreach(flag ${custom_compiler_flags_list})
#    message(STATUS "\tvalue of ${flag} is '${${flag}}'")
# endforeach()



macro(subdirlist result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
        LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()