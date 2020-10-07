/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PetscSuiteSparse.cpp
 */

#include "PetscSuiteSparse.hpp"
#include <petsc.h>
#include <superlu_ddefs.h>

namespace geosx
{

void SuiteSparseSetFromOptions( PetscMatrix const & matrix,
                                LinearSolverParameters const & params )
{
  GEOSX_UNUSED_VAR( matrix );

  // Set options.
  if( params.logLevel > 0 )
  {
    PetscOptionsSetValue( nullptr, "-mat_umfpack_prl", "6" );
  }
  else
  {
    PetscOptionsSetValue( nullptr, "-mat_umfpack_prl", "1" );
  }
  PetscOptionsSetValue( nullptr, "-mat_umfpack_ordering", "BEST" );
}

}

