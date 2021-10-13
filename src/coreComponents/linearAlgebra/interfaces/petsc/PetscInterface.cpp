/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PetscInterface.cpp
 */

#include "PetscInterface.hpp"

#include "linearAlgebra/interfaces/direct/SuiteSparse.hpp"
#include "linearAlgebra/interfaces/direct/SuperLUDist.hpp"
#include "linearAlgebra/interfaces/petsc/PetscPreconditioner.hpp"
#include "linearAlgebra/interfaces/petsc/PetscSolver.hpp"

#include <petscsys.h>

namespace geosx
{

void PetscInterface::initialize()
{
  PetscOptionsSetValue( nullptr, "-no_signal_handler", "" );
  PetscOptionsSetValue( nullptr, "-on_error_abort", "" );
  PETSC_COMM_WORLD = MPI_COMM_GEOSX;
  PetscInitializeNoArguments();
}

void PetscInterface::finalize()
{
  PetscFinalize();
}

std::unique_ptr< LinearSolverBase< PetscInterface > >
PetscInterface::createSolver( LinearSolverParameters params )
{
  if( params.solverType == LinearSolverParameters::SolverType::direct )
  {
    if( params.direct.parallel )
    {
      return std::make_unique< SuperLUDist< PetscInterface > >( std::move( params ) );
    }
    else
    {
      return std::make_unique< SuiteSparse< PetscInterface > >( std::move( params ) );
    }
  }
  else
  {
    return std::make_unique< PetscSolver >( std::move( params ) );
  }
}

std::unique_ptr< PreconditionerBase< PetscInterface > >
PetscInterface::createPreconditioner( LinearSolverParameters params )
{
  return std::make_unique< PetscPreconditioner >( params );
}

std::unique_ptr< PreconditionerBase< PetscInterface > >
PetscInterface::createPreconditioner( LinearSolverParameters params,
                                      array1d< PetscVector > const & nearNullKernel )
{
  return std::make_unique< PetscPreconditioner >( params, nearNullKernel );
}

} //namespace geosx
