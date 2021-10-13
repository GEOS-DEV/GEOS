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
 * @file TrilinosInterface.cpp
 */

#include "TrilinosInterface.hpp"

#include "linearAlgebra/interfaces/direct/SuiteSparse.hpp"
#include "linearAlgebra/interfaces/direct/SuperLUDist.hpp"
#include "linearAlgebra/interfaces/trilinos/TrilinosPreconditioner.hpp"
#include "linearAlgebra/interfaces/trilinos/TrilinosSolver.hpp"

namespace geosx
{

void TrilinosInterface::initialize()
{}

void TrilinosInterface::finalize()
{}

std::unique_ptr< LinearSolverBase< TrilinosInterface > >
TrilinosInterface::createSolver( LinearSolverParameters params )
{
  if( params.solverType == LinearSolverParameters::SolverType::direct )
  {
    if( params.direct.parallel )
    {
      return std::make_unique< SuperLUDist< TrilinosInterface > >( std::move( params ) );
    }
    else
    {
      return std::make_unique< SuiteSparse< TrilinosInterface > >( std::move( params ) );
    }
  }
  else
  {
    return std::make_unique< TrilinosSolver >( std::move( params ) );
  }
}

std::unique_ptr< PreconditionerBase< TrilinosInterface > >
TrilinosInterface::createPreconditioner( LinearSolverParameters params )
{
  return std::make_unique< TrilinosPreconditioner >( params );
}

std::unique_ptr< PreconditionerBase< TrilinosInterface > >
TrilinosInterface::createPreconditioner( LinearSolverParameters params,
                                         array1d< EpetraVector > const & nearNullKernel )
{
  return std::make_unique< TrilinosPreconditioner >( params, nearNullKernel );
}

}
