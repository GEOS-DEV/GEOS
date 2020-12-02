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
 * @file HypreInterface.cpp
 */

#include "HypreInterface.hpp"
#include "linearAlgebra/interfaces/hypre/HyprePreconditioner.hpp"

namespace geosx
{

void HypreInterface::initialize( int & GEOSX_UNUSED_PARAM( argc ),
                                 char * * & GEOSX_UNUSED_PARAM( argv ) )
{}

void HypreInterface::finalize()
{}

std::unique_ptr< PreconditionerBase< HypreInterface > >
geosx::HypreInterface::createPreconditioner( LinearSolverParameters params )
{
  return std::make_unique< HyprePreconditioner >( params );
}

std::unique_ptr< PreconditionerBase< HypreInterface > >
geosx::HypreInterface::createPreconditioner( LinearSolverParameters params,
                                             array1d< HypreVector > const & nearNullKernel )
{
  return std::make_unique< HyprePreconditioner >( params, nearNullKernel );
}

}
