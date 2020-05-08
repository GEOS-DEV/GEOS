/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TrilinosTpetraInterface.cpp
 */

#include "TrilinosTpetraInterface.hpp"

#include "linearAlgebra/interfaces/trilinos/TrilinosTpetraPreconditioner.hpp"

#include <Tpetra_Core.hpp>

namespace geosx
{

void TrilinosTpetraInterface::initialize( int & argc,
                                          char * * & argv )
{
  Tpetra::initialize( &argc, &argv );
}

void TrilinosTpetraInterface::finalize()
{
  Tpetra::finalize();
}

std::unique_ptr< PreconditionerBase< TrilinosTpetraInterface > >
TrilinosTpetraInterface::createPreconditioner( LinearSolverParameters params )
{
  return std::make_unique< TrilinosTpetraPreconditioner >( params );
}

}
