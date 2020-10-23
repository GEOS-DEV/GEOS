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
#include "HYPRE_utilities.h"
#include "_hypre_utilities.h"
#include "_hypre_utilities.hpp"

#include "HypreMatrix.hpp"

namespace geosx
{

void HypreInterface::initialize( int & GEOSX_UNUSED_PARAM( argc ),
                                 char * * & GEOSX_UNUSED_PARAM( argv ) )
{
  HYPRE_Init();
#if defined(OVERRIDE_CREATE)
#if defined(GEOSX_USE_CUDA)
  hypre_HandleDefaultExecPolicy(hypre_handle()) = HYPRE_EXEC_DEVICE;
  hypre_HandleSpgemmUseCusparse(hypre_handle()) = 0;
#endif
#endif
}

void HypreInterface::finalize()
{
  HYPRE_Finalize();
}

std::unique_ptr< PreconditionerBase< HypreInterface > >
geosx::HypreInterface::createPreconditioner( LinearSolverParameters params )
{
  return std::make_unique< HyprePreconditioner >( params );
}

}
