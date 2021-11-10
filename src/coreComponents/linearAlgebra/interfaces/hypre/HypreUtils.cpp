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
 * @file HypreUtils.cpp
 */

#include "HypreUtils.hpp"

#include "linearAlgebra/interfaces/hypre/HypreVector.hpp"

#include <_hypre_parcsr_ls.h>

namespace geosx
{

namespace hypre
{

HYPRE_Int DummySetup( HYPRE_Solver,
                      HYPRE_ParCSRMatrix,
                      HYPRE_ParVector,
                      HYPRE_ParVector )
{
  return 0;
}

HYPRE_Int SuperLUDistSolve( HYPRE_Solver solver,
                            HYPRE_ParCSRMatrix A,
                            HYPRE_ParVector b,
                            HYPRE_ParVector x )
{
  GEOSX_UNUSED_VAR( A );
  return hypre_SLUDistSolve( solver, b, x );
}

HYPRE_Int SuperLUDistDestroy( HYPRE_Solver solver )
{
  return hypre_SLUDistDestroy( solver );
}

/**
 * @brief Holds temporary data needed by hypre's relaxation methods.
 *
 * Many implementations dispatched to by hypre_BoomerAMGRelax crash unless Vtemp,
 * Ztemp (temporary vectors) or l1Norms (row norm vector) parameters are provided.
 * We want these temporary vectors to persist across multiple applications of the relaxation method.
 * In order to make it conform to hypre's solver interface, we allocate these vectors during setup
 * and disguise them behind HYPRE_Solver, hypre's opaque solver struct pointer type.
 * Hence multiple reinterpret_casts (replacing C-style casts used in hypre) are required
 * to convert back and forth between HYPRE_Solver and pointer to actual data struct.
 * We also have to manage memory with raw new/delete here, because a single pointer
 * that is maintained by hypre is all we got to manage the struct with.
 */
struct RelaxationData
{
  HYPRE_Int type = -1;
  HypreVector vtemp;
  HypreVector ztemp;
  HypreVector norms;
};

HYPRE_Int RelaxationCreate( HYPRE_Solver & solver,
                            HYPRE_Int const type )
{
  RelaxationData * const data = new RelaxationData;
  data->type = type;
  solver = reinterpret_cast< HYPRE_Solver >( data );
  return 0;
}

HYPRE_Int getL1NormType( HYPRE_Int const relaxType )
{
  // Extracted from the depths of hypre_BoomerAMGSetup()
  switch( relaxType )
  {
    case 7: return 5;
    case 8:
    case 13:
    case 14: return 4;
    case 18: return 1;
    default: return -1;
  }
}

HYPRE_Int RelaxationSetup( HYPRE_Solver solver,
                           HYPRE_ParCSRMatrix A,
                           HYPRE_ParVector b,
                           HYPRE_ParVector x )
{
  GEOSX_UNUSED_VAR( b, x );

  // Refer to RelaxationData doxygen above for explanation of reinterpret_cast
  RelaxationData * const data = reinterpret_cast< RelaxationData * >( solver );
  data->vtemp.createWithLocalSize( hypre_ParCSRMatrixNumRows( A ), hypre_ParCSRMatrixComm( A ) );
  data->ztemp.createWithLocalSize( hypre_ParCSRMatrixNumRows( A ), hypre_ParCSRMatrixComm( A ) );

  // L1-smoothers need row L1-norms precomputed
  HYPRE_Int const normType = getL1NormType( data->type );
  if( normType > 0 )
  {
    data->norms.createWithLocalSize( hypre_ParCSRMatrixNumRows( A ), hypre_ParCSRMatrixComm( A ) );
    HYPRE_Real * l1Norms = nullptr;
    hypre_ParCSRComputeL1Norms( A, normType, nullptr, &l1Norms );

    // Sadly, hypre's L1 norms method always allocates new storage for result.
    // Therefore we "sneak" it into the vector while deleting old storage.
    hypre_Vector * const localVector = hypre_ParVectorLocalVector( data->norms.unwrapped() );
    hypre_TFree( hypre_VectorData( localVector ), hypre_VectorMemoryLocation( localVector ) );
    hypre_VectorData( localVector ) = l1Norms;
  }
  return 0;
}

HYPRE_Int RelaxationSolve( HYPRE_Solver solver,
                           HYPRE_ParCSRMatrix A,
                           HYPRE_ParVector b,
                           HYPRE_ParVector x )
{
  // Refer to RelaxationData doxygen above for explanation of reinterpret_cast
  RelaxationData * const data = reinterpret_cast< RelaxationData * >( solver );
  HYPRE_Real * const l1Norms = data->norms.ready() ? data->norms.extractLocalVector() : nullptr;
  return hypre_BoomerAMGRelax( A, b, nullptr, data->type, 0, 1.0, 1.0, l1Norms, x, data->vtemp.unwrapped(), data->ztemp.unwrapped() );
}

HYPRE_Int RelaxationDestroy( HYPRE_Solver solver )
{
  // Refer to RelaxationData doxygen above for explanation of reinterpret_cast
  RelaxationData * const data = reinterpret_cast< RelaxationData * >( solver );
  delete data;
  return 0;
}

} // namespace hypre

} // namespace geosx
