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

#include <_hypre_parcsr_mv.h>
#include <_hypre_parcsr_ls.h>

#include <numeric>

namespace geos
{

namespace hypre
{

HYPRE_Vector parVectorToVectorAll( HYPRE_ParVector const vec )
{
  if( hypre_ParVectorMemoryLocation( vec ) == HYPRE_MEMORY_HOST )
  {
    return (HYPRE_Vector)hypre_ParVectorToVectorAll( vec );
  }
  else
  {
    // Explicitly copy data to host before gathering, since hypre_ParVectorToVectorAll relies on UM.
    hypre_ParVector * const hostVec = hypre_ParVectorCloneDeep_v2( vec, HYPRE_MEMORY_HOST );
    hypre_Vector * const fullVec = hypre_ParVectorToVectorAll( hostVec );
    hypre_ParVectorDestroy( hostVec );
    return (HYPRE_Vector)fullVec;
  }
}

HYPRE_Vector parVectorToVector( HYPRE_ParVector const vec, int const targetRank )
{
  // If input vector on device, clone to host so MPI can access its data
  hypre_ParVector * const hostVec = hypre_ParVectorMemoryLocation( vec ) == HYPRE_MEMORY_HOST
                                  ? (hypre_ParVector *)vec
                                  : hypre_ParVectorCloneDeep_v2( vec, HYPRE_MEMORY_HOST );

  MPI_Comm const comm = hypre_ParVectorComm( hostVec );
  int const rank = MpiWrapper::commRank( comm );
  int const numProcs = MpiWrapper::commSize( comm );
  int constexpr tag = 1337;

  HYPRE_Int const globalSize = LvArray::integerConversion< HYPRE_Int >( hypre_ParVectorGlobalSize( hostVec ) );
  HYPRE_Int const firstIndex = LvArray::integerConversion< HYPRE_Int >( hypre_ParVectorFirstIndex( hostVec ) );

  hypre_Vector const * localVec = hypre_ParVectorLocalVector( hostVec );
  HYPRE_Int const localSize = hypre_VectorSize( localVec );
  HYPRE_Real const * const localData = hypre_VectorData( localVec );

  array1d< HYPRE_Int > offsets;
  if( rank == targetRank )
  {
    offsets.resize( numProcs + 1 );
    offsets[numProcs] = globalSize;
  }
  MpiWrapper::gather( &firstIndex, 1, offsets.data(), 1, targetRank, comm );

  hypre_Vector * newVec{};

  if( rank == targetRank )
  {
    newVec = hypre_SeqVectorCreate( globalSize );
    hypre_SeqVectorInitialize_v2( newVec, HYPRE_MEMORY_HOST );

    std::vector< MPI_Request > requests( numProcs, MPI_REQUEST_NULL );
    for( int i = 0; i < numProcs; ++i )
    {
      HYPRE_Real * const data = hypre_VectorData( newVec ) + offsets[i];
      if( i != rank )
      {
        MpiWrapper::iRecv( data, offsets[i + 1] - offsets[i], i, tag, comm, &requests[i] );
      }
      else
      {
        std::copy( localData, localData + localSize, data );
      }
    }
    MpiWrapper::waitAll( numProcs, requests.data(), MPI_STATUSES_IGNORE );
  }
  else
  {
    MpiWrapper::send( localData, localSize, targetRank, tag, comm );
  }

  if( hostVec != vec )
  {
    hypre_ParVectorDestroy( hostVec );
  }

  return (HYPRE_Vector)newVec;
}

HYPRE_Int dummySetup( HYPRE_Solver,
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
  GEOS_UNUSED_VAR( A );
#if defined(GEOSX_USE_SUPERLU_DIST)
  return hypre_SLUDistSolve( solver, b, x );
#else
  GEOS_UNUSED_VAR( solver );
  GEOS_UNUSED_VAR( b );
  GEOS_UNUSED_VAR( x );
  GEOS_ERROR( "GEOSX is configured without support for SuperLU_dist." );
  return -1;
#endif
}

HYPRE_Int SuperLUDistDestroy( HYPRE_Solver solver )
{
#if defined(GEOSX_USE_SUPERLU_DIST)
  return hypre_SLUDistDestroy( solver );
#else
  GEOS_UNUSED_VAR( solver );
  GEOS_ERROR( "GEOSX is configured without support for SuperLU_dist." );
  return -1;
#endif
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
  HYPRE_Real * norms{};
};

HYPRE_Int relaxationCreate( HYPRE_Solver & solver,
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

HYPRE_Int relaxationSetup( HYPRE_Solver solver,
                           HYPRE_ParCSRMatrix A,
                           HYPRE_ParVector b,
                           HYPRE_ParVector x )
{
  GEOS_UNUSED_VAR( b, x );

  // Refer to RelaxationData doxygen above for explanation of reinterpret_cast
  RelaxationData * const data = reinterpret_cast< RelaxationData * >( solver );
  data->vtemp.create( hypre_ParCSRMatrixNumRows( A ), hypre_ParCSRMatrixComm( A ) );
  data->ztemp.create( hypre_ParCSRMatrixNumRows( A ), hypre_ParCSRMatrixComm( A ) );

  // L1-smoothers need row L1-norms precomputed
  HYPRE_Int const normType = getL1NormType( data->type );
  if( normType > 0 )
  {
    if( data->norms )
    {
      hypre_TFree( data->norms, hypre::memoryLocation );
    }
    GEOS_LAI_CHECK_ERROR( hypre_ParCSRComputeL1Norms( A, normType, nullptr, &data->norms ) );
  }
  return 0;
}

HYPRE_Int relaxationSolve( HYPRE_Solver solver,
                           HYPRE_ParCSRMatrix A,
                           HYPRE_ParVector b,
                           HYPRE_ParVector x )
{
  // Refer to RelaxationData doxygen above for explanation of reinterpret_cast
  RelaxationData * const data = reinterpret_cast< RelaxationData * >( solver );
  return hypre_BoomerAMGRelax( A, b, nullptr, data->type, 0, 1.0, 1.0, data->norms, x, data->vtemp.unwrapped(), data->ztemp.unwrapped() );
}

HYPRE_Int relaxationDestroy( HYPRE_Solver solver )
{
  // Refer to RelaxationData doxygen above for explanation of reinterpret_cast
  RelaxationData * const data = reinterpret_cast< RelaxationData * >( solver );
  GEOS_ASSERT( data != nullptr );
  if( data->norms )
  {
    hypre_TFree( data->norms, hypre::memoryLocation );
  }
  delete data;
  return 0;
}

/**
 * @brief Holds data needed by hypre's Chebyshev smoothing.
 *
 * See RelaxationData for explanation of the way this structure is used.
 */
struct ChebyshevData
{
  HYPRE_Int order = 1;       ///< Chebyshev order
  HYPRE_Int eigNumIter = 10; ///< Number of eigenvalue estimation iterations

  HYPRE_Real * diag{};   ///< Scaled diagonal values extracted during setup
  HYPRE_Real * coefs{};  ///< Precomputed coefficients

  // Temporary vectors required by Chebyshev solve
  HypreVector vtemp;
  HypreVector ztemp;
  HypreVector ptemp;
  HypreVector rtemp;
};

HYPRE_Int chebyshevCreate( HYPRE_Solver & solver,
                           HYPRE_Int const order,
                           HYPRE_Int const eigNumIter )
{
  ChebyshevData * const data = new ChebyshevData;
  data->order = order;
  data->eigNumIter = eigNumIter;
  solver = reinterpret_cast< HYPRE_Solver >( data );
  return 0;
}

HYPRE_Int chebyshevSetup( HYPRE_Solver solver,
                          HYPRE_ParCSRMatrix A,
                          HYPRE_ParVector b,
                          HYPRE_ParVector x )
{
  GEOS_UNUSED_VAR( b, x );

  // Refer to ChebyshevData doxygen above for explanation of reinterpret_cast
  ChebyshevData * const data = reinterpret_cast< ChebyshevData * >( solver );
  data->vtemp.create( hypre_ParCSRMatrixNumRows( A ), hypre_ParCSRMatrixComm( A ) );
  data->ztemp.create( hypre_ParCSRMatrixNumRows( A ), hypre_ParCSRMatrixComm( A ) );
  data->ptemp.create( hypre_ParCSRMatrixNumRows( A ), hypre_ParCSRMatrixComm( A ) );
  data->rtemp.create( hypre_ParCSRMatrixNumRows( A ), hypre_ParCSRMatrixComm( A ) );

  HYPRE_Real max_eig, min_eig;
  if( data->eigNumIter > 0 )
  {
    GEOS_LAI_CHECK_ERROR( hypre_ParCSRMaxEigEstimateCG( A, 1, data->eigNumIter, &max_eig, &min_eig ) );
  }
  else
  {
    GEOS_LAI_CHECK_ERROR( hypre_ParCSRMaxEigEstimate( A, 1, &max_eig, &min_eig ) );
  }

  return hypre_ParCSRRelax_Cheby_Setup( A,
                                        max_eig,
                                        min_eig,
                                        0.3, // fraction
                                        data->order,
                                        1, // scale
                                        0, // variant
                                        &data->coefs,
                                        &data->diag );
}

HYPRE_Int chebyshevSolve( HYPRE_Solver solver,
                          HYPRE_ParCSRMatrix A,
                          HYPRE_ParVector b,
                          HYPRE_ParVector x )
{
  // Refer to ChebyshevData doxygen above for explanation of reinterpret_cast
  ChebyshevData * const data = reinterpret_cast< ChebyshevData * >( solver );

  return hypre_ParCSRRelax_Cheby_Solve( A,
                                        b,
                                        data->diag,
                                        data->coefs,
                                        data->order,
                                        1, // scale
                                        0, // variant
                                        x,
                                        data->vtemp.unwrapped(),
                                        data->ztemp.unwrapped(),
                                        data->ptemp.unwrapped(),
                                        data->rtemp.unwrapped() );
}

HYPRE_Int chebyshevDestroy( HYPRE_Solver solver )
{
  // Refer to ChebyshevData doxygen above for explanation of reinterpret_cast
  ChebyshevData * const data = reinterpret_cast< ChebyshevData * >( solver );
  GEOS_ASSERT( data != nullptr );
  if( data->diag )
  {
    hypre_TFree( data->diag, hypre::memoryLocation );
  }
  if( data->coefs )
  {
    hypre_TFree( data->coefs, HYPRE_MEMORY_HOST );
  }
  delete data;
  return 0;
}

} // namespace hypre

} // namespace geos
