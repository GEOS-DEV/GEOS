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
 * @file HypreUtils.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREUTILS_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREUTILS_HPP_

#include "common/DataTypes.hpp"

#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>
#include <_hypre_parcsr_ls.h>

namespace geosx
{

/**
 * @brief Container for hypre preconditioner function pointers.
 *
 * @note: This needs to be here rather than in HyprePreconditioner.cpp,
 *        because HypreSolver needs to access `apply` member.
 */
struct HyprePrecWrapper
{
  /// Alias for setup function type
  using SetupFunc = HYPRE_Int (*)( HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector );

  /// Alias for apply function type
  using SolveFunc = HYPRE_Int (*)( HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector );

  /// Alias for destroy function type
  using DestroyFunc = HYPRE_Int (*)( HYPRE_Solver );

  HYPRE_Solver ptr{};    ///< pointer to preconditioner
  SetupFunc setup{};     ///< pointer to setup function
  SolveFunc solve{};     ///< pointer to apply function
  DestroyFunc destroy{}; ///< pointer to destroy function
};

/**
 * @brief Contains some hypre-specific functions
 */
namespace hypre
{

// Check matching requirements on index/value types between GEOSX and Hypre

// WARNING. We don't have consistent types between HYPRE_Int and localIndex.
//          Decision needs to be made either to use bigint option, or change
//          localIndex to int. We are getting away with this because we do not
//          pass ( localIndex * ) to hypre except when it is on the GPU, in
//          which case we are using int for localIndex.
#if defined(GEOSX_USE_HYPRE_CUDA)
static_assert( sizeof( HYPRE_Int ) == sizeof( geosx::localIndex ),
               "HYPRE_Int and geosx::localIndex must have the same size" );
static_assert( std::is_signed< HYPRE_Int >::value == std::is_signed< geosx::localIndex >::value,
               "HYPRE_Int and geosx::localIndex must both be signed or unsigned" );
#endif

static_assert( sizeof( HYPRE_BigInt ) == sizeof( geosx::globalIndex ),
               "HYPRE_BigInt and geosx::globalIndex must have the same size" );

static_assert( std::is_signed< HYPRE_BigInt >::value == std::is_signed< geosx::globalIndex >::value,
               "HYPRE_BigInt and geosx::globalIndex must both be signed or unsigned" );

static_assert( std::is_same< HYPRE_Real, geosx::real64 >::value,
               "HYPRE_Real and geosx::real64 must be the same type" );

#ifndef HYPRE_NO_GLOBAL_PARTITION
static_assert( false, "Hypre must be built with HYPRE_NO_GLOBAL_PARTITION" )
#endif

/**
 * @brief Converts a non-const array from GEOSX globalIndex type to HYPRE_BigInt
 * @param[in] index the input array
 * @return the converted array
 */
inline HYPRE_BigInt * toHypreBigInt( geosx::globalIndex * const index )
{
  return reinterpret_cast< HYPRE_BigInt * >(index);
}

/**
 * @brief Converts a const array from GEOSX globalIndex type to HYPRE_BigInt
 * @param[in] index the input array
 * @return the converted array
 */
inline HYPRE_BigInt const * toHypreBigInt( geosx::globalIndex const * const index )
{
  return reinterpret_cast< HYPRE_BigInt const * >(index);
}

/**
 * @brief Dummy function that does nothing but conform to hypre's signature for preconditioner setup/apply functions.
 * @return always 0 (success).
 *
 * Typical use is to prevent hypre from calling preconditioner setup when we have already called it on out side.
 */
inline HYPRE_Int HYPRE_DummySetup( HYPRE_Solver,
                                   HYPRE_ParCSRMatrix,
                                   HYPRE_ParVector,
                                   HYPRE_ParVector )
{
  return 0;
}

/**
 * @brief The missing wrapper compatible with hypre solver setup signature.
 * @param solver the solver
 * @param A the matrix
 * @param b the rhs vector (unused)
 * @param x the solution vector (unused)
 * @return hypre error code
 */
inline HYPRE_Int HYPRE_SLUDistSetup( HYPRE_Solver solver,
                                     HYPRE_ParCSRMatrix A,
                                     HYPRE_ParVector b,
                                     HYPRE_ParVector x )
{
  GEOSX_UNUSED_VAR( b, x )
  return hypre_SLUDistSetup( &solver, A, 0 );
}

/**
 * @brief The missing wrapper compatible with hypre solver solve signature.
 * @param solver the solver
 * @param A the matrix (unused)
 * @param b the rhs vector
 * @param x the solution vector
 * @return hypre error code
 */
inline HYPRE_Int HYPRE_SLUDistSolve( HYPRE_Solver solver,
                                     HYPRE_ParCSRMatrix A,
                                     HYPRE_ParVector b,
                                     HYPRE_ParVector x )
{
  GEOSX_UNUSED_VAR( A )
  return hypre_SLUDistSolve( solver, b, x );
}

/**
 * @brief The missing wrapper compatible with hypre solver destroy signature.
 * @param solver the solver
 * @return hypre error code
 */
inline HYPRE_Int HYPRE_SLUDistDestroy( HYPRE_Solver solver )
{
  return hypre_SLUDistDestroy( solver );
}

} // namespace hypre

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREUTILS_HPP_*/
