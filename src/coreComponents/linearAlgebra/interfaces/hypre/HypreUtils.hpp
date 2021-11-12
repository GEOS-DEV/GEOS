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
 * @file HypreUtils.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREUTILS_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREUTILS_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>

#ifdef GEOSX_USE_HYPRE_CUDA
/// Host-device marker for custom hypre kernels
#define GEOSX_HYPRE_DEVICE GEOSX_DEVICE
#else
/// Host-device marker for custom hypre kernels
#define GEOSX_HYPRE_DEVICE
#endif

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

#ifdef GEOSX_USE_HYPRE_CUDA
/// Execution policy for operations on hypre data
using execPolicy = parallelDevicePolicy<>;
/// Memory space used by hypre matrix/vector objects
constexpr LvArray::MemorySpace memorySpace = LvArray::MemorySpace::cuda;
constexpr HYPRE_MemoryLocation memoryLocation = HYPRE_MEMORY_DEVICE;
#else
/// Execution policy for operations on hypre data
using execPolicy = parallelHostPolicy;
/// Memory space used by hypre matrix/vector objects
constexpr LvArray::MemorySpace memorySpace = LvArray::MemorySpace::host;
constexpr HYPRE_MemoryLocation memoryLocation = HYPRE_MEMORY_HOST;
#endif

// Check matching requirements on index/value types between GEOSX and Hypre

// WARNING. We don't have consistent types between HYPRE_Int and localIndex.
//          Decision needs to be made either to use bigint option, or change
//          localIndex to int. We are getting away with this because we do not
//          pass ( localIndex * ) to hypre except when it is on the GPU, in
//          which case we are using int for localIndex.
#ifdef GEOSX_USE_HYPRE_CUDA
static_assert( sizeof( HYPRE_Int ) == sizeof( geosx::localIndex ),
               "HYPRE_Int and geosx::localIndex must have the same size" );
static_assert( std::is_signed< HYPRE_Int >::value == std::is_signed< geosx::localIndex >::value,
               "HYPRE_Int and geosx::localIndex must both be signed or unsigned" );
#endif

/**
 * @brief
 * @param msg
 * @param file
 * @param line
 */
inline void checkDeviceErrors( char const * msg, char const * file, int const line )
{
#ifdef GEOSX_USE_HYPRE_CUDA
  cudaError_t const err = cudaGetLastError();
  GEOSX_ERROR_IF( err != cudaSuccess, GEOSX_FMT( "Previous CUDA errors found: {} ({} at {}:{})", msg, cudaGetErrorString( err ), file, line ) );
#else
  GEOSX_UNUSED_VAR( msg, file, line );
#endif
}

/**
 * @brief Check for previous device errors and report with line information.
 * @param msg custom message to add
 */
#define GEOSX_HYPRE_CHECK_DEVICE_ERRORS( msg ) ::geosx::hypre::checkDeviceErrors( msg, __FILE__, __LINE__ )

static_assert( sizeof( HYPRE_BigInt ) == sizeof( geosx::globalIndex ),
               "HYPRE_BigInt and geosx::globalIndex must have the same size" );

static_assert( std::is_signed< HYPRE_BigInt >::value == std::is_signed< geosx::globalIndex >::value,
               "HYPRE_BigInt and geosx::globalIndex must both be signed or unsigned" );

static_assert( std::is_same< HYPRE_Real, geosx::real64 >::value,
               "HYPRE_Real and geosx::real64 must be the same type" );

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
HYPRE_Int DummySetup( HYPRE_Solver,
                      HYPRE_ParCSRMatrix,
                      HYPRE_ParVector,
                      HYPRE_ParVector );

/**
 * @brief The missing wrapper compatible with hypre solver solve signature.
 * @param solver the solver
 * @param A the matrix (unused)
 * @param b the rhs vector
 * @param x the solution vector
 * @return hypre error code
 */
HYPRE_Int SuperLUDistSolve( HYPRE_Solver solver,
                            HYPRE_ParCSRMatrix A,
                            HYPRE_ParVector b,
                            HYPRE_ParVector x );

/**
 * @brief The missing wrapper compatible with hypre solver destroy signature.
 * @param solver the solver
 * @return hypre error code
 */
HYPRE_Int SuperLUDistDestroy( HYPRE_Solver solver );

/**
 * @brief Create a relaxation-based smoother.
 * @param solver the solver
 * @param type hypre's internal identifier of the relaxation type
 * @return always 0
 */
HYPRE_Int RelaxationCreate( HYPRE_Solver & solver,
                            HYPRE_Int const type );

/**
 * @brief Setup a relaxation-based smoother.
 * @param solver the solver
 * @param A the matrix
 * @param b the rhs vector (unused)
 * @param x the solution vector (unused)
 * @return hypre error code
 */
HYPRE_Int RelaxationSetup( HYPRE_Solver solver,
                           HYPRE_ParCSRMatrix A,
                           HYPRE_ParVector b,
                           HYPRE_ParVector x );

/**
 * @brief Solve with a relaxation-based smoother.
 * @param solver the solver
 * @param A the matrix
 * @param b the rhs vector (unused)
 * @param x the solution vector (unused)
 * @return hypre error code
 */
HYPRE_Int RelaxationSolve( HYPRE_Solver solver,
                           HYPRE_ParCSRMatrix A,
                           HYPRE_ParVector b,
                           HYPRE_ParVector x );

/**
 * @brief Destroy a relaxation-based smoother.
 * @param solver the solver
 * @return always 0
 */
HYPRE_Int RelaxationDestroy( HYPRE_Solver solver );

} // namespace hypre

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREUTILS_HPP_*/
