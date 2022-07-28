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

#include "codingUtilities/Utilities.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>

#if GEOSX_USE_HYPRE_DEVICE == GEOSX_USE_HYPRE_CUDA || GEOSX_USE_HYPRE_DEVICE == GEOSX_USE_HYPRE_HIP
/// Host-device marker for custom hypre kernels
#define GEOSX_HYPRE_DEVICE GEOSX_DEVICE
#define GEOSX_HYPRE_HOST_DEVICE GEOSX_HOST_DEVICE
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
 * @brief Contains some hypre-specific functions.
 */
namespace hypre
{

/**
 * @brief @return Hypre memory location corresponding to a given LvArray memory space.
 * @param space the space
 */
constexpr HYPRE_MemoryLocation getMemoryLocation( LvArray::MemorySpace const space )
{
  switch( space )
  {
    case LvArray::MemorySpace::host: return HYPRE_MEMORY_HOST;
#if GEOSX_USE_HYPRE_DEVICE == GEOSX_USE_HYPRE_CUDA
    case LvArray::MemorySpace::cuda: return HYPRE_MEMORY_DEVICE;
#endif
#if GEOSX_USE_HYPRE_DEVICE == GEOSX_USE_HYPRE_HIP
  case LvArray::MemorySpace::hip: return HYPRE_MEMORY_DEVICE;
#endif
    default: return HYPRE_MEMORY_HOST;
  }
}

/**
 * @brief @return LvArray memory space corresponding to hypre's memory location.
 * @param location the location
 */
constexpr LvArray::MemorySpace getLvArrayMemorySpace( HYPRE_MemoryLocation const location )
{
  switch( location )
  {
    case HYPRE_MEMORY_HOST: return LvArray::MemorySpace::host;
#if GEOSX_USE_HYPRE_DEVICE == GEOSX_USE_HYPRE_CUDA
    case HYPRE_MEMORY_DEVICE: return LvArray::MemorySpace::cuda;
#endif
#if GEOSX_USE_HYPRE_DEVICE == GEOSX_USE_HYPRE_HIP
    case HYPRE_MEMORY_DEVICE: return LvArray::MemorySpace::hip;
#endif
    default: return LvArray::MemorySpace::host;
  }
}

#if GEOSX_USE_HYPRE_DEVICE == GEOSX_USE_HYPRE_CUDA
/// Execution policy for operations on hypre data
using execPolicy = parallelDevicePolicy<>;
/// Memory space used by hypre matrix/vector objects
constexpr LvArray::MemorySpace memorySpace = LvArray::MemorySpace::cuda;
/// Memory location used by hypre matrix/vector objects
constexpr HYPRE_MemoryLocation memoryLocation = HYPRE_MEMORY_DEVICE;
#elif GEOSX_USE_HYPRE_DEVICE == GEOSX_USE_HYPRE_HIP
/// Execution policy for operations on hypre data
using execPolicy = parallelDevicePolicy<>;
/// Memory space used by hypre matrix/vector objects
constexpr LvArray::MemorySpace memorySpace = LvArray::MemorySpace::hip;
/// Memory location used by hypre matrix/vector objects
constexpr HYPRE_MemoryLocation memoryLocation = HYPRE_MEMORY_DEVICE;
#else
/// Execution policy for operations on hypre data
using execPolicy = parallelHostPolicy;
/// Memory space used by hypre matrix/vector objects
constexpr LvArray::MemorySpace memorySpace = LvArray::MemorySpace::host;
/// Memory location used by hypre matrix/vector objects
constexpr HYPRE_MemoryLocation memoryLocation = HYPRE_MEMORY_HOST;
#endif

// Check matching requirements on index/value types between GEOSX and Hypre

// WARNING. We don't have consistent types between HYPRE_Int and localIndex.
//          Decision needs to be made either to use bigint option, or change
//          localIndex to int. We are getting away with this because we do not
//          pass ( localIndex * ) to hypre except when it is on the GPU, in
//          which case we are using int for localIndex.
#if GEOSX_USE_HYPRE_DEVICE == GEOSX_USE_HYPRE_CUDA || GEOSX_USE_HYPRE_DEVICE == GEOSX_USE_HYPRE_HIP
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
#if GEOSX_USE_HYPRE_DEVICE == GEOSX_USE_HYPRE_CUDA
  cudaError_t const err = cudaGetLastError();
  GEOSX_ERROR_IF( err != cudaSuccess, GEOSX_FMT( "Previous CUDA errors found: {} ({} at {}:{})", msg, cudaGetErrorString( err ), file, line ) );
#elif GEOSX_USE_HYPRE_DEVICE == GEOSX_USE_HYPRE_HIP
  hipError_t const err = hipGetLastError();
  GEOSX_ERROR_IF( err != hipSuccess, GEOSX_FMT( "Previous HIP errors found: {} ({} at {}:{})", msg, hipGetErrorString( err ), file, line ) );
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
 * @brief Converts a non-const array from GEOSX globalIndex type to HYPRE_BigInt.
 * @param[in] index the input array
 * @return the converted array
 */
inline HYPRE_BigInt * toHypreBigInt( geosx::globalIndex * const index )
{
  return reinterpret_cast< HYPRE_BigInt * >(index);
}

/**
 * @brief Converts a const array from GEOSX globalIndex type to HYPRE_BigInt.
 * @param[in] index the input array
 * @return the converted array
 */
inline HYPRE_BigInt const * toHypreBigInt( geosx::globalIndex const * const index )
{
  return reinterpret_cast< HYPRE_BigInt const * >(index);
}

/**
 * @brief Gather a parallel vector on a every rank.
 * @param vec the vector to gather
 * @return A newly allocated serial vector (may be null on ranks that don't have any elements)
 *
 * This is a wrapper around hypre_ParVectorToVectorAll() that works for both host-based
 * and device-based vectors without relying on Unified Memory.
 */
HYPRE_Vector parVectorToVectorAll( HYPRE_ParVector const vec );

/**
 * @brief Dummy function that does nothing but conform to hypre's signature for preconditioner setup/apply functions.
 * @return always 0 (success)
 *
 * Typical use is to prevent hypre from calling preconditioner setup when we have already called it on out side.
 */
HYPRE_Int dummySetup( HYPRE_Solver,
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
HYPRE_Int relaxationCreate( HYPRE_Solver & solver,
                            HYPRE_Int const type );

/**
 * @brief Setup a relaxation-based smoother.
 * @param solver the solver
 * @param A the matrix
 * @param b the rhs vector (unused)
 * @param x the solution vector (unused)
 * @return hypre error code
 */
HYPRE_Int relaxationSetup( HYPRE_Solver solver,
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
HYPRE_Int relaxationSolve( HYPRE_Solver solver,
                           HYPRE_ParCSRMatrix A,
                           HYPRE_ParVector b,
                           HYPRE_ParVector x );

/**
 * @brief Destroy a relaxation-based smoother.
 * @param solver the solver
 * @return always 0
 */
HYPRE_Int relaxationDestroy( HYPRE_Solver solver );

/**
 * @brief Returns hypre's identifier of the AMG cycle type.
 * @param type AMG cycle type
 * @return hypre AMG cycle type identifier code
 */
inline HYPRE_Int getAMGCycleType( LinearSolverParameters::AMG::CycleType const & type )
{
  static map< LinearSolverParameters::AMG::CycleType, HYPRE_Int > const typeMap =
  {
    { LinearSolverParameters::AMG::CycleType::V, 1 },
    { LinearSolverParameters::AMG::CycleType::W, 2 },
  };
  return findOption( typeMap, type, "multigrid cycle", "HyprePreconditioner" );
}

/**
 * @brief Returns hypre's identifier of the AMG smoother type.
 * @param type AMG smoother type
 * @return hypre AMG smoother type identifier code
 */
inline HYPRE_Int getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType const & type )
{
  static map< LinearSolverParameters::AMG::SmootherType, HYPRE_Int > const typeMap =
  {
    { LinearSolverParameters::AMG::SmootherType::default_, -1 },
#if GEOSX_USE_HYPRE_DEVICE == GEOSX_USE_HYPRE_CUDA || GEOSX_USE_HYPRE_DEVICE == GEOSX_USE_HYPRE_HIP
    { LinearSolverParameters::AMG::SmootherType::jacobi, 7 },
#else
    { LinearSolverParameters::AMG::SmootherType::jacobi, 0 },
#endif
    { LinearSolverParameters::AMG::SmootherType::fgs, 3 },
    { LinearSolverParameters::AMG::SmootherType::bgs, 4 },
    { LinearSolverParameters::AMG::SmootherType::sgs, 6 },
    { LinearSolverParameters::AMG::SmootherType::l1sgs, 8 },
    { LinearSolverParameters::AMG::SmootherType::chebyshev, 16 },
    { LinearSolverParameters::AMG::SmootherType::l1jacobi, 18 },
  };
  return findOption( typeMap, type, "multigrid relaxation", "HyprePreconditioner" );
}

/**
 * @brief Returns hypre's identifier of the AMG ILU smoother type.
 * @param type AMG ILU smoother type
 * @return hypre AMG ILU smoother type identifier code
 */
inline HYPRE_Int getILUType( LinearSolverParameters::AMG::SmootherType const type )
{
  static map< LinearSolverParameters::AMG::SmootherType, HYPRE_Int > const typeMap =
  {
    { LinearSolverParameters::AMG::SmootherType::ilu0, 0 },
    { LinearSolverParameters::AMG::SmootherType::ilut, 1 },
  };
  return findOption( typeMap, type, "ILU", "HyprePreconditioner" );
}

/**
 * @brief Returns hypre's identifier of the AMG coarse solver type.
 * @param type AMG coarse solver type
 * @return hypre AMG coarse solver type identifier code
 */
inline HYPRE_Int getAMGCoarseType( LinearSolverParameters::AMG::CoarseType const & type )
{
  static map< LinearSolverParameters::AMG::CoarseType, HYPRE_Int > const typeMap =
  {
    { LinearSolverParameters::AMG::CoarseType::default_, -1 },
#if GEOSX_USE_HYPRE_DEVICE == GEOSX_USE_HYPRE_CUDA || GEOSX_USE_HYPRE_DEVICE == GEOSX_USE_HYPRE_HIP
    { LinearSolverParameters::AMG::CoarseType::jacobi, 7 },
#else
    { LinearSolverParameters::AMG::CoarseType::jacobi, 0 },
#endif
    { LinearSolverParameters::AMG::CoarseType::fgs, 3 },
    { LinearSolverParameters::AMG::CoarseType::bgs, 4 },
    { LinearSolverParameters::AMG::CoarseType::sgs, 6 },
    { LinearSolverParameters::AMG::CoarseType::l1sgs, 8 },
    { LinearSolverParameters::AMG::CoarseType::direct, 9 },
    { LinearSolverParameters::AMG::CoarseType::chebyshev, 16 },
    { LinearSolverParameters::AMG::CoarseType::l1jacobi, 18 },
  };
  return findOption( typeMap, type, "multigrid coarse solver", "HyprePreconditioner" );
}

/**
 * @brief Returns hypre's identifier of the AMG coarsening type.
 * @param type AMG coarsening type
 * @return hypre AMG coarsening type identifier code
 */
inline HYPRE_Int getAMGCoarseningType( string const & type )
{
  static map< string, HYPRE_Int > const typeMap =
  {
    { "CLJP", 0 },
    { "Ruge-Stueben", 3 },
    { "Falgout", 6 },
    { "CLJPDebug", 7 },
    { "PMIS", 8 },
    { "PMISDebug", 9 },
    { "HMIS", 10 },
    { "CGC", 21 },
    { "CGC-E", 22 }
  };
  return findOption( typeMap, type, "multigrid coarsening", "HyprePreconditioner" );
}

/**
 * @brief Returns hypre's identifier of the relaxation preconditioner type.
 * @param type relaxation preconditioner type
 * @return hypre relaxation preconditioner type identifier code
 */
inline HYPRE_Int getRelaxationType( LinearSolverParameters::PreconditionerType const type )
{
  static map< LinearSolverParameters::PreconditionerType, HYPRE_Int > const typeMap =
  {
#if GEOSX_USE_HYPRE_DEVICE == GEOSX_USE_HYPRE_CUDA || GEOSX_USE_HYPRE_DEVICE == GEOSX_USE_HYPRE_HIP
    { LinearSolverParameters::PreconditionerType::jacobi, 7 },
#else
    { LinearSolverParameters::PreconditionerType::jacobi, 0 },
#endif
    { LinearSolverParameters::PreconditionerType::fgs, 3 },
    { LinearSolverParameters::PreconditionerType::bgs, 4 },
    { LinearSolverParameters::PreconditionerType::sgs, 6 },
    { LinearSolverParameters::PreconditionerType::l1sgs, 8 },
    { LinearSolverParameters::PreconditionerType::chebyshev, 16 },
    { LinearSolverParameters::PreconditionerType::l1jacobi, 18 },
  };
  return findOption( typeMap, type, "relaxation", "HyprePreconditioner" );
}

/**
 * @brief Returns hypre's identifier of the ILU preconditioner type.
 * @param type ILU preconditioner type
 * @return hypre ILU preconditioner type identifier code
 */
inline HYPRE_Int getILUType( LinearSolverParameters::PreconditionerType const type )
{
  static map< LinearSolverParameters::PreconditionerType, HYPRE_Int > const typeMap =
  {
    { LinearSolverParameters::PreconditionerType::iluk, 0 },
    { LinearSolverParameters::PreconditionerType::ilut, 1 },
  };
  return findOption( typeMap, type, "ILU", "HyprePreconditioner" );
}

/**
 * @enum AMGCoarseningType
 * @brief This enum class specifies the AMG parallel coarsening algorithm.
 */
enum class AMGCoarseningType : HYPRE_Int
{
  CLJP = 0,         //!< Parallel coarsening algorithm using independent sets
  Ruge_Stueben = 3, //!< Classical Ruge-Stueben coarsening on each processor
  Falgout = 6,      //!< Uses @p Ruge_Stueben first, followed by @p CLJP
  CLJPDebug = 7,    //!< Using a fixed random vector, for debugging purposes only
  PMIS = 8,         //!< Parallel coarsening algorithm using independent sets
  PMISDebug = 9,    //!< Using a fixed random vector, for debugging purposes only
  HMIS = 10,        //!< Uses one pass @p Ruge-Stueben on each processor independently, followed by @p PMIS
  CGC = 21,         //!< Coarsening by M. Griebel, B. Metsch and A. Schweitzer
  CGC_E = 22        //!< Coarsening by M. Griebel, B. Metsch and A. Schweitzer
};

/**
 * @enum MGRInterpolationType
 * @brief This enum class specifies the strategy for computing the level intepolation operator in MGR.
 */
enum class MGRInterpolationType : HYPRE_Int
{
  injection = 0,                      //!< Injection \f$[0  I]^{T}\f$
  l1jacobi = 1,                       //!< l1 Jacobi
  jacobi = 2,                         //!< Diagonal scaling
  classicalModifiedInterpolation = 3, //!< Classical modified interpolation
  approximateInverse = 4,             //!< Approximate inverse
  blockJacobi = 12                    //!< Block-Jacobi
};

/**
 * @enum MGRRestrictionType
 * @brief This enum class specifies the strategy for computing the level restriction operator in MGR.
 *
 */
enum class MGRRestrictionType : HYPRE_Int
{
  injection = 0,           //!< Injection \f$[0  I]\f$
  jacobi = 2,              //!< Diagonal scaling
  approximateInverse = 3,  //!< Approximate inverse
  blockJacobi = 12,        //!< Block-Jacobi
  cprLike = 13             //!< CPR-like restriction
};

/**
 * @enum MGRCoarseGridMethod
 * @brief This enum class specifies the strategy for level coarse grid computation in MGR.
 */
enum class MGRCoarseGridMethod : HYPRE_Int
{
  galerkin = 0,          //!< Galerkin coarse grid computation using RAP
  nonGalerkin = 1,       //!< Non-Galerkin coarse grid computation with dropping strategy: inv(A_FF) approximated by its (block) diagonal
                         //!< inverse
  cprLikeDiag = 2,       //!< Non-Galerkin coarse grid computation with dropping strategy: CPR-like approximation with inv(A_FF)
                         //!< approximated by its diagonal inverse
  cprLikeBlockDiag = 3,  //!< Non-Galerkin coarse grid computation with dropping strategy: CPR-like approximation with inv(A_FF)
                         //!< approximated by its block diagonal inverse
  approximateInverse = 4 //!< Non-Galerkin coarse grid computation with dropping strategy: inv(A_FF) approximated by sparse approximate
                         //!< inverse
};

/**
 * @enum MGRFRelaxationMethod
 * @brief This enum class specifies the F-relaxation strategy.
 */
enum class MGRFRelaxationMethod : HYPRE_Int
{
  singleLevel = 0, //!< single-level relaxation
  multilevel = 1,  //!< multilevel relaxation
  amgVCycle = 2    //!< multilevel relaxation
};

/**
 * @enum MGRFRelaxationType
 * @brief This enum class specifies the F-relaxation type.
 */
enum class MGRFRelaxationType : HYPRE_Int
{
  jacobi = 0,                       //!< Jacobi
  forwardHybridGaussSeidel = 3,     //!< hybrid Gauss-Seidel or SOR, forward solve
  backwardHybridGaussSeidel = 4,    //!< hybrid Gauss-Seidel or SOR, backward solve
  hybridSymmetricGaussSeidel = 6,   //!< hybrid symmetric Gauss-Seidel or SSOR
  l1hybridSymmetricGaussSeidel = 8, //!< \f$\ell_1\f$-scaled hybrid symmetric Gauss-Seidel
  l1forwardGaussSeidel = 13,        //!< \f$\ell_1\f$ Gauss-Seidel, forward solve
  l1backwardGaussSeidel = 14,       //!< \f$\ell_1\f$ Gauss-Seidel, backward solve
  l1jacobi = 18,                    //!< \f$\ell_1\f$-scaled Jacobi
};

/**
 * @enum MGRGlobalSmootherType
 * @brief This enum class specifies the global smoother type.
 */
enum class MGRGlobalSmootherType : HYPRE_Int
{
  blockJacobi = 0, //!< block Jacobi (default)
  jacobi = 1,      //!< Jacobi
  ilu0 = 16        //!< incomplete LU factorization
};

} // namespace hypre

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREUTILS_HPP_*/
