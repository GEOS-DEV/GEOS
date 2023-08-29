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
 * @file HypreMGR.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRSTRATEGIES_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRSTRATEGIES_HPP_

#include "linearAlgebra/common/common.hpp"

#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"


#include <_hypre_utilities.h>

namespace geosx
{

/**
 * @brief Container for hypre preconditioner auxiliary data for MGR.
 */
struct HypreMGRData
{
  array1d< HYPRE_Int > pointMarkers;  ///< array1d of unique tags for local degrees of freedom
  HyprePrecWrapper coarseSolver;      ///< MGR coarse solver pointer and functions
  HyprePrecWrapper mechSolver;        ///< MGR mechanics fine solver pointer and functions
};

namespace hypre
{

namespace mgr
{

/**
 * @brief Helper to simplify MGR setup
 * @tparam STRATEGY strategy class (one of structs defined below)
 * @param params MGR parameters
 * @param numComponentsPerField array containing number of components in each field
 * @param precond the preconditioner wrapper
 * @param mgrData additional MGR data struct with marker array populated
 */
template< typename STRATEGY >
void setStrategy( LinearSolverParameters::MGR const & params,
                  arrayView1d< int const > const & numComponentsPerField,
                  HyprePrecWrapper & precond,
                  HypreMGRData & mgrData )
{
  STRATEGY strategy( numComponentsPerField );
  strategy.setup( params, precond, mgrData );
}

/**
 * @brief Helper struct for strategies that provides some basic parameter arrays needed by MGR.
 * @tparam NLEVEL number of reduction levels (not including the coarsest level)
 */
template< int NLEVEL >
class MGRStrategyBase
{
public:

  static constexpr HYPRE_Int numLevels = NLEVEL;       ///< Number of levels

protected:

  HYPRE_Int m_numBlocks{ 0 };                          ///< Number of different matrix blocks treated separately

  std::vector< HYPRE_Int > m_labels[numLevels]{};      ///< Dof labels kept at each level
  HYPRE_Int m_numLabels[numLevels]{ -1 };              ///< Number of dof labels kept
  HYPRE_Int * m_ptrLabels[numLevels]{ nullptr };       ///< Pointers to each level's labels, as consumed by MGR

  MGRFRelaxationMethod m_levelFRelaxMethod[numLevels];    ///< F-relaxation strategy for each level (single or multilevel)
  MGRFRelaxationType m_levelFRelaxType[numLevels];        ///< F-relaxation type for each level
  MGRInterpolationType m_levelInterpType[numLevels];      ///< Interpolation type for each level
  MGRRestrictionType m_levelRestrictType[numLevels];      ///< Restriction type for each level
  MGRCoarseGridMethod m_levelCoarseGridMethod[numLevels]; ///< Coarse grid method for each level
  HYPRE_Int m_levelSmoothType[numLevels]{ -1 };           ///< Smoother type for each level
  HYPRE_Int m_levelSmoothIters[numLevels]{ -1 };          ///< Number of smoother iterations for each level

  HYPRE_Int m_numRestrictSweeps{ -1 }; ///< Number of restrict sweeps
  HYPRE_Int m_numInterpSweeps{ -1 };   ///< Number of interpolation sweeps

  HYPRE_Int m_numRelaxSweeps{ -1 }; ///< F-relaxation number of sweeps

  MGRGlobalSmootherType m_globalSmoothType{ MGRGlobalSmootherType::blockJacobi }; ///< Global smoothing type
  HYPRE_Int m_numGlobalSmoothSweeps{ -1 }; ///< Global smoothing number of iterations

  /**
   * @brief Constructor.
   * @param numBlocks number of blocks
   */
  explicit MGRStrategyBase( HYPRE_Int const numBlocks )
    : m_numBlocks( numBlocks )
  {
    for( HYPRE_Int i = 0; i < numLevels; ++i )
    {
      m_levelFRelaxMethod[i]     = MGRFRelaxationMethod::singleLevel;
      m_levelFRelaxType[i]       = MGRFRelaxationType::jacobi;
      m_levelInterpType[i]       = MGRInterpolationType::jacobi;
      m_levelRestrictType[i]     = MGRRestrictionType::injection;
      m_levelCoarseGridMethod[i] = MGRCoarseGridMethod::galerkin;
    }
  }

  /**
   * @brief Call this after populating lv_cindexes.
   */
  void setupLabels()
  {
    for( HYPRE_Int i = 0; i < numLevels; ++i )
    {
      m_numLabels[i] = m_labels[i].size();
      m_ptrLabels[i] = m_labels[i].data();
    }
  }

  /**
   * @brief Set up BoomerAMG to perform the mechanics F-solve for the first F-relaxation
   * @param precond the preconditioner wrapper
   * @param mgrData auxiliary MGR data
   */
  void setMechanicsFSolver( HyprePrecWrapper & precond,
                            HypreMGRData & mgrData ) const
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &mgrData.mechSolver.ptr ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( mgrData.mechSolver.ptr, 0.0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( mgrData.mechSolver.ptr, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxRowSum( mgrData.mechSolver.ptr, 1.0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetStrongThreshold( mgrData.mechSolver.ptr, 0.8 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( mgrData.mechSolver.ptr, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( mgrData.mechSolver.ptr, 0 ) );

#ifdef GEOSX_USE_HYPRE_CUDA
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCoarsenType( mgrData.mechSolver.ptr, hypre::getAMGCoarseningType( "PMIS" ) ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( mgrData.mechSolver.ptr, hypre::getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::chebyshev ) ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( mgrData.mechSolver.ptr, 2 ) );
#else
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( mgrData.mechSolver.ptr, 1 ) );
#endif

    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumFunctions( mgrData.mechSolver.ptr, 3 ) );

    mgrData.mechSolver.setup = HYPRE_BoomerAMGSetup;
    mgrData.mechSolver.solve = HYPRE_BoomerAMGSolve;
    mgrData.mechSolver.destroy = HYPRE_BoomerAMGDestroy;

    HYPRE_MGRSetFSolver( precond.ptr, mgrData.mechSolver.solve, mgrData.mechSolver.setup, mgrData.mechSolver.ptr );
  }

};

/**
 * @brief Create the MGR preconditioner object.
 * @param params preconditioner parameters
 * @param dofManager pointer to DofManager for the linear system
 * @param precond the preconditioner
 * @param mgrData auxiliary data for MGR
 */
void createMGR( LinearSolverParameters const & params,
                DofManager const * const dofManager,
                HyprePrecWrapper & precond,
                HypreMGRData & mgrData );

} // namespace mgr

} // namespace hypre

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRSTRATEGIES_HPP_*/
