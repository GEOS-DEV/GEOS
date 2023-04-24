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
 * @file CompositionalMultiphaseFVM.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASEFVM_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASEFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief CompositionalMultiphaseFVM strategy.
 *
 * Labels description stored in point_marker_array
 *               0 = pressure
 *               1 = density
 *             ... = densities
 *   numLabels - 1 = density
 *
 * 2-level MGR reduction strategy:
 *   - 1st level: eliminate the reservoir density associated with the volume constraint
 *   - 2nd level: eliminate the other reservoir densities
 *   - The coarse grid (pressure system) is solved with BoomerAMG.
 *
 */
class CompositionalMultiphaseFVM : public MGRStrategyBase< 2 >
{
public:
  /**
   * @brief Constructor.
   * @param numComponentsPerField array with number of components for each field
   */
  explicit CompositionalMultiphaseFVM( arrayView1d< int const > const & numComponentsPerField )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] ) )
  {
    // Level 0: eliminate last density which corresponds to the volume constraint equation
    m_labels[0].resize( m_numBlocks - 1 );
    std::iota( m_labels[0].begin(), m_labels[0].end(), 0 );
    // Level 1: eliminate the remaining the densities
    m_labels[1].push_back( 0 );

    setupLabels();

    m_levelFRelaxMethod[0]     = MGRFRelaxationMethod::singleLevel; //default, i.e. Jacobi
    m_levelInterpType[0]       = MGRInterpolationType::jacobi;
    m_levelRestrictType[0]     = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0] = MGRCoarseGridMethod::galerkin;

    m_levelFRelaxMethod[1]     = MGRFRelaxationMethod::singleLevel; //default, i.e. Jacobi
    m_levelInterpType[1]       = MGRInterpolationType::injection;
    m_levelRestrictType[1]     = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1] = MGRCoarseGridMethod::cprLikeBlockDiag;

    // ILU smoothing for the system made of pressure and densities (except the last one)
    m_levelGlobalSmootherType[1]  = 16;
    m_levelGlobalSmootherIters[1] = 1;
  }

  /**
   * @brief Setup the MGR strategy.
   * @param precond preconditioner wrapper
   * @param mgrData auxiliary MGR data
   */
  void setup( LinearSolverParameters::MGR const &,
              HyprePrecWrapper & precond,
              HypreMGRData & mgrData )
  {
    setReduction( precond, numLevels, mgrData );

    // Note: uncommenting HYPRE_MGRSetLevelFRelaxMethod and commenting HYPRE_MGRSetRelaxType breaks the recipe. This requires further
    // investigation
    //GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxMethod( precond.ptr, toUnderlyingPtr( m_levelFRelaxMethod ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetRelaxType( precond.ptr, 0 ));
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetNumRelaxSweeps( precond.ptr, 1 ));
#ifdef GEOSX_USE_HYPRE_CUDA
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetRelaxType( precond.ptr, getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::l1jacobi ) ) );
#endif
    
    // Configure the BoomerAMG solver used as mgr coarse solver for the pressure reduced system
    setPressureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASEFVM_HPP_*/
