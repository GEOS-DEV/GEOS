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
 * @file SinglePhasePoromechanicsEmbeddedFractures.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEPOROMECHANICSEMBEDDEDFRACTURES_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEPOROMECHANICSEMBEDDEDFRACTURES_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief SinglePhasePoromechanicsEmbeddedFractures strategy.
 *
 * dofLabel: 0 = displacement, x-component
 * dofLabel: 1 = displacement, y-component
 * dofLabel: 2 = displacement, z-component
 * dofLabel: 3 = wn
 * dofLabel: 4 = wt1
 * dofLabel: 5 = wt2
 * dofLabel: 6 = pressure (cell elem + fracture elems)
 *
 * Ingredients:
 * 1. F-points displacement (0,1,2), C-points pressure (3)
 * 2. F-points smoother: AMG, single V-cycle, separate displacement components
 * 3. C-points coarse-grid/Schur complement solver: boomer AMG
 * 4. Global smoother: none
 */
class SinglePhasePoromechanicsEmbeddedFractures : public MGRStrategyBase< 2 >
{
public:

  /**
   * @brief Constructor.
   */
  explicit SinglePhasePoromechanicsEmbeddedFractures( arrayView1d< int const > const & )
    : MGRStrategyBase( 7 )
  {

    /// IDEAL strategy but it seg faults
    // // we keep u and p
    // m_labels[0].push_back( 0 );
    // m_labels[0].push_back( 1 );
    // m_labels[0].push_back( 2 );
    // m_labels[0].push_back( 6 );
    // // we keep p
    // m_labels[1].push_back( 6 );

    // setupLabels();

    // // Level 0
    // m_levelFRelaxMethod[0]     = MGRFRelaxationMethod::singleLevel;
    // m_levelFRelaxType[0]       = MGRFRelaxationType::gsElimWInverse; // gaussian elimination for the dispJump block
    // m_levelInterpType[0]       = MGRInterpolationType::blockJacobi;
    // m_levelRestrictType[0]     = MGRRestrictionType::injection;
    // m_levelCoarseGridMethod[0] = MGRCoarseGridMethod::galerkin;

    // // Level 1
    // m_levelFRelaxMethod[1]     = MGRFRelaxationMethod::amgVCycle;
    // m_levelInterpType[1]       = MGRInterpolationType::jacobi;
    // m_levelRestrictType[1]     = MGRRestrictionType::injection;
    // m_levelCoarseGridMethod[1] = MGRCoarseGridMethod::nonGalerkin;

    // we keep w and p
    m_labels[0].push_back( 3 );
    m_labels[0].push_back( 4 );
    m_labels[0].push_back( 5 );
    m_labels[0].push_back( 6 );
    // we keep p
    m_labels[1].push_back( 6 );

    setupLabels();

    // Level 0
    m_levelFRelaxType[0]         = MGRFRelaxationType::amgVCycle;
    m_levelFRelaxIters[0]        = 1;
    m_levelInterpType[0]         = MGRInterpolationType::jacobi;
    m_levelRestrictType[0]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]   = MGRCoarseGridMethod::nonGalerkin;
    m_levelGlobalSmootherType[0] = MGRGlobalSmootherType::none;

    // Level 1
    m_levelFRelaxType[1]         = MGRFRelaxationType::gsElimWInverse;
    m_levelFRelaxIters[1]        = 1;
    m_levelInterpType[1]         = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[1]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1]   = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[1] = MGRGlobalSmootherType::none;
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
    setReduction( precond, mgrData );

    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetPMaxElmts( precond.ptr, 0 ));

    // Configure the BoomerAMG solver used as F-relaxation for the first level
    setMechanicsFSolver( precond, mgrData );

    // Configure the BoomerAMG solver used as mgr coarse solver for the pressure reduced system
    setPressureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEPOROMECHANICSEMBEDDEDFRACTURES_HPP_*/
