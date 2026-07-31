/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseMixedMFD.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEMIXEDMFD_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEMIXEDMFD_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief SinglePhaseMixedMFD strategy.
 *
 * Two-level MGR reduction for the mixed mimetic finite difference saddle-point system
 *   [ M  -B^T ]
 *   [ B    0  ]
 * with face mass-flux and cell pressure unknowns.
 *
 * Labels description stored in point_marker_array
 *  dofLabel: 0 = face-centered mass flux
 *  dofLabel: 1 = cell-centered pressure
 *
 * Ingredients:
 * 1. F-points face-centered mass flux (the SPD adaptive inner-product block M),
 *    C-points cell-centered pressure
 * 2. F-points smoother: single BoomerAMG V-cycle on the flux block
 * 3. Restriction: injection; prolongation: Jacobi; coarse grid: Galerkin (RAP),
 *    yielding an approximate pressure Schur complement
 * 4. C-points coarse-grid/Schur complement solver: BoomerAMG
 * 5. Global smoother: none
 *
 * The corresponding Krylov solver is (F)GMRES. This mirrors the hypredrive recipe used
 * in the adaptive-MFD reference implementation (two-level MGR with f_dofs = flux labels,
 * AMG F-relaxation and an AMG-solved RAP coarse operator).
 */
class SinglePhaseMixedMFD : public MGRStrategyBase< 1 >
{
public:
  /**
   * @brief Constructor.
   */
  explicit SinglePhaseMixedMFD( arrayView1d< int const > const & )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( 2 ) )
  {
    // Level 0: eliminate the face-centered mass fluxes, keep the cell-centered pressure
    m_labels[0].push_back( 1 );

    setupLabels();

    // Level 0
    m_levelFRelaxType[0]         = MGRFRelaxationType::amgVCycle;
    m_levelFRelaxIters[0]        = 1;
    m_levelInterpType[0]         = MGRInterpolationType::jacobi;
    m_levelRestrictType[0]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]   = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[0] = MGRGlobalSmootherType::none;
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

    // Configure the BoomerAMG solver used as mgr coarse solver for the pressure Schur complement
    setPressureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEMIXEDMFD_HPP_*/
