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
 * @brief SinglePhaseMixedMFD strategy: one-level multigrid reduction of the mixed mimetic
 *        saddle point
 *          [ M   D^T ]
 *          [ D    C  ]
 *        with face mass-flux (F) and cell pressure (C) unknowns.
 *
 * The interpolation is built from the diagonal of the flux block, W_p = -diag(M)^{-1} D^T,
 * the restriction is injection and the Galerkin coarse operator is the cell-centred
 * Laplacian S = C + D diag(M)^{-1} D^T: exact where both cells of a face use the diagonal
 * (TPFA) product, and spectrally equivalent to the pressure Schur complement elsewhere with
 * the constants of diag(M) ~ M. Being an M-matrix, S is what classical AMG is built for;
 * the F-relaxation is one symmetric Gauss-Seidel sweep on M and the coarse solve one BoomerAMG
 * V-cycle.
 *
 * The solver provides custom point markers (LinearSolverParameters::MGR::customPointMarkers):
 *  0 = face flux whose row is exactly diagonal (condensed two-point face, no-flow face)
 *  1 = face flux adjacent to at least one cell with the consistent (MFD) product
 *  2 = cell pressure
 * Both flux labels are F-points of the single reduction level.
 */
class SinglePhaseMixedMFD : public MGRStrategyBase< 1 >
{
public:
  /**
   * @brief Constructor.
   */
  explicit SinglePhaseMixedMFD( arrayView1d< int const > const & )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( 3 ) )
  {
    // Level 0: eliminate every face flux, keep the pressure
    m_labels[0].push_back( 2 );
    setupLabels();

    // l1-scaled relaxations are rejected by hypre here (the pressure rows have an empty C-C block
    // without accumulation) and plain Jacobi has no weight in MGR: symmetric Gauss-Seidel instead
    m_levelFRelaxType[0]         = MGRFRelaxationType::hybridSymmetricGaussSeidel;
    m_levelFRelaxIters[0]        = 1;
    m_levelInterpType[0]         = MGRInterpolationType::jacobi;
    m_levelRestrictType[0]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]   = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[0] = MGRGlobalSmootherType::none;
  }

  /**
   * @brief Setup the MGR strategy.
   * @param mgrParams MGR parameters
   * @param precond preconditioner wrapper
   * @param mgrData auxiliary MGR data
   */
  void setup( LinearSolverParameters::MGR const & mgrParams,
              HyprePrecWrapper & precond,
              HypreMGRData & mgrData )
  {
    GEOS_UNUSED_VAR( mgrParams );
    setReduction( precond, mgrData );

    // one V-cycle on the cell-centred Laplacian
    setPressureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEMIXEDMFD_HPP_*/
