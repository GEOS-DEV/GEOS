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
 * @brief SinglePhaseMixedMFD strategy: two-level MGR reduction of the mixed mimetic
 *        finite difference saddle-point system
 *          [ M  -B^T ]
 *          [ B    0  ]
 *        with face mass-flux and cell pressure unknowns.
 *
 * The solver supplies three labels in point_marker_array:
 *  dofLabel: 0 = face flux whose adjacent cells are all TPFA-compatible (an exactly
 *                diagonal assembled flux row)
 *  dofLabel: 1 = face flux adjacent to at least one MFD-compatible cell
 *  dofLabel: 2 = cell-centered pressure
 *
 * Ingredients:
 * 1. F-points = both face-flux classes (labels 0 and 1), eliminated together in a
 *    single reduction level and relaxed with one Jacobi sweep on the complete flux block.
 * 2. C-points = cell pressures (label 2); Jacobi prolongation, injection restriction,
 *    and a Galerkin (RAP) coarse grid form an approximate pressure Schur complement.
 * 3. The pressure coarse system is solved with BoomerAMG.
 * 4. Global (G) relaxation: none. The Krylov solver is (F)GMRES.
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
    // Eliminate both face-flux labels together and keep cell pressure.
    m_labels[0].push_back( 2 );

    setupLabels();

    // The RT0-like face-flux mass block is well-conditioned, so use one inexpensive
    // Jacobi sweep for F-relaxation.
    m_levelFRelaxType[0]         = MGRFRelaxationType::jacobi;
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

    // Use two BoomerAMG V-cycles for the scalar pressure Schur complement. The stronger
    // coarse solve targets the slowly converging pressure modes without doubling the cost
    // of each FGMRES iteration.
    // No numFunctions override is needed for this scalar coarse system.
    setPressureAMG( mgrData.coarseSolver );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetAggNumLevels( mgrData.coarseSolver.ptr, 0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( mgrData.coarseSolver.ptr,
                                                       hypre::getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::l1sgs ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( mgrData.coarseSolver.ptr, 2 ) );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEMIXEDMFD_HPP_*/
