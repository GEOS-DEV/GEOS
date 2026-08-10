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

#include "common/MpiWrapper.hpp"
#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief SinglePhaseMixedMFD strategy: stencilFlag-guided three-level MGR reduction of the
 *        mixed mimetic finite difference saddle-point system
 *          [ M  -B^T ]
 *          [ B    0  ]
 *        with face mass-flux and cell pressure unknowns.
 *
 * The solver provides custom point markers (LinearSolverParameters::MGR::customPointMarkers)
 * splitting the face-flux dofs by the Global Adaptation classification:
 *  dofLabel: 0 = face flux whose adjacent cells are all TPFA-compatible (the assembled
 *               flux row is exactly diagonal: TPFA/TPFA interfaces, boundary faces of
 *               TPFA cells and no-flow identity rows)
 *  dofLabel: 1 = face flux adjacent to at least one MFD-compatible cell
 *  dofLabel: 2 = cell-centered pressure
 *
 * Ingredients:
 * 1. Level 0: F-points = TPFA-diagonal face fluxes (label 0). The F-block is exactly
 *    diagonal, so a single Jacobi sweep and Jacobi (diagonal) prolongation perform the
 *    elimination exactly, and the Galerkin (RAP) coarse grid is the exact Schur complement.
 *    The cost of this level is proportional to the TPFA fraction selected by the
 *    residual tolerance, mirroring the sparsity-reduction metric of the adaptive scheme.
 * 2. Level 1: F-points = MFD face fluxes (label 1), relaxed with symmetric Gauss-Seidel
 *    sweeps on the (SPD, well-conditioned) MFD flux block; Jacobi prolongation,
 *    injection restriction, Galerkin (RAP) coarse grid approximating the pressure
 *    Schur complement.
 * 3. Coarsest level: cell-pressure system solved with BoomerAMG.
 * 4. Global smoother: none. The Krylov solver is (F)GMRES.
 */
class SinglePhaseMixedMFD : public MGRStrategyBase< 2 >
{
public:
  /**
   * @brief Constructor.
   */
  explicit SinglePhaseMixedMFD( arrayView1d< int const > const & )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( 3 ) )
  {
    // Level 0: eliminate the TPFA-diagonal face fluxes, keep the MFD fluxes and the pressure
    m_labels[0].push_back( 1 );
    m_labels[0].push_back( 2 );

    // Level 1: eliminate the MFD face fluxes, keep the pressure
    m_labels[1].push_back( 2 );

    setupLabels();

    // Level 0: the F-block is exactly diagonal - one Jacobi sweep is an exact solve
    m_levelFRelaxType[0]         = MGRFRelaxationType::jacobi;
    m_levelFRelaxIters[0]        = 1;
    m_levelInterpType[0]         = MGRInterpolationType::jacobi;
    m_levelRestrictType[0]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]   = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[0] = MGRGlobalSmootherType::none;

    // Level 1: SGS sweeps on the well-conditioned SPD flux block outperform an AMG
    // V-cycle; Jacobi prolongation avoids the setup cost of approximateInverse
    m_levelFRelaxType[1]         = MGRFRelaxationType::hybridSymmetricGaussSeidel;
    m_levelFRelaxIters[1]        = 3;
    m_levelInterpType[1]         = MGRInterpolationType::jacobi;
    m_levelRestrictType[1]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1]   = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[1] = MGRGlobalSmootherType::none;
  }

  /**
   * @brief Setup the MGR strategy.
   * @param precond preconditioner wrapper
   * @param mgrData auxiliary MGR data
   */
  void setup( LinearSolverParameters::MGR const & mgrParams,
              HyprePrecWrapper & precond,
              HypreMGRData & mgrData )
  {
    // if the label-1 set is empty the level-1 reduction is the identity: use a single
    // level eliminating all (diagonal) fluxes at once
    integer localAnyLive = 0;
    for( localIndex i = 0; i < mgrParams.customPointMarkers.size(); ++i )
    {
      if( mgrParams.customPointMarkers[i] == 1 )
      {
        localAnyLive = 1;
        break;
      }
    }
    bool const anyLiveFaces = MpiWrapper::max( localAnyLive ) == 1;

    m_labels[0].clear();
    m_labels[1].clear();
    if( anyLiveFaces )
    {
      m_labels[0].push_back( 1 );
      m_labels[0].push_back( 2 );
      m_labels[1].push_back( 2 );
    }
    else
    {
      m_labels[0].push_back( 2 );
    }
    setupLabels();

    setReduction( precond, mgrData, anyLiveFaces ? numLevels : 1 );

    // Configure the BoomerAMG solver used as mgr coarse solver for the pressure Schur
    // complement. Two V-cycles per MGR application: the coarse solve accuracy governs the
    // outer FGMRES iteration count (a single cycle leaves the reduction quality unused)
    setPressureAMG( mgrData.coarseSolver );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( mgrData.coarseSolver.ptr, 2 ) );
    // strength threshold suited to 3D anisotropic permeability fields (hypre default 0.25
    // over-couples the weak direction and degrades coarsening quality)
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetStrongThreshold( mgrData.coarseSolver.ptr, 0.6 ) );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEMIXEDMFD_HPP_*/
