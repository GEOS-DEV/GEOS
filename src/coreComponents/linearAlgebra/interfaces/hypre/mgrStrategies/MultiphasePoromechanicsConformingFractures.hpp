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
 * @file SinglePhasePoromechanicsConformingFractures.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSMULIPHASEPROMECHANICSCONFORMINGFRACTURES_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSMULIPHASEPROMECHANICSCONFORMINGFRACTURES_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief MultiPhasePoromechanicsConformingFractures strategy.
 *
 * dofLabel: 0 = displacement, x-component
 * dofLabel: 1 = displacement, y-component
 * dofLabel: 2 = displacement, z-component
 * dofLabel: 3 = face-centered lagrange multiplier (tn)
 * dofLabel: 4 = face-centered lagrange multiplier (tt1)
 * dofLabel: 5 = face-centered lagrange multiplier (tt2)
 * dofLabel: 6 = pressure (cell elem + fracture elems)
 * dofLabel: 7 = density (cell elem + fracture elems)
 * dofLabel: ... = densities
 * dofLabel: numLabels - 1 = largest density / vol constaint
 *
 * Ingredients:
 * 1. Level 0: F-points lag mult (3,4,5), C-points displacements, pressure and densities (0,1,2,6,7...)
 * 2. Level 1: F-points displacement (0,1,2), C-points pressure and densities (6,7...)
 * 3. Level 2: F-points volume constaint (n) and C-points pressure (6,7...,n-1)
 * 4. Level 3: F-points densities (7,..,n-1) and C-points pressure (6)
 * 6. Global smoother: none
 */
class MultiphasePoromechanicsConformingFractures : public MGRStrategyBase< 4 >
{
public:

  /**
   * @brief Constructor.
   */
  explicit MultiphasePoromechanicsConformingFractures( arrayView1d< int const > const & numComponentsPerField)
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] + numComponentsPerField[1] + numComponentsPerField[2] ) )
  {
    //totalDisplacement (3) comes first, then contact fields traction (3) and eventually compositional variables (3+)

    // we keep u, dens and p - elim lag mult
    m_labels[0].resize( m_numBlocks - numComponentsPerField[1] );
    HYPRE_Int const numResLabels = LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[2] );

    m_labels[0][0]= 0;
    m_labels[0][1]= 1;
    m_labels[0][2]= 2;
    m_labels[0][3]= 6;//pressure
    for(int i=4, c=7; i<int(m_labels[0].size()); i++, c++)
      m_labels[0][i] = c;
  
    // we keep dens and p - elim u 
    m_labels[1].resize( numResLabels );
    m_labels[1][0] = 6;
    for(int i=1, c=7; i<int(m_labels[1].size()); i++, c++)
      m_labels[1][i] = c;
    // we keep p - elim total density
    m_labels[2].resize( numResLabels - 1);
    m_labels[2][0] = 6;
    for(int i=1, c=7; i<int(m_labels[2].size()); i++, c++)
      m_labels[2][i] = c;

    // we keep only p
    m_labels[3].push_back( 6 );


    setupLabels();

    // Level 0 - lag mult
    // m_levelFRelaxType[0]          = MGRFRelaxationType::l1jacobi;
    // m_levelFRelaxIters[0]         = 1;
    m_levelFRelaxType[0]          = MGRFRelaxationType::none;
    m_levelFRelaxIters[0]         = 0;

    m_levelInterpType[0]          = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[0]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[0]  = MGRGlobalSmootherType::none;
    // m_levelGlobalSmootherType[0]  = MGRGlobalSmootherType::ilu0;
    // m_levelGlobalSmootherIters[0] = 1;
    
    // Level 1 - displacement
    m_levelFRelaxType[1]          = MGRFRelaxationType::amgVCycle;
    m_levelFRelaxIters[1]         = 1;
    m_levelInterpType[1]          = MGRInterpolationType::jacobi;
    m_levelRestrictType[1]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1]    = MGRCoarseGridMethod::nonGalerkin;
    m_levelGlobalSmootherType[1]  = MGRGlobalSmootherType::none;
    
    // Level 2 - total densities 
    m_levelFRelaxType[2]          = MGRFRelaxationType::jacobi;
    m_levelFRelaxIters[2]         = 1;
    m_levelInterpType[2]          = MGRInterpolationType::jacobi;
    m_levelRestrictType[2]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[2]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[2]  = MGRGlobalSmootherType::none;

    // Level 3 - remaining densities
    m_levelFRelaxType[3]          = MGRFRelaxationType::none;
    m_levelInterpType[3]          = MGRInterpolationType::injection;
    m_levelRestrictType[3]        = MGRRestrictionType::blockColLumped; // True-IMPES
    m_levelCoarseGridMethod[3]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[3]  = MGRGlobalSmootherType::ilu0;
    m_levelGlobalSmootherIters[3] = 1;

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
    setReduction( precond, mgrData );
    setDisplacementAMG( mgrData.mechSolver, mgrParams.separateComponents );
    setPressureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre


}
#endif /* GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSMULIPHASEPROMECHANICSCONFORMINGFRACTURES_HPP_ */
