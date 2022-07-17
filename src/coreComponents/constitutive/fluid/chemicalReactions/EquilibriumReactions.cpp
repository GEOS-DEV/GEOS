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
 * @file EquilibriumReactions.cpp
 */

#include "EquilibriumReactions.hpp"

#include "functions/FunctionManager.hpp"
#include "linearAlgebra/interfaces/dense/BlasLapackLA.hpp"

namespace geosx
{

using namespace stringutilities;

namespace constitutive
{

namespace chemicalReactions
{

EquilibriumReactions::EquilibriumReactions( string const & name, integer const numPrimarySpecies, integer const numSecSpecies ):
  ReactionsBase( name, numPrimarySpecies, numSecSpecies )
{
  // Here we should either read the database or the input values.

  // Hardcoding values for now

  // Activity coefficient related constants
  m_DebyeHuckelA = 0.5465;
  m_DebyeHuckelB = 0.3346;
  m_WATEQBDot = 0.0438;

  m_ionSizePrimary.resize( m_numPrimarySpecies );
  m_ionSizePrimary[0] = 9.00;
  m_ionSizePrimary[1] = 4.00;
  m_ionSizePrimary[2] = 6.00;
  m_ionSizePrimary[3] = 4.00;
  m_ionSizePrimary[4] = 3.00;
  m_ionSizePrimary[5] = 8.00;
  m_ionSizePrimary[6] = 4.00;

  m_ionSizeSec.resize( m_numSecondarySpecies );
  m_ionSizeSec[0] = 3.50;
  m_ionSizeSec[1] = 3.00;
  m_ionSizeSec[2] = 4.50;
  m_ionSizeSec[3] = 3.00;
  m_ionSizeSec[4] = 4.00;
  m_ionSizeSec[5] = 3.00;
  m_ionSizeSec[6] = 3.00;
  m_ionSizeSec[7] = 4.00;
  m_ionSizeSec[8] = 3.00;
  m_ionSizeSec[9] = 3.00;
  m_ionSizeSec[10] = 4.00;

  m_chargePrimary.resize( m_numPrimarySpecies );
  m_chargePrimary[0] = 1;
  m_chargePrimary[1] = -1;
  m_chargePrimary[2] = 2;
  m_chargePrimary[3] = -2;
  m_chargePrimary[4] = -1;
  m_chargePrimary[5] = 2;
  m_chargePrimary[6] = 1;

  m_chargeSec.resize( m_numSecondarySpecies );
  m_chargeSec[0] = -1;
  m_chargeSec[1] = 0;
  m_chargeSec[2] = -2;
  m_chargeSec[3] = 0;
  m_chargeSec[4] = 1;
  m_chargeSec[5] = 0;
  m_chargeSec[6] = 0;
  m_chargeSec[7] = 1;
  m_chargeSec[8] = 0;
  m_chargeSec[9] = 0;
  m_chargeSec[10] = -1;

  // Stochiometric Matrix
  // First index: 0 = OH-, 1 = CO2, 2 = CO3-2, 3 = H2CO3, 4 = CaHCO3+, 5 = CaCO3, 6 = CaSO4, 7 = CaCl+, 8 = CaCl2, 9 = MgSO4, 10 = NaSO4-
  // Second index: 0 = H+, 1 = HCO3-, 2 = Ca+2, 3 = SO4-2, 4 = Cl-, 5 = Mg+2, 6 = Na+1
  m_stoichMatrix.resize( m_numSecondarySpecies, m_numPrimarySpecies );
  m_stoichMatrix[0][0] = -1;
  m_stoichMatrix[1][0] = 1;
  m_stoichMatrix[1][1] = 1;
  m_stoichMatrix[2][0] = -1;
  m_stoichMatrix[2][1] = 1;
  m_stoichMatrix[3][0] = 1;
  m_stoichMatrix[3][1] = 1;
  m_stoichMatrix[4][1] = 1;
  m_stoichMatrix[4][2] = 1;
  m_stoichMatrix[5][0] = -1;
  m_stoichMatrix[5][1] = 1;
  m_stoichMatrix[5][2] = 1;
  m_stoichMatrix[6][2] = 1;
  m_stoichMatrix[6][3] = 1;
  m_stoichMatrix[7][2] = 1;
  m_stoichMatrix[7][4] = 1;
  m_stoichMatrix[8][2] = 1;
  m_stoichMatrix[8][4] = 2;
  m_stoichMatrix[9][5] = 1;
  m_stoichMatrix[9][3] = 1;
  m_stoichMatrix[10][6] = 1;
  m_stoichMatrix[10][3] = 1;

  // Equilibrium Constant
  m_log10EqConst.resize( m_numSecondarySpecies );
  m_log10EqConst[0] = 13.99;
  m_log10EqConst[1] = -6.36;
  m_log10EqConst[2] = 10.33;
  m_log10EqConst[3] = -3.77;
  m_log10EqConst[4] = -1.09;
  m_log10EqConst[5] = 7.07;
  m_log10EqConst[6] = -2.16;
  m_log10EqConst[7] = 0.67;
  m_log10EqConst[8] = 0.60;
  m_log10EqConst[9] = -2.43;
  m_log10EqConst[10] = -0.82;
}

EquilibriumReactions::KernelWrapper EquilibriumReactions::createKernelWrapper() const
{
  return KernelWrapper( m_numPrimarySpecies,
                        m_numSecondarySpecies,
                        m_log10EqConst,
                        m_stoichMatrix,
                        m_chargePrimary,
                        m_chargeSec,
                        m_ionSizePrimary,
                        m_ionSizeSec,
                        m_DebyeHuckelA,
                        m_DebyeHuckelB,
                        m_WATEQBDot );
}


GEOSX_HOST_DEVICE
void EquilibriumReactions::KernelWrapper::updateConcentrations( real64 const temperature,
                                                                arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                                                                arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                                                arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration ) const

{
  array2d<real64> matrix( m_numPrimarySpecies, m_numPrimarySpecies );
  array1d<real64> rhs( m_numPrimarySpecies );
  array1d<real64> solution( m_numPrimarySpecies );

  bool converged = false;
  for( int iteration = 0; iteration < m_maxNumIterations; iteration++ )
  {
    matrix.zero();
    rhs.zero();
    solution.zero();

    assembleEquilibriumReactionSystem( temperature,
                                       primarySpeciesTotalConcentration,
                                       primarySpeciesConcentration,
                                       secondarySpeciesConcentration,
                                       matrix,
                                       rhs );

    real64 const residualNorm = BlasLapackLA::vectorNorm2( rhs.toSliceConst() );

    if( residualNorm < m_newtonTol && iteration >= 1 )
    {
      converged = true;
      break;
    }

    #if defined(GEOSX_USE_CUDA)
    // for now we do this. We ll need a gpu solver for dense matrices.
    GEOSX_ERROR( "Geochem only works on CPUs for now." );
    #else
    BlasLapackLA::solveLinearSystem( matrix, rhs, solution );
    #endif

    updatePrimarySpeciesConcentrations( solution, primarySpeciesConcentration );
  }
}

GEOSX_HOST_DEVICE
void EquilibriumReactions::KernelWrapper::assembleEquilibriumReactionSystem( real64 const temperature,
                                                                             arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                                                                             arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                                                             arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                                                                             arrayView2d<real64> const matrix,
                                                                             arrayView1d<real64> const rhs ) const
{

  stackArray1d< real64, ReactionsBase::maxNumPrimarySpecies > log10PrimaryActCoeff( m_numPrimarySpecies );
  stackArray1d< real64, ReactionsBase::maxNumSecondarySpecies > log10SecActCoeff( m_numSecondarySpecies );
  stackArray2d< real64, ReactionsBase::maxNumSecondarySpecies * ReactionsBase::maxNumPrimarySpecies > dLog10SecConc_dLog10PrimaryConc( m_numPrimarySpecies, m_numSecondarySpecies );
  stackArray1d< real64, ReactionsBase::maxNumPrimarySpecies > totalConcentration( m_numPrimarySpecies );
  stackArray2d< real64, ReactionsBase::maxNumPrimarySpecies * ReactionsBase::maxNumPrimarySpecies > dTotalConc_dLog10PrimaryConc( m_numPrimarySpecies, m_numPrimarySpecies );

  real64 ionicStrength = 0.0;
  stackArray1d< real64, ReactionsBase::maxNumPrimarySpecies > dIonicStrength_dPrimaryConcentration( m_numPrimarySpecies );
  stackArray1d< real64, ReactionsBase::maxNumPrimarySpecies > dLog10PrimaryActCoeff_dIonicStrength( m_numPrimarySpecies );
  stackArray1d< real64, ReactionsBase::maxNumSecondarySpecies > dLog10SecActCoeff_dIonicStrength( m_numSecondarySpecies );

  /// activity coefficients
  computeIonicStrength( primarySpeciesConcentration,
                        secondarySpeciesConcentration,
                        ionicStrength );

  computeLog10ActCoefBDotModel( temperature,
                                ionicStrength,
                                log10PrimaryActCoeff,
                                dLog10PrimaryActCoeff_dIonicStrength,
                                log10SecActCoeff,
                                dLog10SecActCoeff_dIonicStrength );

  computeSeondarySpeciesConcAndDerivative( temperature,
                                           log10PrimaryActCoeff,
                                           dLog10PrimaryActCoeff_dIonicStrength,
                                           log10SecActCoeff,
                                           dLog10SecActCoeff_dIonicStrength,
                                           primarySpeciesConcentration,
                                           secondarySpeciesConcentration,
                                           dLog10SecConc_dLog10PrimaryConc );

  computeTotalConcAndDerivative( temperature,
                                 primarySpeciesConcentration,
                                 secondarySpeciesConcentration,
                                 dLog10SecConc_dLog10PrimaryConc,
                                 totalConcentration,
                                 dTotalConc_dLog10PrimaryConc );

  //Matteo: I assume we want to solve this to find the primary and the secondary species concentrations?
  for( int i=0; i<m_numPrimarySpecies; i++ )
  {
    rhs[i] = 1 - totalConcentration[i] / primarySpeciesTotalConcentration[i];
    for( int j=0; j<m_numPrimarySpecies; j++ )
    {
      matrix[i][j] = -dTotalConc_dLog10PrimaryConc[i][j] / primarySpeciesTotalConcentration[i];
    }
  }
}

// function to compute the derivative of the concentration of secondary species with respect to the concentration of the primary species.
GEOSX_HOST_DEVICE
void EquilibriumReactions::KernelWrapper::computeSeondarySpeciesConcAndDerivative( real64 const temperature,
                                                                                   arraySlice1d< real64 const > const & log10PrimaryActCoeff,
                                                                                   arraySlice1d< real64 const > const & dLog10PrimaryActCoeff_dIonicStrength,
                                                                                   arraySlice1d< real64 const > const & log10SecActCoeff,
                                                                                   arraySlice1d< real64 const > const & dLog10SecActCoeff_dIonicStrength,
                                                                                   arraySlice1d< real64 const > const & primarySpeciesConcentration,
                                                                                   arraySlice1d< real64 > const & secondarySpeciesConectration,
                                                                                   arraySlice2d< real64 > const & dLog10SecConc_dLog10PrimaryConc ) const
{
  // Compute d(concentration of dependent species)/d(concentration of basis species)
  for( int iSec = 0; iSec < m_numSecondarySpecies; ++iSec )
  {
    real64 log10SecConc = -m_log10EqConst[iSec] - log10SecActCoeff[iSec];

    for( int jPri = 0; jPri < m_numPrimarySpecies; ++jPri )
    {
      real64 const dIonicStrength_dPrimaryConc = log( 10 ) * 0.5 * m_chargePrimary[jPri] * m_chargePrimary[jPri];

      log10SecConc += m_stoichMatrix[iSec][jPri] * ( log10( primarySpeciesConcentration[jPri] ) + log10PrimaryActCoeff[jPri] );
      dLog10SecConc_dLog10PrimaryConc[iSec][jPri] += m_stoichMatrix[iSec][jPri] - dLog10SecActCoeff_dIonicStrength[iSec] * primarySpeciesConcentration[jPri] *
                                                     dIonicStrength_dPrimaryConc;
      for( int kDerivative = 0; kDerivative < m_numPrimarySpecies; ++kDerivative )
      {
        // add contribution to the derivtive from all primary activity coefficients
        dLog10SecConc_dLog10PrimaryConc[iSec][jPri] += m_stoichMatrix[iSec][kDerivative] * dLog10PrimaryActCoeff_dIonicStrength[kDerivative] * primarySpeciesConcentration[jPri] *
                                                       dIonicStrength_dPrimaryConc;
      }

    }
    secondarySpeciesConectration[iSec] = pow( 10, log10SecConc );
  }

}

GEOSX_HOST_DEVICE
void EquilibriumReactions::KernelWrapper::computeTotalConcAndDerivative( real64 const temperature,
                                                                         arraySlice1d< real64 const > const & primarySpeciesConcentration,
                                                                         arraySlice1d< real64 const > const & secondarySpeciesConectration,
                                                                         arraySlice2d< real64 const > const & dLog10SecConc_dLog10PrimaryConc,
                                                                         arraySlice1d< real64 > const & totalConc,
                                                                         arraySlice2d< real64 > const & dTotalConc_dLog10PrimaryConc ) const


{
  // This function computes the total concentration and its derivative with respect to log10(basis species concentrations).
  for( int iPri = 0; iPri < m_numPrimarySpecies; ++iPri )
  {
    totalConc[iPri] = primarySpeciesConcentration[iPri];
    // d(total concentration)/d(log10(concentration))
    dTotalConc_dLog10PrimaryConc[iPri][iPri] = log( 10.0 ) * primarySpeciesConcentration[iPri];
    // contribution from all dependent species
    for( int jSec = 0; jSec < m_numSecondarySpecies; ++jSec )
    {
      totalConc[iPri] += m_stoichMatrix[jSec][iPri] * secondarySpeciesConectration[jSec];
      for( int kDerivative = 0; kDerivative < m_numPrimarySpecies; ++kDerivative )
      {
        // add contribution to the derivtive from dependent species via the chain rule
        dTotalConc_dLog10PrimaryConc[iPri][kDerivative] += m_stoichMatrix[jSec][iPri] * log( 10.0 ) *
                                                           secondarySpeciesConectration[jSec] * dLog10SecConc_dLog10PrimaryConc[jSec][kDerivative];
      }
    }
  }
}

GEOSX_HOST_DEVICE
void EquilibriumReactions::KernelWrapper::
  updatePrimarySpeciesConcentrations( arrayView1d<real64 const> const solution,
                                      arraySlice1d< real64 > const & primarySpeciesConcentration ) const
{
  for( integer i = 0; i < m_numPrimarySpecies; i++ )
  {
    primarySpeciesConcentration[i] *= pow( 10, solution[i] );
  }
}

} // end namespace chemicalReactions

} // namespace constitutive

} // end namespace geosx


// // Not sure if this works
//   m_inputTotalConc[0:6] = 1E-20	// Does this notation work?
//   m_inputTotalConc[0] = pow( 10.0, -7 )
/*
   // Accurate value for H+ concentration. Hopefully it is not negative
   m_inputTotalConc[0] =
      0+m_inputTotalConc[1]-2*m_inputTotalConc[2]+2*m_inputTotalConc[3]+m_inputTotalConc[4]-2*m_inputTotalConc[5]-m_inputTotalConc[6];
 */

/*
   // Not sure if this is the right place to give the initial guess
   // Mismatch between lhs and rhs as one is log10 and the other isn't
   // Fixed it with the lines below
   m_log10PrimaryConc = log10(m_inputTotalConc);
   // If for some reason the total concentration of H+ is negative (not sure if this can happen)
   if (2*m_inputTotalConc[2]-2*m_inputTotalConc[3]-m_inputTotalConc[4]+2*m_inputTotalConc[5]+m_inputTotalConc[6]<0)
   {
    m_log10PrimaryConc[0] =
       log10(-2*m_inputTotalConc[2]+2*m_inputTotalConc[3]+m_inputTotalConc[4]-2*m_inputTotalConc[5]-m_inputTotalConc[6]);
   }
   else
   {
    m_log10PrimaryConc[0] = -7;
   }
 */
