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

#include "constitutive/fluid/PVTFunctions/EquilibriumReactions.hpp"

#include "constitutive/fluid/PVTFunctions/CO2EOSSolver.hpp"
#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"
#include "functions/FunctionManager.hpp"

namespace geosx
{

using namespace stringutilities;

namespace constitutive
{

namespace chemicalReactions
{

EquilibriumReactions::EquilibriumReactions( string const & name ):
  ReactionBase( name )
{
  // Hardcoding values for now

  // Equilibrium constants
  m_log10EqConst.resize( m_numSecSpecies );	
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

  // Activity coefficient related constants
  real64 m_DebyeHuckelA, m_DebyeHuckelB, m_WATEQBDot 
  m_DebyeHuckelA = 0.5465;
  m_DebyeHuckelB = 0.3346;
  m_WATEQBDot = 0.0438;
  arrayView1d<real64>  m_ionSizePrimary, m_ionSizeSec;  
  m_ionSizePrimary.resize( m_numPrimarySepcies )
  m_ionSizePrimary[0] = 9.00;
  m_ionSizePrimary[1] = 4.00;
  m_ionSizePrimary[2] = 6.00;
  m_ionSizePrimary[3] = 4.00;
  m_ionSizePrimary[4] = 3.00;
  m_ionSizePrimary[5] = 8.00;
  m_ionSizePrimary[6] = 4.00;

  m_ionSizeSec.resize( m_numSecSepcies )
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

  arrayView1d<real64>  m_chargePrimary, m_chargeSec;  // should be an integer and not real
  m_chargePrimary.resize( m_numPrimarySepcies )
  m_chargePrimary[0] = 1;
  m_chargePrimary[1] = -1;
  m_chargePrimary[2] = 2;
  m_chargePrimary[3] = -2;
  m_chargePrimary[4] = -1;
  m_chargePrimary[5] = 2;
  m_chargePrimary[6] = 1;

  m_chargeSec.resize( m_numSecSepcies )
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


  m_log10SecActCoeff.resize( m_numSecSpecies );	
  m_log10SecActCoeff[0:10] = 0	//Assume dilute solution for first pass

  m_log10PrimaryActCoeff.resize( m_numPrimarySpecies );	
  m_log10PrimaryActCoeff[0:6] = 0	//Assume dilute solution for first pass

  m_inputTotalConc.resize( m_numPrimarySpecies );	
  m_inputTotalConc[0:6] = 1E-20	
  // Accurate value for H+ concentration. Hopefully it is not negative. Have to check about the non-negative requirement for this
  m_inputTotalConc[0] = 0+m_inputTotalConc[1]-2*m_inputTotalConc[2]+2*m_inputTotalConc[3]+m_inputTotalConc[4]-2*m_inputTotalConc[5]-m_inputTotalConc[6];

  arrayView1d<real64>  m_log10PrimaryConc;	
  m_log10PrimaryConc.resize( m_numPrimarySpecies );	

  // Not sure if this is the right place to give the initial guess
  m_log10PrimaryConc = log10(m_inputTotalConc);
  // If for some reason the total concentration of H+ is negative (not sure if this can happen)
  if (2*m_inputTotalConc[2]-2*m_inputTotalConc[3]-m_inputTotalConc[4]+2*m_inputTotalConc[5]+m_inputTotalConc[6]<0)
  {
    m_log10PrimaryConc[0] = log10(-2*m_inputTotalConc[2]+2*m_inputTotalConc[3]+m_inputTotalConc[4]-2*m_inputTotalConc[5]-m_inputTotalConc[6]);
  }
  else
  {
    m_log10PrimaryConc[0] = -7;
  }

  arrayView1d<real64>  m_log10SecConc;	
  m_log10SecConc.resize( m_numSecSpecies );
  m_log10SecConc[0:10] = -20;

  arrayView2d<real64>  dLog10SecConc_dLog10PrimaryConc;	
  dLog10SecConc_dLog10PrimaryConc.resize( m_numSecSpecies, m_numPrimarySpecies );	
  dLog10SecConc_dLog10PrimaryConc[0:10][0:6] = 0;

  arrayView2d<real64>  dTotalConc_dLog10PrimaryConc;	
  dTotalConc_dLog10PrimaryConc.resize( m_numPrimarySpecies, m_numPrimarySpecies );
  dTotalConc_dLog10PrimaryConc[0:6][0:6] = 0;


  // Stochiometric Matrix
  // First index: 0 = OH-, 1 = CO2, 2 = CO3-2, 3 = H2CO3, 4 = CaHCO3+, 5 = CaCO3, 6 = CaSO4, 7 = CaCl+, 8 = CaCl2, 9 = MgSO4, 10 = NaSO4-
  // Second index: 0 = H+, 1 = HCO3-, 2 = Ca+2, 3 = SO4-2, 4 = Cl-, 5 = Mg+2, 6 = Na+1
  m_stoichMatrix.resize( m_numSecSpecies, m_numPrimarySpecies );
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
  m_log10EqConst.resize( m_numSecSpecies );	
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
  return KernelWrapper(  );
}


GEOSX_HOST_DEVICE 
void EquilibriumReactions::KernelWrapper::updateConcentrations( real64 const & temperature,
                                                                arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration
                                                                arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                                                arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration ) const

{ 
  DenseMatrix matrix(m_numPrimarySpecies, m_numPrimarySpecies);
  DenseVector rhs(m_numPrimarySpecies), solution(m_numPrimarySpecies);

  bool converged = false;

  for( int iteration = 0; iteration < m_maxNumIterations; iteration++ )
  {
    matrix.zero();
    rhs.zero();

    assembleEquilibriumReactionSystem( temperature,
                                       primarySpeciesTotalConcentration,
                                       primarySpeciesConcentration,
                                       secondarySpeciesConcentration );

    real64 const residualNorm = LvArray::tensorOps::l2norm(rhs); // may need the size as constexpr

    if( residualNorm < newtonTol && iteration >= 1 )
    {
      converged = true;
      break;
    }

    m_denseLinearSolver.solve( matrix, rhs, solution );

    updateSolution( solution );
 }
}

GEOSX_HOST_DEVICE 
void EquilibriumReactions::KernelWrapper::assembleEquilibriumReactionSystem( real64 const & temperature,
                                                                             arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration
                                                                             arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                                                             arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                                                                             DenseMatrix & matrix,
                                                                             DenseVector & rhs ) const
{

  stackArray1d< ReactionBase::maxNumPrimarySpecies > log10PrimaryConc(m_numPrimarySpecies);
  stackArray1d< ReactionBase::maxNumPrimarySpecies > log10PrimaryActCoeff(m_numPrimarySpecies);
  stackArray1d< ReactionBase::maxNumSecondarySpecies > log10SecConc(m_numSecSpecies);
  stackArray1d< ReactionBase::maxNumSecondarySpecies > log10SecActCoeff(m_numSecSpecies);

  stackArray2d< ReactionBase::maxNumSecondarySpecies > dLog10SecConc_dLog10PrimaryConc(m_numPrimarySpecies, m_numSecSpecies);
  stackArray1d<ReactionBase::maxNumPrimarySpecies> totalConcentration( m_numPrimarySpecies );

  // activity coefficients
  ComputeLog10ActCoefBDotModel( temperature,
                                ionicStrength,
                                log10PrimaryActCoeff,
                                dLog10PrimaryActCoeff_dIonicStrength,
                                log10SecActCoeff,
                                dLog10SecActCoeff_dIonicStrength);

  computeLog10SecConcAndDerivative( temperature, 
                                    log10PrimaryConc, 
                                    log10SecConc, 
                                    dLog10SecConc_dLog10PrimaryConc );
  
  computeTotalConcAndDerivative( temperature, 
                                 log10PrimaryConc, 
                                 log10SecConc, 
                                 dLog10SecConc_dLog10PrimaryConc, 
                                 totalConcentration, 
                                 dTotalConc_dLog10PrimaryConc );

  //Matteo: I assume we want to solve this to find the primary and the secondary species concentrations?
  for(int i=0; i<m_numPrimarySpecies; i++)
  {
    rhs[i] = totalConcentration[i] - primarySpeciesTotalConcentration[i];
    for ( int j=0; j<m_numPrimarySpecies; j++)
    {
      matrix[i][j] = dTotalConc_dLog10PrimaryConc[i][j];
    }
  }
}

// function to compute the derivative of the concentration of dependent species with respect to the concentration of the basis species.
GEOSX_HOST_DEVICE 
void EquilibriumReactions::KernelWrapper::computeLog10SecConcAndDerivative( real64 const temperature,
                                                                            arraySlice1d< real64 const > const & log10PrimaryConc,
                                                                            arraySlice1d< real64 > & log10SecConc,
                                                                            arraySlice2d< real64 > & dLog10SecConc_dLog10PrimaryConc ) const 
{
  // Compute d(concentration of dependent species)/d(concentration of basis species)
  for( int iSec = 0; iSec < m_numSecSpecies; ++iSec )
  {
    log10SecConc[iSec] = -m_log10EqConst[iSec] - m_log10SecActCoeff[iSec];

    for( int jPri = 0; j < m_numPrimarySpecies; ++j )
    {
      log10SecConc[iSec] += m_stoichMatrix[iSec][jPri] * (log10PrimaryConc[jPri] + m_log10PrimaryActCoeff[jPri]);
      dLog10SecConc_dLog10PrimaryConc[iSec][jPri] += m_stoichMatrix[iSec][jPri];
    }
  }

}

GEOSX_HOST_DEVICE 
void EquilibriumReactions::KernelWrapper::computeTotalConcAndDerivative( real64 const & temperature,
                                                                         arraySlice1d< real64 const > const & log10PrimaryConc,
                                                                         arraySlice1d< real64 const > const & log10SecConc,
                                                                         arraySlice2d< real64 const > const & dLog10SecConc_dLog10PrimaryConc,
                                                                         arraySlice1d< real64 > const & totalConc,
                                                                         arraySlice2d< real64 > const & dTotalConc_dLog10PrimaryConc)


{
  // This function computes the total concentration and its derivative with respect to log10(basis species concentrations). 
  for( int iPri = 0; iPri < m_numPrimarySpecies; ++iPri )
  {
    real64 const primaryConc = pow( 10.0, log10PrimaryConc[iPri] );		
    totalConc[iPri] = primaryConc;			
    dTotalConc_dLog10PrimaryConc[iPri][iPri] = log( 10.0 ) * primaryConc;		// d(total concentration)/d(log10(concentration))
    // contribution from all dependent species
    for( int jSec = 0;  jSec < m_numSecSpecies; ++jSec )
    {
      real64 const concSec = pow( 10.0, log10SecConc[jSec] );
      totalConc[iPri] += m_stoichMatrix[jSec][iPri] * concSec;	
      for( int kDerivative = 0; kDerivative < m_numPrimarySpecies; ++kDerivative )		// add contribution to the derivtive from dependent species via the chain rule
      {
        dTotalConc_dLog10PrimaryConc[iPri][kDerivative] += m_stoichMatrix[jSec][iPri] * log( 10.0 ) * concSec * dLog10SecConc_dLog10PrimaryConc[jSec][kDerivative];
      }
    }
  }
}

REGISTER_CATALOG_ENTRY( ReactionBase, EquilibriumReactions, string const &, string_array const &, string_array const &, string_array const &, array1d< real64 > const & )

} // end namespace chemicalReactions

} // namespace constitutive

} // end namespace geosx
