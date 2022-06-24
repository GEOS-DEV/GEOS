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
 * @file EquilibriumReaction.cpp
 */

#include "constitutive/fluid/PVTFunctions/EquilibriumReaction.hpp"

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

EquilibriumReaction::EquilibriumReaction( string const & name ):
  ReactionBase( name )
{
  // Hardcoding values for now

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
  m_log10EqConst.resize( m_numSecSpecies );	// Not sure if this is the correct way of allocating the size
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

EquilibriumReaction::KernelWrapper EquilibriumReaction::createKernelWrapper() const
{
  return KernelWrapper(  );
}


GEOSX_HOST_DEVICE 
void EquilibriumReactionUpdate::updateConcentrations( real64 const & temperature,
                                                      arraySlice1d< real64 > const & totalConc,
                                                      arraySlice2d< real64 > const & dTotalConc_dLog10PrimaryConc )

{ 
  // This computes the value of the function we are trying to solve (mass balance equation)
  stackArray1d< ReactionBase::maxNumPrimarySpecies > log10PrimaryConc(m_numPrimarySpecies), log10PrimaryActCoeff(m_numPrimarySpecies);
  stackArray1d< ReactionBase::maxNumSecondarySpecies > log10SecConc(m_numSecSpecies);
  stackArray2d< ReactionBase::maxNumSecondarySpecies > dLog10SecConc_dLog10PrimaryConc(m_numPrimarySpecies, m_numSecSpecies);

  computeLog10SecConcAndDerivative( temperature, log10PrimaryConc, log10SecConc, dLog10SecConc_dLog10PrimaryConc );
  
  computeTotalConcAndDerivative( temperature, log10PrimaryConc, log10SecConc, dLog10SecConc_dLog10PrimaryConc, totalConc, dTotalConc_dLog10PrimaryConc );
  
  //Matteo: I assume we want to solve this to find the total concentration?
  stackArray1d<ReactionBase::maxNumPrimarySpecies> funValue (m_numPrimarySpecies);
  for(int i=0; i<funValue.size(); i++)
  {
    funValue[i] = totalConc[i] - inputTotalConc[i];
  }
}

// function to compute the derivative of the concentration of dependent species with respect to the concentration of the basis species.

GEOSX_HOST_DEVICE 
void EquilibriumReactionUpdate::computeLog10SecConcAndDerivative( real64 const temperature,
                                                                  arraySlice1d< real64 const > const & log10PrimaryConc,
                                                                  arraySlice1d< real64 > & log10SecConc,
                                                                  arraySlice2d< real64 > & dLog10SecConc_dLog10PrimaryConc )
{
  // Compute d(concentration of dependent species)/d(concentration of basis species)
  for( localIndex iSec = 0; iSec < m_numSecSpecies; ++iSec )
  {
    log10SecConc[iSec] = -m_log10EqConst[iSec] - m_log10SecActCoeff[iSec];

    for( localIndex jPri = 0; j < m_numPrimarySpecies; ++j )
    {
      log10SecConc[iSec] += m_stoichMatrix[iSec][jPri] * (log10PrimaryConc[jPri] + m_log10PrimaryActCoeff[jPri]);
      dLog10SecConc_dLog10PrimaryConc[iSec][jPri] += m_stoichMatrix[iSec][jPri];
    }
    
  }

}

GEOSX_HOST_DEVICE 
void EquilibriumReactionUpdate::computeTotalConcAndDerivative( real64 const & temperature,
                                                               arraySlice1d< real64 const > const & log10PrimaryConc,
                                                               arraySlice1d< real64 const > const & log10SecConc,
                                                               arraySlice2d< real64 const > const & dLog10SecConc_dLog10PrimaryConc,
                                                               arraySlice1d< real64 > const & totalConc,
                                                               arraySlice2d< real64 > const & dTotalConc_dLog10PrimaryConc)


{
  // This function computes the total concentration and its derivative with respect to log10(basis species concentrations). 
  for( localIndex iPri = 0; iPri < m_numPrimarySpecies; ++iPri )
  {
    real64 primaryConc = pow( 10.0, log10PrimaryConc[iPri] );		
    totalConc[iPri] = primaryConc;			
    dTotalConc_dLog10PrimaryConc[iPri][iPri] = log( 10.0 ) * primaryConc;		// d(total concentration)/d(log10(concentration))
    // contribution from all dependent species
    for( localIndex jSec = 0;  jSec < m_numSecSpecies; ++jSec )
    {
      real64 concSec = pow( 10.0, log10SecConc[jSec] );
      totalConc[iPri] += m_stoichMatrix[jSec][iPri] * concSec;	
      for( localIndex kDerivative = 0; kDerivative < m_numPrimarySpecies; ++kDerivative )		// add contribution to the derivtive from dependent species via the chain rule
      {
        dTotalConc_dLog10PrimaryConc[iPri][kDerivative] += m_stoichMatrix[jSec][iPri] * log( 10.0 ) * concSec * dLog10SecConc_dLog10PrimaryConc[jSec][kDerivative];
      }
    }
  }
}


/*
void EquilibriumReaction::allocateConstitutiveData( dataRepository::Group & const parent,
                                                    localIndex const numConstitutivePointsPerParentIndex )
{
  ReactionBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_databaseType = stringToDatabaseType( m_databaseTypeString );  // Have to define this datatype somewhere
  m_activityCoefModel = stringToActivityCoefModel( m_activityCoefModelString );   // Have to define this datatype somewhere

// Have to write this function somewhere
  ReadDatabase();

// Find whether H+ and H2O are specified as basis species
  bool HplusNotFound = 1;
  bool H2OFound = 0;
  localIndex const NBasis = numBasisSpecies();
  m_isHplus.resize( NBasis );
  for( localIndex id = 0; id < NBasis; ++id )
  {
    if( m_basisSpeciesNames[id] == "H+" )
    {
      HplusNotFound = 0;
      m_isHplus[id] = 1;
    }
    else
    {
      m_isHplus[id] = 0;
    }
    if( m_basisSpeciesNames[id] == "H2O" )
      H2OFound = 1;
  }

  GEOSX_ERROR_IF( HplusNotFound, "ReactionBase: H+ is not specified in basisSpeciesNames" );
  GEOSX_ERROR_IF( H2OFound, "ReactionBase: H2O cannot be specified in basisSpeciesNames" );

  resizeFields( parent.size() );
}

// initialize the concentration of the basis and dependent species and their derivatives.
// assignemnts associated with kinetic reaction rates have been commented out.
void EquilibriumReaction::resizeFields( localIndex const size )
{
  localIndex const NBasis = numBasisSpecies();    //These functions need to be defined somehwere
  localIndex const NDependent = numDependentSpecies();

  m_concentrationAct.resize( size, NBasis );    // I think concentration*ActivityCoefficient for the basis species
  m_dependentConc.resize( size, NDependent + 1 );
  m_dDependentConc_dConc.resize( size, NDependent + 1, NBasis + 1 );

  m_totalConc.resize( size, NBasis );     // I think this is the total concentration for mass balance. Not 100% sure
  m_dTotalConc_dConc.resize( size, NBasis, NBasis );
}

// Compute log10(ActivityCoefficient) for basis and dependent species along with their derivatives with respect of Ionic strength.
// Only B-Dot model has been impletemented.
// If there are plans to implement other models, maybe it is better to create multiple functions, one each for each activity coefficient
// model.
void EquilibriumReactionUpdate::computeLogActCoef( real64 const pressure,
                                                   real64 const temperature,
                                                   real64 const ionicStrength,
                                                   array1d< real64 > & logActCoef1,
                                                   array1d< real64 > & dLogActCoef1,
                                                   array1d< real64 > & logActCoef2,
                                                   array1d< real64 > & dLogActCoef2 )
{
  GEOSX_UNUSED_VAR( pressure );
  localIndex const NBasis = numBasisSpecies();
  localIndex const NDependent = numDependentSpecies();
  const array1d< Species > & dependentSpecies = m_thermoDatabase->GetDependentSpecies();  // Have to define the Species datatype somewhere
  const array1d< Species > & basisSpecies = m_thermoDatabase->GetBasisSpecies();
  const array1d< localIndex > & basisSpeciesIndices = m_thermoDatabase->GetBasisSpeciesIndices();
  const ActCoefParameters & actCoefParameters = m_thermoDatabase->GetActCoefParameters(); // Have to define the ActCoefParameters datatype
                                                                                          // somewhere

  logActCoef1.resize( NBasis );
  dLogActCoef1.resize( NBasis );
  logActCoef2.resize( NDependent );
  dLogActCoef2.resize( NDependent );

  real64 DHA, DHB, Bdot, DHazero, charge;
  interpolation( "DHA", actCoefParameters.temperatures, actCoefParameters.DHAs, temperature, DHA ); // Have to define this function
                                                                                                    // somehwere. Maybe there is already one
                                                                                                    // defined somewhere
  interpolation( "DHB", actCoefParameters.temperatures, actCoefParameters.DHBs, temperature, DHB );

  // constants for the BDot model
  static real64 C=-1.0312;
  static real64 F=0.0012806;
  static real64 G=255.9;
  static real64 E=0.4445;
  static real64 H=-0.001606;
  real64 TK = temperature + 273.15;

  if( m_activityCoefModel == ActivityCoefModel::BDot )    // Have to define ActivityCoefModel somewhere
  {
    interpolation( "Bdot", actCoefParameters.temperatures, actCoefParameters.BDots, temperature, Bdot );
    for( localIndex i = 0; i < NBasis; ++i )
    {
      localIndex ic = basisSpeciesIndices[i];
      DHazero = basisSpecies[ic].DHazero;
      charge = basisSpecies[ic].charge;
      logActCoef1[i] = Bdot * ionicStrength - DHA * charge * charge * sqrt( ionicStrength ) / (1.0 + DHazero * DHB * sqrt( ionicStrength ));
      dLogActCoef1[i] = Bdot - DHA * charge * charge *
                        (0.5 / sqrt( ionicStrength ) / (1.0 + DHazero * DHB * sqrt( ionicStrength )) - 0.5 * DHazero * DHB / (1.0 + DHazero * DHB * sqrt( ionicStrength )) /
                         (1.0 + DHazero * DHB * sqrt( ionicStrength ))); // dlog10(activityCoefficient)/dI
    }
    for( localIndex i = 0; i < dependentSpecies.size(); ++i )
    {
      DHazero = dependentSpecies[i].DHazero;
      charge = dependentSpecies[i].charge;
      // special consideration for dissolved CO2, O2 and H2
      if( dependentSpecies[i].name == "CO2(aq)" || dependentSpecies[i].name == "H2(aq)" || dependentSpecies[i].name == "O2(aq)" )
      {
        logActCoef2[i] = ((C + F * TK + G / TK) * ionicStrength -(E + H * TK) * (ionicStrength / (ionicStrength+1) )) / 2.303;
        dLogActCoef2[i] = ((C + F * TK + G / TK) - (E + H * TK) * (1.0 / (ionicStrength + 1) - ionicStrength / (ionicStrength + 1) / (ionicStrength + 1))) / 2.303;
      }
      else if( strstr( dependentSpecies[i].name.c_str(), "(aq)" ) || dependentSpecies[i].type == SpeciesType::Gas || dependentSpecies[i].type == SpeciesType::Solid )
      {
        logActCoef2[i] = 0; //Activity coefficients of gases, solid, and soluble gases other than CO2, O2, and H2 is 1 - Log10(1) = 0.
        dLogActCoef2[i] = 0;
      }
      else
      {
        logActCoef2[i] = Bdot * ionicStrength - DHA * charge * charge * sqrt( ionicStrength ) / (1.0 + DHazero * DHB * sqrt( ionicStrength ));
        dLogActCoef2[i] = Bdot - DHA * charge * charge *
                          (0.5 / sqrt( ionicStrength ) / (1.0 + DHazero * DHB * sqrt( ionicStrength )) - 0.5 * DHazero * DHB / (1.0 + DHazero * DHB * sqrt( ionicStrength )) /
                           (1.0 + DHazero * DHB * sqrt( ionicStrength )));
      }
    }
  }
  else
    GEOSX_ERROR( "wrong activity coef model" );
}

// function to compute the derivative of the concentration of dependent species with respect to the concentration of the basis species.
// also seems to update the concentration of the dependent species
// computes ionic strength using the input concentration of dependent species. In other words, ionic strength is computed using old
// dependent concentrations not the ones computed later in the function
// calls the function that computes acitivity coefficient and its derivatives.
// Perhaps a more descriptive function name is useful.
void EquilibriumReaction::Compute( real64 const pressure,
                                   real64 const temperature,
                                   arraySlice1d< real64 const > const & concentration,    // Looks like it is log10(concentration)
                                   arraySlice1d< real64 > const & dependentConc,    // Why are these listed as const when they are getting updated?
                                   arraySlice2d< real64 > const & dDependentConc_dConc,
                                   arraySlice1d< real64 > const & concentrationAct,
                                   ThermoDatabase & thermoDatabase )
{
  localIndex const NBasis = numBasisSpecies();
  localIndex const NDependent = numDependentSpecies();
  const array1d< Species > & dependentSpecies = thermoDatabase->GetDependentSpecies();
  const array1d< Species > & basisSpecies = m_thermoDatabase->GetBasisSpecies();
  const array1d< localIndex > & basisSpeciesIndices = m_thermoDatabase->GetBasisSpeciesIndices();
  const ActCoefParameters & actCoefParameters = m_thermoDatabase->GetActCoefParameters();

  static bool firstTime = 1;

  //get ionic strength
  dependentConc[NDependent] = 0;
  for( localIndex i = 0; i < NBasis; ++i )
  {
    localIndex id = basisSpeciesIndices[i];
    real64 charge = basisSpecies[id].charge;
    real64 conc = pow( 10.0, concentration[i] );
    dependentConc[NDependent] += 0.5 * charge * charge * conc;
  }
  for( localIndex i = 0; i < dependentSpecies.size(); ++i )
  {
    real64 charge = dependentSpecies[i].charge;
    real64 conc = pow( 10.0, dependentConc[i] );
    dependentConc[NDependent] += 0.5 * charge * charge * conc;
  }

  if( firstTime )
  {
    dependentConc[NDependent] = 1.0; // Initial guess, need to fix it later
    firstTime = 0;
  }

  //compute log(Activity coefficient) and d(log(Activity Coefficient))/dI
  array1d< real64 > logActCoef1;
  array1d< real64 > dLogActCoef1;
  array1d< real64 > logActCoef2;
  array1d< real64 > dLogActCoef2;
  ComputeLogActCoef( pressure, temperature, dependentConc[NDependent], logActCoef1, dLogActCoef1, logActCoef2, dLogActCoef2 );

// Compute d(concentration of dependent species)/d(concentration of basis species)
  real64 logK;
  for( localIndex i = 0; i < NDependent + 1; ++i )
    for( localIndex j = 0; j < NBasis + 1; ++j )
      dDependentConc_dConc[i][j] = 0;

  for( localIndex i = 0; i < dependentSpecies.size(); ++i )
  {
    // special consideation for solid species (derivative and concentration = 0?)
    if( dependentSpecies[i].type == SpeciesType::Solid )
    {
      dependentConc[i] = 0.0;
      continue;
    }
    real64 nu1 = dependentSpecies[i].stochs[0];   // Get the stoichiometric coefficient for the dependent species currently being reviewed.
    interpolation( "logK", actCoefParameters.temperatures, dependentSpecies[i].logKs, temperature, logK );

    dependentConc[i] = logK / nu1 - logActCoef2[i] * m_useActCoef;  // I think m_useActCoef is a flag to determine whether to use activity
                                                                    // coefficients or not.
    for( localIndex j = 1; j < dependentSpecies[i].speciesIndices.size(); ++j )
    {
      real64 nu2 = dependentSpecies[i].stochs[j];
      localIndex ic = dependentSpecies[i].speciesIndices[j];
      if( ic == NBasis )
      {
        dependentConc[i] -= nu2 / nu1 * m_logActH2O;  // Looks like the second last entry is meant for water
      }
      else if( ic == NBasis + 1 )
        dependentConc[i] -= nu2 / nu1 * m_logFO2g;  // Looksl like the last entry is meant for O2 gas
      else
      {
        dependentConc[i] -= nu2 / nu1 * (concentration[ic] + logActCoef1[ic] * m_useActCoef);
        dDependentConc_dConc[i][ic] -= nu2 / nu1;
      }
    }
  }

  // compute log10(activity) = log10(concentration) + log10(activity coefficient)
  for( localIndex i = 0; i < NBasis; ++i )
    concentrationAct[i] = concentration[i] + logActCoef1[i] * m_useActCoef;

}

// currently this function computes the total concentration and its derivative with respect to log10(basis species concentrations).
void EquilibriumReaction::ComputeChemistry( real64 const & temperature,
                                            arraySlice1d< real64 const > const & concentration,
                                            arraySlice1d< real64 const > const & concAct,
                                            arraySlice1d< real64 const > const & dependentConc,
                                            arraySlice2d< real64 const > const & dDependentConc_dConc,
                                            arraySlice1d< real64 > const & totalConc,
                                            arraySlice2d< real64 > const & dTotalConc_dConc )
{
  static const real64 RConst = 8.314;

  localIndex const NBasis = numBasisSpecies();
  localIndex const NDependent = numDependentSpecies();

  for( localIndex ic = 0; ic < NBasis; ++ic )
    for( localIndex jc = 0; jc < NBasis; ++jc )
      dTotalConc_dConc[ic][jc] = 0.0;

  for( localIndex ic = 0; ic < NBasis; ++ic )
  {
    real64 concBasis = pow( 10.0, concentration[ic] );    // the variable concentration here is defined as log10(concentration)
    totalConc[ic] = concBasis;      // for each basis species, total concentration = basis concentration + sum of the contribution of all
                                    // dependent species to the basis species. Not defined on a log10 basis
    dTotalConc_dConc[ic][ic] = log( 10.0 ) * concBasis;   // d(total concentration)/d(log10(concentration))
    // contribution from all dependent species
    for( localIndex id = 0; id < NDependent; ++id )
    {
      real64 concDependent = pow( 10.0, dependentConc[id] );
      totalConc[ic] -= m_stochMatrix[ic][id] * concDependent; // not entirely sure why the negative sign is introduced. m_stochMatrix can be
                                                              // defined such that the negative sign makes sense as long as we are
                                                              // consistent everywhere
      for( localIndex idc = 0; idc < NBasis; ++idc )    // add contribution to the derivtive from dependent species via the chain rule
      {
        dTotalConc_dConc[ic][idc] -= m_stochMatrix[ic][id] * log( 10.0 ) * concDependent * dDependentConc_dConc[id][idc];
      }
    }
  }

}
*/

REGISTER_CATALOG_ENTRY( ReactionBase, EquilibriumReaction, string const &, string_array const &, string_array const &, string_array const &, array1d< real64 > const & )

} // end namespace chemicalReactions

} // namespace constitutive

} // end namespace geosx
