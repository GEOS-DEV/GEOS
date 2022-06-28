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
 * @file KineticReactions.cpp
 */

#include "constitutive/fluid/PVTFunctions/KineticReactions.hpp"

#include "constitutive/fluid/PVTFunctions/CO2EOSSolver.hpp"
#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"
#include "functions/FunctionManager.hpp"

namespace geosx
{

using namespace stringutilities;

namespace constitutive
{

namespace chemicalReactions
{} // namespace

KineticReactions::KineticReactions( string const & name,
                                          string_array const & inputParams,
                                          string_array const & phaseNames,
                                          string_array const & componentNames,
                                          array1d< real64 > const & componentMolarWeight ):
  ReactionBase( name,
                componentNames,
                componentMolarWeight )
{}

KineticReactions::KernelWrapper KineticReactions::createKernelWrapper() const
{
  return KernelWrapper(  );
}

// initialize the concentration of the basis and dependent species and their derivatives.
void KineticReactions::resizeFields( localIndex const size )
{
  localIndex const NReaction = numKineticReaction();

  m_concentrationAct.resize( size, NBasis );    // I think concentration*ActivityCoefficient for the basis species
  m_dependentConc.resize( size, NDependent + 1 );
  m_dDependentConc_dConc.resize( size, NDependent + 1, NBasis + 1 );

  m_kineticReactionRate.resize( size, NReaction );	
  m_dKineticReactionRate_dConc.resize( size, NReaction, NBasis );

  m_kineticSpeciesReactionRate.resize( size, NBasis );
  m_dKineticSpeciesReactionRate_dConc.resize( size, NBasis, NBasis );
}

// function to compute the derivative of the concentration of dependent species with respect to the concentration of the basis species.
// also seems to update the concentration of the dependent species
// computes ionic strength using the input concentration of dependent species. In other words, ionic strength is computed using old
// dependent concentrations not the ones computed later in the function
// calls the function that computes acitivity coefficient and its derivatives.
// Perhaps a more descriptive function name is useful.
void KineticReactions::Compute( real64 const pressure,
                                   real64 const temperature,
                                   arraySlice1d< real64 const > const & concentration,    // Looks like it is log10(concentration)
                                   arraySlice1d< real64 > const & dependentConc,    // Why are these listed as const when they are getting
                                                                                    // updated?
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
// for each basis species, total concentration = basis concentration + sum of the contribution of all dependent species to the basis
// species.
// The kinetic reactions aspect of the code has been commented out
void KineticReactions::ComputeChemistry( real64 const & temperature,
                                            arraySlice1d< real64 const > const & concentration,
                                            arraySlice1d< real64 const > const & surfaceArea0,		// All the commented out enteries are
                                            arraySlice1d< real64 const > const & volumeFraction0,
                                            arraySlice1d< real64 const > const & volumeFraction,
                                            real64 const & porosity0,
                                            real64 const & porosity,
                                            arraySlice1d< real64 const > const & concAct,
                                            arraySlice1d< real64 const > const & dependentConc,
                                            arraySlice2d< real64 const > const & dDependentConc_dConc,
                                            arraySlice1d< real64 > const & totalConc,
                                            arraySlice2d< real64 > const & dTotalConc_dConc )
                                           arraySlice1d< real64 > const & kineticReactionRate,
                                           arraySlice2d< real64 > const & dKineticReactionRate_dConc, 
                                           arraySlice1d< real64 > const & kineticSpeciesReactionRate,
                                           arraySlice2d< real64 > const & dKineticSpeciesReactionRate_dConc,
                                           KineticReactions & kineticReactions )

{
  static const real64 RConst = 8.314;

  localIndex const NBasis = numBasisSpecies();
  localIndex const NDependent = numDependentSpecies();
//  localIndex const NReaction = numKineticReaction();

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
   const array1d< KineticReaction > & kineticReactionArray = kineticReactions->GetKineticReactions();

   for( localIndex ir = 0; ir < NReaction; ++ir )
   {

    const KineticReaction & kineticReaction = kineticReactionArray[ir];
    const array1d< localIndex > & basisSpeciesIndices = kineticReaction.basisSpeciesIndices;

    real64 SIndex = -kineticReaction.logK;

    array1d< real64 > dSIndex( NBasis );
    for( localIndex ic = 0; ic < NBasis; ++ic )
      dSIndex[ic] = 0.0;

    for( localIndex ic = 0; ic < kineticReaction.stochs.size(); ++ic )
    {
      localIndex id = basisSpeciesIndices[ic];
      SIndex += kineticReaction.stochs[ic] * concAct[id];
      dSIndex[id] = kineticReaction.stochs[ic];

    }

    real64 S = surfaceArea0[ir] * pow( volumeFraction[ir] / volumeFraction0[ir], 2.0/3.0 ) * pow( porosity / porosity0, 2.0/3.0 );

    real64 rateTemp = exp( -kineticReaction.E / RConst * (1.0 / (temperature + 273.15) - 1.0 / 298.15));

    real64 SS = (pow( 10.0, SIndex ) - 1.0);

    kineticReactionRate[ir] = S * kineticReaction.rateConst * rateTemp * SS;

    for( localIndex ic = 0; ic < NBasis; ++ic )
    {

      dKineticReactionRate_dConc[ir][ic] = S * kineticReaction.rateConst * rateTemp * pow( 10.0, SIndex ) * log( 10.0 ) * dSIndex[ic];

    }

   }

   for( localIndex ic = 0; ic < NBasis; ++ic )
   {
    kineticSpeciesReactionRate[ic] = 0.0;

    for( localIndex id = 0; id < NBasis; ++id )
      dKineticSpeciesReactionRate_dConc[ic][id] = 0.0;
   }

   for( localIndex ir = 0; ir < NReaction; ++ir )
   {

    const KineticReaction & kineticReaction = kineticReactionArray[ir];
    const array1d< localIndex > & basisSpeciesIndices = kineticReaction.basisSpeciesIndices;

    for( localIndex i = 0; i < (kineticReaction.stochs).size(); ++i )
    {

      localIndex ic = basisSpeciesIndices[i];

      kineticSpeciesReactionRate[ic] += -(kineticReaction.stochs)[i] * kineticReactionRate[ir];

      for( localIndex id = 0; id < NBasis; ++id )
      {

        dKineticSpeciesReactionRate_dConc[ic][id] += -(kineticReaction.stochs)[i] * dKineticReactionRate_dConc[ir][id];

      }

    }

   }
}

// function to calculation the reaction rate. Includes impact of temperature, concentration, surface area, volume fraction and porosity
void KineticReactionUpdate::ComputeKineticReactionRate( real64 const & temperature,
                                                      arraySlice1d< real64 const > const & concentration,
                                                      arraySlice1d< real64 const > const & surfaceArea0,
                                                      arraySlice1d< real64 const > const & volumeFraction0,
                                                      arraySlice1d< real64 const > const & volumeFraction,
                                                      real64 const & porosity0,
                                                      real64 const & porosity,
                                                      arraySlice1d< real64 > const & kineticReactionRate,
                                                      KineticReactions & kineticReactions )
{
  static const real64 RConst = 8.314;
  // Get reaction rate information
  const array1d< KineticReaction > & kineticReactionArray = kineticReactions->GetKineticReactions();  // Have to read this from input
                                                                                                      // file/database/somewhere before
                                                                                                      // using here
  for( localIndex ir = 0; ir < kineticReactionArray.size(); ++ir )
  {
    const KineticReaction & kineticReaction = kineticReactionArray[ir];
    const array1d< localIndex > & basisSpeciesIndices = kineticReaction.basisSpeciesIndices;
    // calculate saturation index
    real64 SIndex = -kineticReaction.logK;
    for( localIndex ic = 0; ic < kineticReaction.stochs.size(); ++ic )
    {
      SIndex += kineticReaction.stochs[ic] * concentration[basisSpeciesIndices[ic] ];     // Check that the input "concentration" is
                                                                                          // actually ln(activity coefficient*concentration)
    }
    // surface area is assumed to scale with the volume fraction. Check whether this is the volume fraction of the mineral
    // dissolving/precipitating. Not sure why porosity is included.
    real64 S = surfaceArea0[ir] * pow( volumeFraction[ir] / volumeFraction0[ir], 2.0/3.0 ) * pow( porosity / porosity0, 2.0/3.0 );
    // computing the rate at the correct temperature. Looks like EQ36 database has it at 298.15 K
    real64 rateTemp = exp( -kineticReaction.E / RConst * (1.0 / (temperature + 273.15) - 1.0 / 298.15));
    real64 SS = (pow( 10.0, SIndex ) - 1.0);
    kineticReactionRate[ir] = S * kineticReaction.rateConst * rateTemp * SS;
  }
}



REGISTER_CATALOG_ENTRY( ReactionBase, KineticReactions, string const &, string_array const &, string_array const &, string_array const &, array1d< real64 > const & )

} // end namespace chemicalReactions

} // namespace constitutive

} // end namespace geosx
