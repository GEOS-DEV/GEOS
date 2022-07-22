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

#include "constitutive/fluid/chemicalReactions/KineticReactions.hpp"
#include "functions/FunctionManager.hpp"

namespace geosx
{

using namespace stringutilities;

namespace constitutive
{

namespace chemicalReactions
{

KineticReactions::KineticReactions( string const & name, integer const numPrimarySpecies, integer const numSecSpecies ):
  ReactionsBase( name, numPrimarySpecies, numSecSpecies )
{
  // Stochiometric Matrix for the kinetic reactions (in terms of primary species only)
  // First index: 0 = Ca(OH)2 dissolution, 1 = CaCO3 dissolution
  // Second index: 0 = H+, 1 = HCO3-, 2 = Ca+2, 3 = SO4-2, 4 = Cl-, 5 = Mg+2, 6 = Na+1
  m_stoichMatrix.resize( m_numKineticReaction, m_numPrimarySpecies );
  m_stoichMatrix[0][0] = -2;
  m_stoichMatrix[0][2] = 1;
  m_stoichMatrix[1][0] = -1;
  m_stoichMatrix[1][1] = 1;
  m_stoichMatrix[1][2] = 1;

  // Equilibrium Constant
  m_log10EqConst.resize( m_numKineticReactions );
  m_log10EqConst[0] = 20.19;
  m_log10EqConst[1] = 1.32;

  // Rate Constant
  // have to check the values as the functional form of the rate equation here differs from the implementation in GEOS
  m_reactionRateConstant.resize( m_numKineticReactions );
  m_reactionRateConstant[0] = 9.95e-1;
  m_reactionRateConstant[1] = 9.95e-3;
  // Here we should either read the database or the input values.
}

KineticReactions::KernelWrapper KineticReactions::createKernelWrapper() const
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
                        m_WATEQBDot,
                        m_reactionRateConstant );
}

// function to calculation the reaction rate. Includes impact of temperature, concentration, surface area, volume fraction and porosity
void KineticReactions::KernelWrapper::computeReactionsRate( real64 const & temperature,
                                                            arraySlice1d< geosx::real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                                            arraySlice1d< geosx::real64, compflow::USD_COMP - 1 > const & log10PrimaryActCoeff,
                                                            real64 const & specificSurfaceArea,
                                                            arraySlice1d< real64 const > & reactionRates ) const
{

  for( int iRxn = 0; iRxn < m_numKineticReactions; iRxn++ )
  {
    real64 saturationIndex = -m_log10EqConst[iRxn];

    for( int iPri = 0; iPri < m_numPrimarySpecies; ++iPri )
    {
      saturationIndex += m_stoichMatrix[iRxn][iPri] * log10( primarySpeciesConcentration[iPri] );
      saturationIndex += m_stoichMatrix[iRxn][iPri] * log10PrimaryActCoeff[iPri];
    }

    reactionRates[iRxn] = specificSurfaceArea * (1.0 - pow( 10, saturationIndex ) ) * m_reactionRateConstant[ir];
  }
}

} // end namespace chemicalReactions

} // namespace constitutive

} // end namespace geosx

/*
  for( localIndex ir = 0; ir < NReaction; ++ir )
  {

    for( localIndex ic = 0; ic < kineticReaction.stochs.size(); ++ic )
    {
      dSIndex[id] = kineticReaction.stochs[ic];

    }


    for( localIndex ic = 0; ic < NBasis; ++ic )
    {

      dKineticReactionRate_dConc[ir][ic] = S * kineticReaction.rateConst * rateTemp * pow( 10.0, SIndex ) * log( 10.0 ) * dSIndex[ic];

    }

  }


  for( localIndex ir = 0; ir < NReaction; ++ir )
  {

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
// arraySlice1d< real64 const > const & concentration,
// arraySlice1d< real64 const > const & surfaceArea0,
// arraySlice1d< real64 const > const & volumeFraction0,
// arraySlice1d< real64 const > const & volumeFraction,
// real64 const & porosity0,
// real64 const & porosity,
// arraySlice1d< real64 > const & kineticReactionRate )
// {
  // for( localIndex ir = 0; ir < kineticReactionArray.size(); ++ir )
  // {
  //   const KineticReaction & kineticReaction = kineticReactionArray[ir];
  //   const array1d< localIndex > & basisSpeciesIndices = kineticReaction.basisSpeciesIndices;
  //   // calculate saturation index
  //   real64 SIndex = -kineticReaction.logK;
  //   for( localIndex ic = 0; ic < kineticReaction.stochs.size(); ++ic )
  //   {
  //     SIndex += kineticReaction.stochs[ic] * concentration[basisSpeciesIndices[ic] ];     // Check that the input "concentration" is
  //                                                                                         // actually ln(activity
  // coefficient*concentration)
  //   }
  //   // surface area is assumed to scale with the volume fraction. Check whether this is the volume fraction of the mineral
  //   // dissolving/precipitating. Not sure why porosity is included.
  //   real64 S = surfaceArea0[ir] * pow( volumeFraction[ir] / volumeFraction0[ir], 2.0/3.0 ) * pow( porosity / porosity0, 2.0/3.0 );
  //   // computing the rate at the correct temperature. Looks like EQ36 database has it at 298.15 K
  //   real64 rateTemp = exp( -kineticReaction.E / RConst * (1.0 / (temperature + 273.15) - 1.0 / 298.15));
  //   real64 SS = (pow( 10.0, SIndex ) - 1.0);
  //   kineticReactionRate[ir] = S * kineticReaction.rateConst * rateTemp * SS;
//}


*/