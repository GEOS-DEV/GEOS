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
{ // namespace

KineticReactions::KineticReactions( string const & name, integer const numPrimarySpecies, integer const numSecSpecies ):
  ReactionBase( name, numPrimarySpecies, numSecSpecies )
{
  // Here we should either read the database or the input values.
}

KineticReactions::KernelWrapper KineticReactions::createKernelWrapper() const
{
  return KernelWrapper( m_log10EqConst,
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
void KineticReactions::KernelWrapper::computeReactionRates( real64 const & temperature,
                                                            arraySlice1d< real64 const > const & concentration,
                                                            arraySlice1d< real64 const > const & surfaceArea0,
                                                            arraySlice1d< real64 const > const & volumeFraction0,
                                                            arraySlice1d< real64 const > const & volumeFraction,
                                                            real64 const & porosity0,
                                                            real64 const & porosity,
                                                            arraySlice1d< real64 > const & kineticReactionRate )
{
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

} // end namespace chemicalReactions

} // namespace constitutive

} // end namespace geosx
