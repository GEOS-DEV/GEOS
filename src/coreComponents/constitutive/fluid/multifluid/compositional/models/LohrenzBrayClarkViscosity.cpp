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
 * @file LohrenzBrayClarkViscosity.cpp
 */

#include "LohrenzBrayClarkViscosity.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

LohrenzBrayClarkViscosity::LohrenzBrayClarkViscosity( string const & name,
                                                      ComponentProperties const & componentProperties ):
  FunctionBase( name, componentProperties )
{}

LohrenzBrayClarkViscosity::KernelWrapper
LohrenzBrayClarkViscosity::createKernelWrapper() const
{
  return KernelWrapper( );
}

GEOS_HOST_DEVICE
void LohrenzBrayClarkViscosityUpdate::compute( ComponentProperties::KernelWrapper const & componentProperties,
                                               real64 const & pressure,
                                               real64 const & temperature,
                                               arraySlice1d< real64 const > const & phaseComposition,
                                               arraySlice2d< real64 const > const & dPhaseComposition,
                                               real64 const & density,
                                               arraySlice1d< real64 const > const & dDensity,
                                               real64 & viscosity,
                                               arraySlice1d< real64 > const & dViscosity,
                                               bool useMass ) const
{
  GEOS_UNUSED_VAR( pressure );   // No-direct pressure dependence (instead through density)

  GEOS_UNUSED_VAR( componentProperties );
  GEOS_UNUSED_VAR( temperature );
  GEOS_UNUSED_VAR( phaseComposition, dPhaseComposition );
  GEOS_UNUSED_VAR( density, dDensity );
  GEOS_UNUSED_VAR( viscosity, dViscosity );
  GEOS_UNUSED_VAR( useMass );

#ifdef HAHAHA
  integer constexpr maxNumComps = Mul
                                  integer const numComponents = componentProperties.getNumberOfComponents();

  // Estimate pure component properties at dilute-gas conditions (pressure near atmospheric) using
  // Stiel and Thodos [1961] correlation: https://doi.org/10.1002/aic.690070416
  // Dilute viscosity is solely temperature dependent
  // Units are converted so componentViscosity is in centipoise to match original reference

  stackArray1d< real64 > componentDiluteViscosity( numComponents );
  stackArray1d< real64 > dComponentDiluteViscosity_dTemperature( numComponents );

  computeComponentDiluteViscosity_StielThodos( numComponents,
                                               temperature,
                                               componentDiluteViscosity,
                                               dComponentDiluteViscosity_dTemperature );

  // Estimate phase viscosity (in cp) at dilute gas conditions using either
  // the Herning and Zipperer [1936], Wilke [1950], or Brokaw[1968] mixture rule.
  // Store in final phaseViscosity array. The classic LBC model uses Herning-Zipperer,
  // but the other two may be more accurate.

  //computePhaseDiluteViscosity_HerningZipperer( numComponents,
  //computePhaseDiluteViscosity_Wilke( numComponents,
  computePhaseDiluteViscosity_Brokaw( numComponents,
                                      numPhases,
                                      componentDiluteViscosity,
                                      dComponentDiluteViscosity_dTemperature,
                                      phaseComposition,
                                      phaseViscosity );


  // Estimate phase viscosity at given (P,T) conditions using LBC [1964] correlation.
  // This is an additional term added to the dilute gas estimate above.

  computePhaseViscosity_LohrenzBrayClark( numComponents,
                                          numPhases,
                                          phaseDensity,
                                          phaseComposition,
                                          phaseViscosity );


  // scale centipoise to pascal.seconds

  for( integer p=0; p<numPhases; ++p )
  {
    phaseViscosity.value[p] *= 1e-3;
    for( integer k=0; k<phaseViscosity.derivs[p].size(); ++k )
    {
      phaseViscosity.derivs[p][k] *= 1e-3;
    }
  }
#endif
}

} // end namespace compositional

} // namespace constitutive

} // end namespace geos
