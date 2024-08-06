/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HydrogenFluidUpdate.cpp
 */

#include "HydrogenFluidUpdate.hpp"

namespace geos
{

namespace constitutive
{

HydrogenFluidUpdate::HydrogenFluidUpdate( arrayView1d< geos::real64 const > componentMolarWeight,
                                          bool const useMass,
                                          bool const isThermal,
                                          PhaseProp::ViewType phaseFraction,
                                          PhaseProp::ViewType phaseDensity,
                                          PhaseProp::ViewType phaseMassDensity,
                                          PhaseProp::ViewType phaseViscosity,
                                          PhaseProp::ViewType phaseEnthalpy,
                                          PhaseProp::ViewType phaseInternalEnergy,
                                          PhaseComp::ViewType phaseCompFraction,
                                          FluidProp::ViewType totalDensity )
  : MultiFluidBase::KernelWrapper( std::move( componentMolarWeight ),
                                   useMass,
                                   std::move( phaseFraction ),
                                   std::move( phaseDensity ),
                                   std::move( phaseMassDensity ),
                                   std::move( phaseViscosity ),
                                   std::move( phaseEnthalpy ),
                                   std::move( phaseInternalEnergy ),
                                   std::move( phaseCompFraction ),
                                   std::move( totalDensity ) )
{
  GEOS_UNUSED_VAR( isThermal );
}

} //namespace constitutive

} //namespace geos
