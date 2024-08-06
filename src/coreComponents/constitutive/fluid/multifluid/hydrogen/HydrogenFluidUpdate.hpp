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
 * @file HydrogenFluidUpdate.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_HYDROGEN_HYDROGENFLUIDUPDATE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_HYDROGEN_HYDROGENFLUIDUPDATE_HPP_

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @brief Kernel wrapper class for HydrogenFluid.
 */
class HydrogenFluidUpdate final : public MultiFluidBase::KernelWrapper
{
public:
  using PhaseProp = MultiFluidBase::PhaseProp;
  using PhaseComp = MultiFluidBase::PhaseComp;
  using FluidProp = MultiFluidBase::FluidProp;

public:
  HydrogenFluidUpdate( arrayView1d< real64 const > componentMolarWeight,
                       bool const useMass,
                       bool const isThermal,
                       PhaseProp::ViewType phaseFraction,
                       PhaseProp::ViewType phaseDensity,
                       PhaseProp::ViewType phaseMassDensity,
                       PhaseProp::ViewType phaseViscosity,
                       PhaseProp::ViewType phaseEnthalpy,
                       PhaseProp::ViewType phaseInternalEnergy,
                       PhaseComp::ViewType phaseCompFraction,
                       FluidProp::ViewType totalDensity );

  GEOS_HOST_DEVICE
  virtual void compute( real64 const pressure,
                        real64 const temperature,
                        arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                        PhaseProp::SliceType const phaseFraction,
                        PhaseProp::SliceType const phaseDensity,
                        PhaseProp::SliceType const phaseMassDensity,
                        PhaseProp::SliceType const phaseViscosity,
                        PhaseProp::SliceType const phaseEnthalpy,
                        PhaseProp::SliceType const phaseInternalEnergy,
                        PhaseComp::SliceType const phaseCompFraction,
                        FluidProp::SliceType const totalDensity ) const override;

  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure,
                       real64 const temperature,
                       arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const override;
};

GEOS_HOST_DEVICE
void HydrogenFluidUpdate::compute( real64 const pressure,
                                   real64 const temperature,
                                   arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                                   PhaseProp::SliceType const phaseFraction,
                                   PhaseProp::SliceType const phaseDensity,
                                   PhaseProp::SliceType const phaseMassDensity,
                                   PhaseProp::SliceType const phaseViscosity,
                                   PhaseProp::SliceType const phaseEnthalpy,
                                   PhaseProp::SliceType const phaseInternalEnergy,
                                   PhaseComp::SliceType const phaseCompFraction,
                                   FluidProp::SliceType const totalDensity ) const
{
  GEOS_UNUSED_VAR( pressure );
  GEOS_UNUSED_VAR( temperature );
  GEOS_UNUSED_VAR( composition );
  GEOS_UNUSED_VAR( phaseFraction );
  GEOS_UNUSED_VAR( phaseDensity );
  GEOS_UNUSED_VAR( phaseMassDensity );
  GEOS_UNUSED_VAR( phaseViscosity );
  GEOS_UNUSED_VAR( phaseEnthalpy );
  GEOS_UNUSED_VAR( phaseInternalEnergy );
  GEOS_UNUSED_VAR( phaseCompFraction );
  GEOS_UNUSED_VAR( totalDensity );
}

GEOS_HOST_DEVICE
void HydrogenFluidUpdate::update( localIndex const k,
                                  localIndex const q,
                                  real64 const pressure,
                                  real64 const temperature,
                                  arraySlice1d< geos::real64 const, compflow::USD_COMP - 1 > const & composition ) const
{
  compute( pressure,
           temperature,
           composition,
           m_phaseFraction( k, q ),
           m_phaseDensity( k, q ),
           m_phaseMassDensity( k, q ),
           m_phaseViscosity( k, q ),
           m_phaseEnthalpy( k, q ),
           m_phaseInternalEnergy( k, q ),
           m_phaseCompFraction( k, q ),
           m_totalDensity( k, q ) );
}

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_HYDROGEN_HYDROGENFLUIDUPDATE_HPP_
