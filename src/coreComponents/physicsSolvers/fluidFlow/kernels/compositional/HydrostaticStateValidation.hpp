/* SPDX-License-Identifier: LGPL-2.1-only */
#ifndef GEOS_HYDROSTATIC_STATE_VALIDATION_HPP
#define GEOS_HYDROSTATIC_STATE_VALIDATION_HPP
#include "HydrostaticMobility.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"
#include "mesh/ElementSubRegionBase.hpp"

namespace geos
{
namespace isothermalCompositionalMultiphaseBaseKernels
{
/** Temporary records are neither input, plot output, nor restart state. */
void recordHydrostaticState( ElementSubRegionBase & subRegion, localIndex element,
                             integer numPhases, HydrostaticPhaseVector const & phasePressure, real64 tolerance );
bool hasHydrostaticStateRecords( ElementSubRegionBase const & subRegion );
/** Reflash selected rows with normal solver code, preserving all unselected state. */
void realizeSelectedHydrostaticState( ElementSubRegionBase & subRegion, std::function< void() > const & update );
/** All ranks call this in the ordinary subregion initialization order. */
void validateRealizedHydrostaticState( ElementSubRegionBase & subRegion,
                                      constitutive::CapillaryPressureBase const & capillary );
}
}
#endif
