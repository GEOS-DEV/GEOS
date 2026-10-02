/* SPDX-License-Identifier: LGPL-2.1-only */
#ifndef GEOS_HYDROSTATIC_MOBILITY_HPP
#define GEOS_HYDROSTATIC_MOBILITY_HPP
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include <array>
#include <functional>

namespace geos
{
namespace isothermalCompositionalMultiphaseBaseKernels
{
using HydrostaticPhaseVector = std::array< real64, 3 >;
using HydrostaticMobility = std::function< void( HydrostaticPhaseVector const &, HydrostaticPhaseVector & ) >;
/** Host evaluation of the actual non-hysteretic relative permeability law. */
HydrostaticMobility makeHydrostaticMobility( constitutive::RelativePermeabilityBase & model );
}
}
#endif
