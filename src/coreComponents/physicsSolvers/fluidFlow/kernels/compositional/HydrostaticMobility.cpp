/* SPDX-License-Identifier: LGPL-2.1-only */
#include "HydrostaticMobility.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/relativePermeability/RelativePermeabilitySelector.hpp"

namespace geos
{
namespace isothermalCompositionalMultiphaseBaseKernels
{
HydrostaticMobility makeHydrostaticMobility( constitutive::RelativePermeabilityBase & model )
{
  using namespace constitutive;
  integer const numPhases = model.numFluidPhases();
  GEOS_THROW_IF( numPhases < 2 || numPhases > 3,
                 "Hydrostatic capillary initialization requires two or three relative-permeability phases", InputError, model.getDataContext() );
  GEOS_THROW_IF( dynamicCast< TableRelativePermeabilityHysteresis * >( &model ) != nullptr,
                 "Hydrostatic capillary initialization does not support history-dependent relative permeability without an initialized-history coupling contract",
                 InputError, model.getDataContext() );
  model.forWrappers( []( dataRepository::WrapperBase const & wrapper ) { wrapper.move( hostMemorySpace, false ); } );
  HydrostaticMobility evaluate;
  constitutiveUpdatePassThru( model, [&]( auto & concrete )
  {
    if constexpr ( !std::is_same_v< TYPEOFREF( concrete ), TableRelativePermeabilityHysteresis > )
    {
      auto const wrapper = concrete.createKernelWrapper();
      evaluate = [wrapper, numPhases]( HydrostaticPhaseVector const & input, HydrostaticPhaseVector & output )
      {
        StackArray< real64, 2, 3, compflow::LAYOUT_PHASE > saturation( 1, numPhases );
        StackArray< real64, 3, 3, relperm::LAYOUT_RELPERM > trapped( 1, 1, numPhases ), values( 1, 1, numPhases );
        StackArray< real64, 4, 9, relperm::LAYOUT_RELPERM_DS > derivatives( 1, 1, numPhases, numPhases );
        for( integer ip = 0; ip < numPhases; ++ip ) saturation[0][ip] = input[ip];
        wrapper.compute( saturation[0].toSliceConst(), trapped[0][0], values[0][0], derivatives[0][0] );
        for( integer ip = 0; ip < numPhases; ++ip ) output[ip] = values[0][0][ip];
      };
    }
  } );
  return evaluate;
}
}
}
