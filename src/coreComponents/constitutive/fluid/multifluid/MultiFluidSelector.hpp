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
 * @file multiFluidSelector.hpp
 */
#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDSELECTOR_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDSELECTOR_HPP_

#include "constitutive/ConstitutivePassThruHandler.hpp"
#include "constitutive/fluid/multifluid/blackOil/DeadOilFluid.hpp"
#include "constitutive/fluid/multifluid/blackOil/BlackOilFluid.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/CO2BrineFluid.hpp"

#include "common/GeosxConfig.hpp"
#ifdef GEOSX_USE_PVTPackage
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluidPVTPackage.hpp"
#endif
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluid.hpp"

namespace geos
{

namespace constitutive
{

template< typename LAMBDA >
void constitutiveUpdatePassThru( MultiFluidBase const & fluid,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< DeadOilFluid,
                               BlackOilFluid,
#ifdef GEOSX_USE_PVTPackage
                               CompositionalMultiphaseFluidPVTPackage,
#endif
                               CO2BrinePhillipsFluid,
                               CO2BrineEzrokhiFluid,
                               CO2BrinePhillipsThermalFluid,
// Including these in a CUDA build will lead to compiler segfault.
// Need to split compilation units for all the options
#if !defined(GEOS_DEVICE_COMPILE)
                               CO2BrineEzrokhiThermalFluid,
                               CompositionalTwoPhaseLohrenzBrayClarkViscosity,
#endif
                               CompositionalTwoPhaseConstantViscosity
                               >::execute( fluid, std::forward< LAMBDA >( lambda ) );
}

template< typename LAMBDA >
void constitutiveUpdatePassThru( MultiFluidBase & fluid,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< DeadOilFluid,
                               BlackOilFluid,
#ifdef GEOSX_USE_PVTPackage
                               CompositionalMultiphaseFluidPVTPackage,
#endif
                               CO2BrinePhillipsFluid,
                               CO2BrineEzrokhiFluid,
                               CO2BrinePhillipsThermalFluid,
// Including these in a CUDA build will lead to compiler segfault.
// Need to split compilation units for all the options"
#if !defined(GEOS_DEVICE_COMPILE)
                               CO2BrineEzrokhiThermalFluid,
                               CompositionalTwoPhaseLohrenzBrayClarkViscosity,
#endif
                               CompositionalTwoPhaseConstantViscosity
                               >::execute( fluid, std::forward< LAMBDA >( lambda ) );
}

namespace internal
{
// Lists the possible numbers of components supported by a fluid type
// Most of the fluid types support only 2 components so this is the default
template< typename FluidType >
struct Components
{
  using type = camp::int_seq< integer, 2 >;
};

// Blackoiub fluid models support anything from 2 or 3 components
template<>
struct Components< DeadOilFluid >
{
  using type = camp::int_seq< integer, 2, 3 >;
};
template<>
struct Components< BlackOilFluid >
{
  using type = camp::int_seq< integer, 2, 3 >;
};

// Compositional fluid models support anything from 2 to 5 components
template<>
struct Components< CompositionalMultiphaseFluidPVTPackage >
{
  using type = camp::int_seq< integer, 2, 3, 4, 5 >;
};
template<>
struct Components< CompositionalTwoPhaseLohrenzBrayClarkViscosity >
{
  using type = camp::int_seq< integer, 2, 3, 4, 5 >;
};
template<>
struct Components< CompositionalTwoPhaseConstantViscosity >
{
  using type = camp::int_seq< integer, 2, 3, 4, 5 >;
};

template< typename Components >
struct ComponentSelector
{};

template< integer ... Is >
struct ComponentSelector< camp::int_seq< integer, Is ... > >
{
  template< typename FluidType, typename LAMBDA >
  static void select( int numComps, FluidType & fluid, LAMBDA && lambda )
  {
    bool notSupported = false;
// With gcc8, the fold expression below issues a spurious warning
// warning: suggest parentheses around '&&' within '||' [-Wparentheses]
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=94505
#if (defined(__GNUC__) && (__GNUC__ < 10))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"
#endif
    ( ((numComps == Is) && (lambda( fluid, std::integral_constant< integer, Is >() ), true)) || ...) || (notSupported = true, false);
#if (defined(__GNUC__) && (__GNUC__ < 10))
#pragma GCC diagnostic pop
#endif
    if( notSupported )
    {
      GEOS_THROW( "Unsupported number of components: " << numComps << " for fluid " << FluidType::catalogName(), SimulationError );
    }
  }
};

}

template< bool THERMAL = false, typename LAMBDA = NoOpFunc >
void constitutiveComponentUpdatePassThru( MultiFluidBase & fluidBase,
                                          integer const numComps,
                                          LAMBDA && lambda )
{
  ConstitutivePassThruHandler< DeadOilFluid,
                               BlackOilFluid,
#ifdef GEOSX_USE_PVTPackage
                               CompositionalMultiphaseFluidPVTPackage,
#endif
                               CO2BrinePhillipsFluid,
                               CO2BrineEzrokhiFluid,
                               CO2BrinePhillipsThermalFluid,
                               CO2BrineEzrokhiThermalFluid,
                               CompositionalTwoPhaseLohrenzBrayClarkViscosity,
                               CompositionalTwoPhaseConstantViscosity
                               >::execute( fluidBase, [&]( auto & fluid )
  {
    using FluidType = TYPEOFREF( fluid );
    if constexpr (!THERMAL || FluidType::thermal())
    {
      using Components = typename internal::Components< FluidType >::type;
      internal::ComponentSelector< Components >::select( numComps, fluid, std::forward< LAMBDA >( lambda ));
    }
  } );
}

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDSELECTOR_HPP_
