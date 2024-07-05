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

namespace detail
{

// Lists the possible numbers of components supported by a fluid type
// Most of the fluid types support only 2 components so this is the default
template< typename FLUID_TYPE >
struct Components
{
  using type = camp::int_seq< integer, 2 >;
};

// Dead oil fluid model supports 2 or 3 components
template<>
struct Components< DeadOilFluid >
{
  using type = camp::int_seq< integer, 2, 3 >;
};

// Blackoil fluid model supports 3 components
template<>
struct Components< BlackOilFluid >
{
  using type = camp::int_seq< integer, 3 >;
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

/**
 * @brief Structure to execute a lambda on a fluid model depending on the number of components
 * @tparam THERMAL_ACTIVE Determines if the call has been deactivated due to the lambda being thermal
 *         but the fluid model not being thermal
 * @tparam COMPONENT_LIST Type listing the components to expand over
 */
template< bool THERMAL_ACTIVE, typename COMPONENT_LIST >
struct ComponentSelector
{};

/** @brief If the call is deactivated due to the fluid not being thermal
 */
template< typename COMPONENT_LIST >
struct ComponentSelector< false, COMPONENT_LIST >
{
  template< typename FluidType, typename LAMBDA >
  static void execute( int GEOS_UNUSED_PARAM( numComps ), FluidType & GEOS_UNUSED_PARAM( fluid ), LAMBDA && GEOS_UNUSED_PARAM( lambda ) )
  {}
};

template< integer ... Is >
struct ComponentSelector< true, camp::int_seq< integer, Is ... > >
{
  template< typename FluidType, typename LAMBDA >
  static void execute( int numComps, FluidType & fluid, LAMBDA && lambda )
  {
// With gcc8, the fold expression below issues a spurious warning
// warning: suggest parentheses around '&&' within '||' [-Wparentheses]
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=94505
#if (defined(__GNUC__) && (__GNUC__ < 10))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"
#endif
    bool const supported = ( ((numComps == Is) && (lambda( fluid, std::integral_constant< integer, Is >() ), true)) || ...) || false;
#if (defined(__GNUC__) && (__GNUC__ < 10))
#pragma GCC diagnostic pop
#endif
    if( !supported )
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
  constitutiveUpdatePassThru( fluidBase, [&]( auto & fluid )
  {
    using FluidType = TYPEOFREF( fluid );
    using Components = typename detail::Components< FluidType >::type;
    constexpr bool active = !THERMAL || FluidType::thermal();
    detail::ComponentSelector< active, Components >::execute( numComps, fluid, std::forward< LAMBDA >( lambda ));
  } );
}

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDSELECTOR_HPP_
