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

struct CO2BrineFluidSelector
{};
struct CompositionalMultiphaseFluidSelector
{};

template< typename ... TYPES >
struct ConstitutivePassThruHandler< CompositionalMultiphaseFluidSelector, TYPES... >
{
  template< typename LAMBDA >
  static void execute( MultiFluidBase const & fluid, LAMBDA && lambda )
  {
    if( dynamicCast< CompositionalTwoPhasePengRobinsonConstantViscosity const * >( &fluid ) )
    {
      lambda( static_cast< CompositionalTwoPhasePengRobinsonConstantViscosity const & >( fluid ) );
    }
    else if( dynamicCast< CompositionalTwoPhaseSoaveRedlichKwongConstantViscosity const * >( &fluid ) )
    {
      lambda( static_cast< CompositionalTwoPhaseSoaveRedlichKwongConstantViscosity const & >( fluid ) );
    }
#ifdef GEOSX_USE_PVTPackage
    else if( dynamicCast< CompositionalMultiphaseFluidPVTPackage const * >( &fluid ) )
    {
      lambda( static_cast< CompositionalMultiphaseFluidPVTPackage const & >( fluid ) );
    }
#endif
    else
    {
      ConstitutivePassThruHandler< TYPES... >::execute( fluid, std::forward< LAMBDA >( lambda ) );
    }
  }

  template< typename LAMBDA >
  static void execute( MultiFluidBase & fluid, LAMBDA && lambda )
  {
    if( dynamicCast< CompositionalTwoPhasePengRobinsonConstantViscosity * >( &fluid ) )
    {
      lambda( static_cast< CompositionalTwoPhasePengRobinsonConstantViscosity & >( fluid ) );
    }
    else if( dynamicCast< CompositionalTwoPhaseSoaveRedlichKwongConstantViscosity * >( &fluid ) )
    {
      lambda( static_cast< CompositionalTwoPhaseSoaveRedlichKwongConstantViscosity & >( fluid ) );
    }
#ifdef GEOSX_USE_PVTPackage
    else if( dynamicCast< CompositionalMultiphaseFluidPVTPackage * >( &fluid ) )
    {
      lambda( static_cast< CompositionalMultiphaseFluidPVTPackage & >( fluid ) );
    }
#endif
    else
    {
      ConstitutivePassThruHandler< TYPES... >::execute( fluid, std::forward< LAMBDA >( lambda ) );
    }
  }
};

template< typename ... TYPES >
struct ConstitutivePassThruHandler< CO2BrineFluidSelector, TYPES... >
{
  template< typename LAMBDA >
  static void execute( MultiFluidBase const & fluid, LAMBDA && lambda )
  {
    if( dynamicCast< CO2BrinePhillipsFluid const * >( &fluid ) )
    {
      lambda( static_cast< CO2BrinePhillipsFluid const & >( fluid ) );
    }
    else if( dynamicCast< CO2BrineEzrokhiFluid const * >( &fluid ) )
    {
      lambda( static_cast< CO2BrineEzrokhiFluid const & >( fluid ) );
    }
    else if( dynamicCast< CO2BrinePhillipsThermalFluid const * >( &fluid ) )
    {
      lambda( static_cast< CO2BrinePhillipsThermalFluid const & >( fluid ) );
    }
    else if( dynamicCast< CO2BrineEzrokhiThermalFluid const * >( &fluid ) )
    {
      lambda( static_cast< CO2BrineEzrokhiThermalFluid const & >( fluid ) );
    }
    else
    {
      ConstitutivePassThruHandler< TYPES... >::execute( fluid, std::forward< LAMBDA >( lambda ) );
    }
  }

  template< typename LAMBDA >
  static void execute( MultiFluidBase & fluid, LAMBDA && lambda )
  {
    if( dynamicCast< CO2BrinePhillipsFluid * >( &fluid ) )
    {
      lambda( static_cast< CO2BrinePhillipsFluid & >( fluid ) );
    }
    else if( dynamicCast< CO2BrineEzrokhiFluid * >( &fluid ) )
    {
      lambda( static_cast< CO2BrineEzrokhiFluid & >( fluid ) );
    }
    else if( dynamicCast< CO2BrinePhillipsThermalFluid * >( &fluid ) )
    {
      lambda( static_cast< CO2BrinePhillipsThermalFluid & >( fluid ) );
    }
    else if( dynamicCast< CO2BrineEzrokhiThermalFluid * >( &fluid ) )
    {
      lambda( static_cast< CO2BrineEzrokhiThermalFluid & >( fluid ) );
    }
    else
    {
      ConstitutivePassThruHandler< TYPES... >::execute( fluid, std::forward< LAMBDA >( lambda ) );
    }
  }
};

template< typename LAMBDA >
void constitutiveUpdatePassThru( MultiFluidBase const & fluid,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< DeadOilFluid,
                               BlackOilFluid,
                               CompositionalMultiphaseFluidSelector,
                               CO2BrineFluidSelector
                               >::execute( fluid, std::forward< LAMBDA >( lambda ) );
}

template< typename LAMBDA >
void constitutiveUpdatePassThru( MultiFluidBase & fluid,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< DeadOilFluid,
                               BlackOilFluid,
                               CompositionalMultiphaseFluidSelector,
                               CO2BrineFluidSelector
                               >::execute( fluid, std::forward< LAMBDA >( lambda ) );
}

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDSELECTOR_HPP_
