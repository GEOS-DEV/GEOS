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
#include "constitutive/fluid/DeadOilFluid.hpp"
#include "constitutive/fluid/BlackOilFluid.hpp"
#include "constitutive/fluid/CO2BrineFluid.hpp"

#include "common/GeosxConfig.hpp"
#ifdef GEOSX_USE_PVTPackage
#include "constitutive/fluid/CompositionalMultiphaseFluid.hpp"
#endif

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
                               CompositionalMultiphaseFluid,
#endif
                               CO2BrinePhillipsFluid,
                               CO2BrineEzrokhiFluid,
                               CO2BrinePhillipsThermalFluid
//                               ,CO2BrineEzrokhiThermalFluid   "Uncommenting this will lead to compiler segfault. Need to split compilation
// units for all the options"
                               >::execute( fluid, std::forward< LAMBDA >( lambda ) );
}

template< typename LAMBDA >
void constitutiveUpdatePassThru( MultiFluidBase & fluid,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< DeadOilFluid,
                               BlackOilFluid,
#ifdef GEOSX_USE_PVTPackage
                               CompositionalMultiphaseFluid,
#endif
                               CO2BrinePhillipsFluid,
                               CO2BrineEzrokhiFluid,
                               CO2BrinePhillipsThermalFluid
                               //,CO2BrineEzrokhiThermalFluid   "Uncommenting this will lead to compiler segfault. Need to split compilation
// units for all the options"
                               >::execute( fluid, std::forward< LAMBDA >( lambda ) );
}

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDSELECTOR_HPP_
