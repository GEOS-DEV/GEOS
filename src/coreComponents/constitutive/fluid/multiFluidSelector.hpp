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
#ifndef GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDSELECTOR_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDSELECTOR_HPP_

#include "constitutive/ConstitutivePassThruHandler.hpp"
#include "constitutive/fluid/DeadOilFluid.hpp"
#include "constitutive/fluid/BlackOilFluid.hpp"
#include "constitutive/fluid/MultiPhaseMultiComponentFluid.hpp"

#include "common/GeosxConfig.hpp"
#ifdef GEOSX_USE_PVTPackage
#include "constitutive/fluid/CompositionalMultiphaseFluid.hpp"
#endif

namespace geosx
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
                               CO2BrineEzrokhiFluid >::execute( fluid, std::forward< LAMBDA >( lambda ) );
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
                               CO2BrineEzrokhiFluid >::execute( fluid, std::forward< LAMBDA >( lambda ) );
}

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDSELECTOR_HPP_
