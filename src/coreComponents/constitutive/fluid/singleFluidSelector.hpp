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
 * @file singlePhaseSelector.hpp
 */
#ifndef GEOSX_CONSTITUTIVE_FLUID_SINGLEPHASESELECTOR_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_SINGLEPHASESELECTOR_HPP_

#include "constitutive/ConstitutivePassThruHandler.hpp"
#include "constitutive/fluid/CompressibleSinglePhaseFluid.hpp"

namespace geosx
{

namespace constitutive
{

template< typename LAMBDA >
void constitutiveUpdatePassThru( SingleFluidBase const & fluid,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< CompressibleSinglePhaseFluid >::execute( fluid, std::forward< LAMBDA >( lambda ) );
}

template< typename LAMBDA >
void constitutiveUpdatePassThru( SingleFluidBase & fluid,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< CompressibleSinglePhaseFluid >::execute( fluid, std::forward< LAMBDA >( lambda ) );
}

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_SINGLEPHASESELECTOR_HPP_
