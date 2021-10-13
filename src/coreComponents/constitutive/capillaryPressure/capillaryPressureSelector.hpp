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
 * @file capillaryPressureSelector.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSURESELECTOR_HPP
#define GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSURESELECTOR_HPP

#include "constitutive/ConstitutivePassThruHandler.hpp"
#include "constitutive/capillaryPressure/BrooksCoreyCapillaryPressure.hpp"
#include "constitutive/capillaryPressure/TableCapillaryPressure.hpp"
#include "constitutive/capillaryPressure/VanGenuchtenCapillaryPressure.hpp"

namespace geosx
{

namespace constitutive
{

template< typename LAMBDA >
void constitutiveUpdatePassThru( CapillaryPressureBase const & capPres,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< BrooksCoreyCapillaryPressure,
                               TableCapillaryPressure,
                               VanGenuchtenCapillaryPressure >::execute( capPres, std::forward< LAMBDA >( lambda ) );
}

template< typename LAMBDA >
void constitutiveUpdatePassThru( CapillaryPressureBase & capPres,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< BrooksCoreyCapillaryPressure,
                               TableCapillaryPressure,
                               VanGenuchtenCapillaryPressure >::execute( capPres, std::forward< LAMBDA >( lambda ) );
}

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSURESELECTOR_HPP
