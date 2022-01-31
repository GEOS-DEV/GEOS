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
 * @file thermalConductivitySelector.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYSELECTOR_HPP
#define GEOSX_CONSTITUTIVE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYSELECTOR_HPP

#include "constitutive/ConstitutivePassThruHandler.hpp"
#include "constitutive/thermalConductivity/ConstantThermalConductivity.hpp"
#include "constitutive/thermalConductivity/VolumeWeightedThermalConductivity.hpp"

namespace geosx
{

namespace constitutive
{

template< typename LAMBDA >
void constitutiveUpdatePassThru( ThermalConductivityBase const & thermalConductivity,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< ConstantThermalConductivity,
                               VolumeWeightedThermalConductivity >::execute( thermalConductivity, std::forward< LAMBDA >( lambda ) );
}

template< typename LAMBDA >
void constitutiveUpdatePassThru( ThermalConductivityBase & thermalConductivity,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< ConstantThermalConductivity,
                               VolumeWeightedThermalConductivity >::execute( thermalConductivity, std::forward< LAMBDA >( lambda ) );
}

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYSELECTOR_HPP
