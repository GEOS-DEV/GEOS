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
 * @file multiPhaseThermalConductivitySelector.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_MULTIPHASE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYSELECTOR_HPP
#define GEOSX_CONSTITUTIVE_MULTIPHASE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYSELECTOR_HPP

#include "constitutive/ConstitutivePassThruHandler.hpp"
#include "constitutive/thermalConductivity/MultiPhaseConstantThermalConductivity.hpp"
#include "constitutive/thermalConductivity/MultiPhaseVolumeWeightedThermalConductivity.hpp"

namespace geos
{

namespace constitutive
{

template< typename LAMBDA >
void constitutiveUpdatePassThru( MultiPhaseThermalConductivityBase const & thermalConductivity,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< MultiPhaseConstantThermalConductivity,
                               MultiPhaseVolumeWeightedThermalConductivity >::execute( thermalConductivity, std::forward< LAMBDA >( lambda ) );
}

template< typename LAMBDA >
void constitutiveUpdatePassThru( MultiPhaseThermalConductivityBase & thermalConductivity,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< MultiPhaseConstantThermalConductivity,
                               MultiPhaseVolumeWeightedThermalConductivity >::execute( thermalConductivity, std::forward< LAMBDA >( lambda ) );
}

} // namespace constitutive

} // namespace geos

#endif //GEOSX_CONSTITUTIVE_MULTIPHASE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYSELECTOR_HPP
