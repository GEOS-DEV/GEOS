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
 * @file inducedSeismicityFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_INDUCEDSEISMICITYFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_INDUCEDSEISMICITYFIELDS_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{
/**
 * A scope for field traits.
 */
namespace fields
{

namespace inducedSeismicity
{

DECLARE_FIELD( t_a,
               "t_a",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Dieterich constitutive relaxation time of seismicity rate" );

DECLARE_FIELD( aSigma,
               "aSigma",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Product of direct effect parameter a and initial normal stress" );

DECLARE_FIELD( initialPressure,
               "initialPressure",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Initial pore pressure");

DECLARE_FIELD( pressure,
               "pressure",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Pore pressure");

DECLARE_FIELD( pressure_n,
               "pressure_n",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Pore pressure at the previous converged time step");

DECLARE_FIELD( pressureRate,
               "pressureRate",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Pore pressure rate");

DECLARE_FIELD( initialNormalStress,
               "initialNormalStress",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Initial normal stress acting on the fault");

DECLARE_FIELD( normalStress,
               "normalStress",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Normal stress acting on the fault");

DECLARE_FIELD( normalStress_n,
               "normalStress_n",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Normal stress acting on the fault at the previous converged time step");

DECLARE_FIELD( normalStressRate,
               "normalStressRate",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Normal stress rate acting on the fault");
               
DECLARE_FIELD( initialShearStress,
               "initialShearStress",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Initial shear stress acting on the fault");

DECLARE_FIELD( shearStress,
               "shearStress",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Shear stress acting on the fault");

DECLARE_FIELD( shearStress_n,
               "shearStress_n",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Shear stress acting on the fault at the previous converged time step");

DECLARE_FIELD( shearStressRate,
               "shearStressRate",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Shear stress rate acting on the fault");

DECLARE_FIELD( seismicityRate,
               "seismicityRate",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Seismicity rate");

DECLARE_FIELD( h,
               "h",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Log of the seismicity rate");

DECLARE_FIELD( h_n,
               "h_n",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Log of the seismicity rate at the previous converged time step");

DECLARE_FIELD( logDenom,
               "logDenom",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Log of the denominator of the integral form of the seismicity rate");

DECLARE_FIELD( logDenom_n,
               "logDenom_n",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Log of the denominator of the integral form of the seismicity rate at the previous converged time step");

}

}

}

#endif // GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_INDUCEDSEISMICITYFIELDS_HPP_
