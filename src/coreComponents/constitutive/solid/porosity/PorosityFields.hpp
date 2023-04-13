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
 * @file PorosityFields.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_POROSITY_POROSITYFIELDS_HPP
#define GEOSX_CONSTITUTIVE_SOLID_POROSITY_POROSITYFIELDS_HPP

#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace porosity
{

DECLARE_FIELD( porosity,
               "porosity",
               array2d< real64 >,
               1.0, // important for newly created face elements
               LEVEL_0,
               WRITE_AND_READ,
               "Rock porosity" );

DECLARE_FIELD( porosity_n,
               "porosity_n",
               array2d< real64 >,
               1.0, // important for newly created face elements
               NOPLOT,
               WRITE_AND_READ,
               "Rock porosity at the previous converged time step" );

DECLARE_FIELD( dPorosity_dPressure,
               "dPorosity_dPressure",
               array2d< real64 >,
               0.0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of rock porosity with respect to pressure" );

DECLARE_FIELD( dPorosity_dTemperature,
               "dPorosity_dTemperature",
               array2d< real64 >,
               0.0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of rock porosity with respect to temperature" );

DECLARE_FIELD( initialPorosity,
               "initialPorosity",
               array2d< real64 >,
               0.0,
               NOPLOT,
               WRITE_AND_READ,
               "Initial porosity" );

DECLARE_FIELD( referencePorosity,
               "referencePorosity",
               array1d< real64 >,
               1.0,
               LEVEL_0,
               WRITE_AND_READ,
               "Reference porosity" );

DECLARE_FIELD( biotCoefficient,
               "biotCoefficient",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Biot coefficient" );

DECLARE_FIELD( thermalExpansionCoefficient,
               "thermalExpansionCoefficient",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Thermal expansion coefficient" );

}

}

}

#endif // GEOSX_CONSTITUTIVE_SOLID_POROSITY_POROSITYFIELDS_HPP_
