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
 * @file WellSolverBaseFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLSOLVERBASEFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLSOLVERBASEFIELDS_HPP_

#include "mesh/MeshFields.hpp"

namespace geos
{
/**
 * A scope for field traits.
 */
namespace fields
{

namespace well
{

DECLARE_FIELD( pressure,
               "pressure",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Pressure" );

DECLARE_FIELD( pressure_n,
               "pressure_n",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Pressure at the previous converged time step" );

DECLARE_FIELD( temperature,
               "wellTemperature",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Temperature" );

DECLARE_FIELD( temperature_n,
               "wellTemperature_n",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Temperature at the previous converged time step" );

DECLARE_FIELD( gravityCoefficient,
               "gravityCoefficient",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Gravity coefficient (dot product of gravity acceleration by gravity vector)" );

}

}

}

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLSOLVERBASEFIELDS_HPP_
