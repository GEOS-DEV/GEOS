/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FlowSolverBaseFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERBASEFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERBASEFIELDS_HPP_

#include "mesh/MeshFields.hpp"

namespace geos
{
/**
 * A scope for field traits.
 */
namespace fields
{

namespace flow
{

DECLARE_FIELD_WITH_NAMESPACE( pressure,
                              "pressure",
                              array1d< real64 >,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Pressure",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( pressure_n,
                              "pressure_n",
                              array1d< real64 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Pressure at the previous converged time step",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( pressure_k,
                              "pressure_k",
                              array1d< real64 >,
                              0,
                              NOPLOT,
                              NO_WRITE,
                              "Pressure at the previous sequential iteration",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( initialPressure,
                              "initialPressure",
                              array1d< real64 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Initial pressure",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( deltaPressure,
                              "deltaPressure",
                              array1d< real64 >,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Delta pressure: current pressure - initial pressure",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( facePressure,
                              "facePressure",
                              array1d< real64 >,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Face pressure",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( facePressure_n,
                              "facePressure_n",
                              array1d< real64 >,
                              0,
                              NOPLOT,
                              NO_WRITE,
                              "Face pressure at the previous converged time step",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( pressureGradient,
                              "pressureGradient",
                              array2d< real64 >,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Pressure gradient",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( temperature,
                              "temperature",
                              array1d< real64 >,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Temperature",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( temperature_n,
                              "temperature_n",
                              array1d< real64 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Temperature at the previous converged time step",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( temperature_k,
                              "temperature_k",
                              array1d< real64 >,
                              0,
                              NOPLOT,
                              NO_WRITE,
                              "Temperature at the previous sequential iteration",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( initialTemperature,
                              "initialTemperature",
                              array1d< real64 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Initial temperature",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( faceTemperature,
                              "faceTemperature",
                              array1d< real64 >,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Face temperature",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( netToGross,
                              "netToGross",
                              array1d< real64 >,
                              1,
                              NOPLOT,
                              NO_WRITE,
                              "Net to gross",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( deltaVolume,
                              "deltaVolume",
                              array1d< real64 >,
                              0,
                              NOPLOT,
                              NO_WRITE,
                              "Accumulated volume updates",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( aperture0,
                              "aperture_n",
                              array1d< real64 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Initial aperture",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( hydraulicAperture,
                              "hydraulicAperture",
                              array1d< real64 >,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Hydraulic aperture",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( gravityCoefficient,
                              "gravityCoefficient",
                              array1d< real64 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Gravity coefficient (dot product of gravity acceleration by gravity vector)",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( mimGravityCoefficient,
                              "mimGravityCoefficient",
                              array1d< real64 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Mimetic gravity coefficient",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( macroElementIndex,
                              "macroElementIndex",
                              array1d< integer >,
                              -1,
                              LEVEL_1,
                              WRITE_AND_READ,
                              "Index of the macroelement for a given element",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( bcPressure,
                              "bcPressure",
                              array1d< real64 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Boundary condition pressure",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( bcTemperature,
                              "bcTemperature",
                              array1d< real64 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Boundary condition temperature",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( elementStabConstant,
                              "elementStabConstant",
                              array1d< real64 >,
                              0,
                              LEVEL_1,
                              WRITE_AND_READ,
                              "Stabilization constant for pressure jump stabilization",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( transMultiplier,
                              "permeabilityTransMultiplier",
                              array1d< real64 >,
                              1,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Permeability transmissibility multipliers",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( pressureScalingFactor,
                              "pressureScalingFactor",
                              array1d< real64 >,
                              1,
                              NOPLOT,
                              NO_WRITE,
                              "Scaling factors for pressure",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( temperatureScalingFactor,
                              "temperatureScalingFactor",
                              array1d< real64 >,
                              1,
                              NOPLOT,
                              NO_WRITE,
                              "Scaling factors for temperature",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( energy,
                              "energy",
                              array1d< real64 >,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Energy",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( energy_n,
                              "energy_n",
                              array1d< real64 >,
                              0,
                              NOPLOT,
                              NO_WRITE,
                              "Energy at the previous converged time step",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( fractureCreationTime,
                              "fractureCreationTime",
                              array1d< real64 >,
                              0,
                              LEVEL_1,
                              WRITE_AND_READ,
                              "The creation time for the fracture cell.",
                              "flow" );
}

}

}

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERBASEFIELDS_HPP_
