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
 * @file FlowSolverBaseExtrinsicData.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERBASEEXTRINSICDATA_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERBASEEXTRINSICDATA_HPP_

#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{
/**
 * A scope for extrinsic mesh data traits.
 */
namespace extrinsicMeshData
{

namespace flow
{

EXTRINSIC_MESH_DATA_TRAIT( pressure,
                           "pressure",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Pressure" );

EXTRINSIC_MESH_DATA_TRAIT( pressure_n,
                           "pressure_n",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Pressure at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( initialPressure,
                           "initialPressure",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Initial pressure" );

EXTRINSIC_MESH_DATA_TRAIT( deltaPressure,
                           "deltaPressure",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Delta pressure: current pressure - initial pressure" );

EXTRINSIC_MESH_DATA_TRAIT( facePressure,
                           "facePressure",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Face pressure" );

EXTRINSIC_MESH_DATA_TRAIT( facePressure_n,
                           "facePressure_n",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Face pressure at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( temperature,
                           "temperature",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Temperature" );

EXTRINSIC_MESH_DATA_TRAIT( faceTemperature,
                           "faceTemperature",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Face temperature" );

EXTRINSIC_MESH_DATA_TRAIT( temperature_n,
                           "temperature_n",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Temperature at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( deltaVolume,
                           "deltaVolume",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Accumulated volume updates" );

EXTRINSIC_MESH_DATA_TRAIT( aperture0,
                           "aperture_n",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Initial aperture" );

EXTRINSIC_MESH_DATA_TRAIT( hydraulicAperture,
                           "hydraulicAperture",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Hydraulic aperture" );

EXTRINSIC_MESH_DATA_TRAIT( gravityCoefficient,
                           "gravityCoefficient",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Gravity coefficient (dot product of gravity acceleration by gravity vector)" );

EXTRINSIC_MESH_DATA_TRAIT( mimGravityCoefficient,
                           "mimGravityCoefficient",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Mimetic gravity coefficient" );

EXTRINSIC_MESH_DATA_TRAIT( bcPressure,
                           "bcPressure",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Boundary condition pressure" );

EXTRINSIC_MESH_DATA_TRAIT( bcTemperature,
                           "bcTemperature",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Boundary temperature" );

}

}

}

#endif // GEOSX_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERBASEEXTRINSICDATA_HPP_
