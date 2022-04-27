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
 * @file ReactiveCompositionalMultiphaseOBLExtrinsicData.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_REACTIVECOMPOSITIONALMULTIPHASEOBLEXTRINSICDATA_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_REACTIVECOMPOSITIONALMULTIPHASEOBLEXTRINSICDATA_HPP_

#include "common/DataLayouts.hpp"
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

using array2dLayoutComp = array2d< real64, compflow::LAYOUT_COMP >;
using array2dLayoutOBLOpVals = array2d< real64, compflow::LAYOUT_OBL_OPERATOR_VALUES >;
using array3dLayoutOBLOpDers = array3d< real64, compflow::LAYOUT_OBL_OPERATOR_DERIVATIVES >;

EXTRINSIC_MESH_DATA_TRAIT( globalCompFraction_n,
                           "globalCompFraction_n",
                           array2dLayoutComp,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Global component fraction at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( bcGlobalCompFraction,
                           "bcGlobalCompFraction",
                           array2dLayoutComp,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Boundary condition global component fraction" );

EXTRINSIC_MESH_DATA_TRAIT( initialTemperature,
                           "initialTemperature",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Initial temperature" );

EXTRINSIC_MESH_DATA_TRAIT( referencePorosity,
                           "referencePorosity",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Reference porosity" );

EXTRINSIC_MESH_DATA_TRAIT( referencePoreVolume,
                           "referencePoreVolume",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Reference volume of porous space" );

EXTRINSIC_MESH_DATA_TRAIT( referenceRockVolume,
                           "referenceRockVolume",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Reference volume of rock" );

EXTRINSIC_MESH_DATA_TRAIT( rockVolumetricHeatCapacity,
                           "rockVolumetricHeatCapacity",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Volumetric rock heat capacity" );

EXTRINSIC_MESH_DATA_TRAIT( rockThermalConductivity,
                           "rockThermalConductivity",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Rock thermal conductivity" );

EXTRINSIC_MESH_DATA_TRAIT( rockKineticRateFactor,
                           "rockKineticRateFactor",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Rock kinetic rate factor" );

EXTRINSIC_MESH_DATA_TRAIT( OBLOperatorValues,
                           "OBLOperatorValues",
                           array2dLayoutOBLOpVals,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Values of OBL operators" );


EXTRINSIC_MESH_DATA_TRAIT( OBLOperatorDerivatives,
                           "OBLOperatorDerivatives",
                           array3dLayoutOBLOpDers,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivatives of OBL operators with respect to all primary variables (Dofs)" );

EXTRINSIC_MESH_DATA_TRAIT( OBLOperatorValues_n,
                           "OBLOperatorValues_n",
                           array2dLayoutOBLOpVals,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Values of OBL operators at the previous converged time step" );

}

}

}

#endif // GEOSX_PHYSICSSOLVERS_FLUIDFLOW_REACTIVECOMPOSITIONALMULTIPHASEOBLEXTRINSICDATA_HPP_
