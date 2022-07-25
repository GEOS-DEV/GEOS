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
 * @file CompositionalMultiphaseBaseExtrinsicData.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEEXTRINSICDATA_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEEXTRINSICDATA_HPP_

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

using array2dLayoutPhase = array2d< real64, compflow::LAYOUT_PHASE >;
using array3dLayoutPhase_dC = array3d< real64, compflow::LAYOUT_PHASE_DC >;
using array2dLayoutComp = array2d< real64, compflow::LAYOUT_COMP >;
using array3dLayoutComp_dC = array3d< real64, compflow::LAYOUT_COMP_DC >;
using array3dLayoutPhaseComp = array3d< real64, compflow::LAYOUT_PHASE_COMP >;

EXTRINSIC_MESH_DATA_TRAIT( globalCompDensity,
                           "globalCompDensity",
                           array2dLayoutComp,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Global component density" );

EXTRINSIC_MESH_DATA_TRAIT( globalCompDensity_n,
                           "globalCompDensity_n",
                           array2dLayoutComp,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Global component density updates at the previous converged time step " );

EXTRINSIC_MESH_DATA_TRAIT( globalCompFraction,
                           "globalCompFraction",
                           array2dLayoutComp,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Global component fraction" );

EXTRINSIC_MESH_DATA_TRAIT( faceGlobalCompFraction,
                           "faceGlobalCompFraction",
                           array2dLayoutComp,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Face global component fraction" );

EXTRINSIC_MESH_DATA_TRAIT( dGlobalCompFraction_dGlobalCompDensity,
                           "dGlobalCompFraction_dGlobalCompDensity",
                           array3dLayoutComp_dC,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of global component fraction with respect to component density" );

EXTRINSIC_MESH_DATA_TRAIT( phaseVolumeFraction,
                           "phaseVolumeFraction",
                           array2dLayoutPhase,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Phase volume fraction" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseVolumeFraction_dPressure,
                           "dPhaseVolumeFraction_dPressure",
                           array2dLayoutPhase,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase volume fraction with respect to pressure" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseVolumeFraction_dTemperature,
                           "dPhaseVolumeFraction_dTemperature",
                           array2dLayoutPhase,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase volume fraction with respect to temperature" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseVolumeFraction_dGlobalCompDensity,
                           "dPhaseVolumeFraction_dGlobalCompDensity",
                           array3dLayoutPhase_dC,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase volume fraction with respect to global component density" );

EXTRINSIC_MESH_DATA_TRAIT( phaseMobility,
                           "phaseMobility",
                           array2dLayoutPhase,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Phase mobility" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseMobility_dPressure,
                           "dPhaseMobility_dPressure",
                           array2dLayoutPhase,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase mobility with respect to pressure" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseMobility_dTemperature,
                           "dPhaseMobility_dTemperature",
                           array2dLayoutPhase,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase mobility with respect to temperature" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseMobility_dGlobalCompDensity,
                           "dPhaseMobility_dGlobalCompDensity",
                           array3dLayoutPhase_dC,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase volume fraction with respect to global component density" );

EXTRINSIC_MESH_DATA_TRAIT( phaseVolumeFraction_n,
                           "phaseVolumeFraction_n",
                           array2dLayoutPhase,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Phase volume fraction at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( phaseMobility_n,
                           "phaseMobility_n",
                           array2dLayoutPhase,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Phase mobility at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( phaseOutflux,
                           "phaseOutflux",
                           array2dLayoutPhase,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Phase outflux" );

EXTRINSIC_MESH_DATA_TRAIT( componentOutflux,
                           "componentOutflux",
                           array2dLayoutComp,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Component outflux" );

EXTRINSIC_MESH_DATA_TRAIT( phaseCFLNumber,
                           "phaseCFLNumber",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           NO_WRITE,
                           "Phase CFL number" );

EXTRINSIC_MESH_DATA_TRAIT( componentCFLNumber,
                           "componentCFLNumber",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           NO_WRITE,
                           "Component CFL number" );


}

}

}

#endif // GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEEXTRINSICDATA_HPP_
