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
 * @file CompositionalMultiphaseWellExtrinsicData.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELLEXTRINSICDATA_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELLEXTRINSICDATA_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{
/**
 * A scope for extrinsic mesh data traits.
 */
namespace extrinsicMeshData
{

namespace well
{

using array2dLayoutFluid_dC = array2d< real64, compflow::LAYOUT_FLUID_DC >;
using array2dLayoutPhase = array2d< real64, compflow::LAYOUT_PHASE >;
using array3dLayoutPhase_dC = array3d< real64, compflow::LAYOUT_PHASE_DC >;
using array2dLayoutComp = array2d< real64, compflow::LAYOUT_COMP >;
using array3dLayoutComp_dC = array3d< real64, compflow::LAYOUT_COMP_DC >;
using array3dLayoutPhaseComp = array3d< real64, compflow::LAYOUT_PHASE_COMP >;

EXTRINSIC_MESH_DATA_TRAIT( temperature,
                           "wellTemperature",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Temperature" );

EXTRINSIC_MESH_DATA_TRAIT( temperatureOld,
                           "wellTemperatureOld",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "TemperatureOld" );

EXTRINSIC_MESH_DATA_TRAIT( globalCompDensity,
                           "globalCompDensity",
                           array2dLayoutComp,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Global component density" );

EXTRINSIC_MESH_DATA_TRAIT( globalCompDensityOld,
                           "globalCompDensityOld",
                           array2dLayoutComp,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Global component density at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( mixtureConnectionRate,
                           "wellElementMixtureConnectionRate",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Mixture connection rate" );

EXTRINSIC_MESH_DATA_TRAIT( mixtureConnectionRateOld,
                           "wellElementMixtureConnectionRateOld",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Mixture connection rate at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( globalCompFraction,
                           "globalCompFraction",
                           array2dLayoutComp,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Global component fraction" );

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

EXTRINSIC_MESH_DATA_TRAIT( dPhaseVolumeFraction_dGlobalCompDensity,
                           "dPhaseVolumeFraction_dGlobalCompDensity",
                           array3dLayoutPhase_dC,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase volume fraction with respect to global component density" );

EXTRINSIC_MESH_DATA_TRAIT( phaseVolumeFractionOld,
                           "phaseVolumeFractionOld",
                           array2dLayoutPhase,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Phase volume fraction at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( phaseDensityOld,
                           "phaseDensityOld",
                           array2dLayoutPhase,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Phase density at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( totalDensityOld,
                           "totalDensityOld",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Total density at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( phaseComponentFractionOld,
                           "phaseComponentFractionOld",
                           array3dLayoutPhaseComp,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Phase component fraction at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( totalMassDensity,
                           "totalMassDensity",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Total mass density" );

EXTRINSIC_MESH_DATA_TRAIT( dTotalMassDensity_dPressure,
                           "dTotalMassDensity_dPressure",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of total mass density with respect to pressure" );

EXTRINSIC_MESH_DATA_TRAIT( dTotalMassDensity_dGlobalCompDensity,
                           "dTotalMassDensity_dComp", // to avoid a rebaseline
                           array2dLayoutFluid_dC,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of total mass density with respect to global component density" );

EXTRINSIC_MESH_DATA_TRAIT( compPerforationRate,
                           "compPerforationRate",
                           array2d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Component perforation rate" );

EXTRINSIC_MESH_DATA_TRAIT( dCompPerforationRate_dPres,
                           "dCompPerforationRate_dPres",
                           array3d< real64 >,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of component perforation rate with respect to pressure" );

EXTRINSIC_MESH_DATA_TRAIT( dCompPerforationRate_dComp,
                           "dCompPerforationRate_dComp",
                           array4d< real64 >,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of component perforation rate with respect to global component density" );


}

}

}

#endif // GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELLEXTRINSICDATA_HPP_
