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
 * @file MultiFluidExtrinsicData.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDEXTRINSICDATA_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDEXTRINSICDATA_HPP_

#include "constitutive/fluid/layouts.hpp"
#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{

namespace extrinsicMeshData
{

namespace multifluid
{

using array2dLayoutFluid = array2d< real64, constitutive::multifluid::LAYOUT_FLUID >;
using array3dLayoutFluid_dC = array3d< real64, constitutive::multifluid::LAYOUT_FLUID_DC >;
using array3dLayoutPhase = array3d< real64, constitutive::multifluid::LAYOUT_PHASE >;
using array4dLayoutPhase_dC = array4d< real64, constitutive::multifluid::LAYOUT_PHASE_DC >;
using array4dLayoutPhaseComp = array4d< real64, constitutive::multifluid::LAYOUT_PHASE_COMP >;
using array5dLayoutPhaseComp_dC = array5d< real64, constitutive::multifluid::LAYOUT_PHASE_COMP_DC >;

EXTRINSIC_MESH_DATA_TRAIT( phaseFraction,
                           "phaseFraction",
                           array3dLayoutPhase,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Phase fraction" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseFraction_dPressure,
                           "dPhaseFraction_dPressure",
                           array3dLayoutPhase,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase fraction with respect to pressure" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseFraction_dTemperature,
                           "dPhaseFraction_dTemperature",
                           array3dLayoutPhase,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase fraction with respect to temperature" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseFraction_dGlobalCompFraction,
                           "dPhaseFraction_dGlobalCompFraction",
                           array4dLayoutPhase_dC,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase fraction with respect to global component fraction" );

EXTRINSIC_MESH_DATA_TRAIT( phaseDensity,
                           "phaseDensity",
                           array3dLayoutPhase,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Phase density" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseDensity_dPressure,
                           "dPhaseDensity_dPressure",
                           array3dLayoutPhase,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase density with respect to pressure" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseDensity_dTemperature,
                           "dPhaseDensity_dTemperature",
                           array3dLayoutPhase,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase density with respect to temperature" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseDensity_dGlobalCompFraction,
                           "dPhaseDensity_dGlobalCompFraction",
                           array4dLayoutPhase_dC,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase density with respect to global component fraction" );

EXTRINSIC_MESH_DATA_TRAIT( phaseMassDensity,
                           "phaseMassDensity",
                           array3dLayoutPhase,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Phase mass density" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseMassDensity_dPressure,
                           "dPhaseMassDensity_dPressure",
                           array3dLayoutPhase,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase mass density with respect to pressure" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseMassDensity_dTemperature,
                           "dPhaseMassDensity_dTemperature",
                           array3dLayoutPhase,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase mass density with respect to temperature" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseMassDensity_dGlobalCompFraction,
                           "dPhaseMassDensity_dGlobalCompFraction",
                           array4dLayoutPhase_dC,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase mass density with respect to global component fraction" );

EXTRINSIC_MESH_DATA_TRAIT( phaseViscosity,
                           "phaseViscosity",
                           array3dLayoutPhase,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Phase viscosity" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseViscosity_dPressure,
                           "dPhaseViscosity_dPressure",
                           array3dLayoutPhase,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase viscosity with respect to pressure" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseViscosity_dTemperature,
                           "dPhaseViscosity_dTemperature",
                           array3dLayoutPhase,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase viscosity with respect to temperature" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseViscosity_dGlobalCompFraction,
                           "dPhaseViscosity_dGlobalCompFraction",
                           array4dLayoutPhase_dC,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase viscosity with respect to global component fraction" );

EXTRINSIC_MESH_DATA_TRAIT( phaseCompFraction,
                           "phaseCompFraction",
                           array4dLayoutPhaseComp,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Phase component fraction" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseCompFraction_dPressure,
                           "dPhaseCompFraction_dPressure",
                           array4dLayoutPhaseComp,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase component fraction with respect to pressure" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseCompFraction_dTemperature,
                           "dPhaseCompFraction_dTemperature",
                           array4dLayoutPhaseComp,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase component fraction with respect to temperature" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseCompFraction_dGlobalCompFraction,
                           "dPhaseCompFraction_dGlobalCompFraction",
                           array5dLayoutPhaseComp_dC,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase component fraction with respect to global component fraction" );

EXTRINSIC_MESH_DATA_TRAIT( totalDensity,
                           "totalDensity",
                           array2dLayoutFluid,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Total density" );

EXTRINSIC_MESH_DATA_TRAIT( initialTotalMassDensity,
                           "initialTotalMassDensity",
                           array2dLayoutFluid,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Initial total mass density" );

EXTRINSIC_MESH_DATA_TRAIT( dTotalDensity_dPressure,
                           "dTotalDensity_dPressure",
                           array2dLayoutFluid,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of total density with respect to pressure" );

EXTRINSIC_MESH_DATA_TRAIT( dTotalDensity_dTemperature,
                           "dTotalDensity_dTemperature",
                           array2dLayoutFluid,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of total density with respect to temperature" );

EXTRINSIC_MESH_DATA_TRAIT( dTotalDensity_dGlobalCompFraction,
                           "dTotalDensity_dGlobalCompFraction",
                           array3dLayoutFluid_dC,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of total density with respect to global component fraction" );


}

}

}

#endif // GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDEXTRINSICDATA_HPP_
