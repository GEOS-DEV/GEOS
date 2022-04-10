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

EXTRINSIC_MESH_DATA_TRAIT( dPhaseFraction,
                           "dPhaseFraction",
                           array4dLayoutPhase_dC,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase fraction with respect to pressure, temperature, and global component fractions" );

EXTRINSIC_MESH_DATA_TRAIT( phaseDensity,
                           "phaseDensity",
                           array3dLayoutPhase,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Phase density" );

EXTRINSIC_MESH_DATA_TRAIT( phaseDensityOld,
                           "phaseDensityOld",
                           array3dLayoutPhase,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Phase density at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseDensity,
                           "dPhaseDensity",
                           array4dLayoutPhase_dC,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase density with respect to pressure, temperature, and global component fractions" );

EXTRINSIC_MESH_DATA_TRAIT( phaseMassDensity,
                           "phaseMassDensity",
                           array3dLayoutPhase,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Phase mass density" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseMassDensity,
                           "dPhaseMassDensity",
                           array4dLayoutPhase_dC,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase mass density with respect to pressure, temperature, and global component fractions" );

EXTRINSIC_MESH_DATA_TRAIT( phaseViscosity,
                           "phaseViscosity",
                           array3dLayoutPhase,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Phase viscosity" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseViscosity,
                           "dPhaseViscosity",
                           array4dLayoutPhase_dC,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase viscosity with respect to pressure, temperature, and global component fractions" );

EXTRINSIC_MESH_DATA_TRAIT( phaseEnthalpy,
                           "phaseEnthalpy",
                           array3dLayoutPhase,
                           0,
                           NOPLOT, // default behavior overridden by thermal models
                           NO_WRITE,
                           "Phase enthalpy" );

EXTRINSIC_MESH_DATA_TRAIT( phaseEnthalpyOld,
                           "phaseEnthalpyOld",
                           array3dLayoutPhase,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Phase enthalpy at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseEnthalpy,
                           "dPhaseEnthalpy",
                           array4dLayoutPhase_dC,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase enthalpy with respect to pressure, temperature, and global component fractions" );

EXTRINSIC_MESH_DATA_TRAIT( phaseInternalEnergy,
                           "phaseInternalEnergy",
                           array3dLayoutPhase,
                           0,
                           NOPLOT, // default behavior overridden by thermal models
                           NO_WRITE,
                           "Phase internal energy" );

EXTRINSIC_MESH_DATA_TRAIT( phaseInternalEnergyOld,
                           "phaseInternalEnergyOld",
                           array3dLayoutPhase,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Phase internal energy at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseInternalEnergy,
                           "dPhaseInternalEnergy",
                           array4dLayoutPhase_dC,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase internal energy with respect to pressure, temperature, and global component fractions" );

EXTRINSIC_MESH_DATA_TRAIT( phaseCompFraction,
                           "phaseCompFraction",
                           array4dLayoutPhaseComp,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Phase component fraction" );

EXTRINSIC_MESH_DATA_TRAIT( phaseCompFractionOld,
                           "phaseCompFractionOld",
                           array4dLayoutPhaseComp,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Phase component fraction at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseCompFraction,
                           "dPhaseCompFraction",
                           array5dLayoutPhaseComp_dC,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of phase component fraction with respect to pressure, temperature, and global component fractions" );

EXTRINSIC_MESH_DATA_TRAIT( totalDensity,
                           "totalDensity",
                           array2dLayoutFluid,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Total density" );

EXTRINSIC_MESH_DATA_TRAIT( totalDensityOld,
                           "totalDensityOld",
                           array2dLayoutFluid,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Total density at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( initialTotalMassDensity,
                           "initialTotalMassDensity",
                           array2dLayoutFluid,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Initial total mass density" );

EXTRINSIC_MESH_DATA_TRAIT( dTotalDensity,
                           "dTotalDensity",
                           array3dLayoutFluid_dC,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of total density with respect to pressure, temperature, and global component fractions" );

}

}

}

#endif // GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDEXTRINSICDATA_HPP_
