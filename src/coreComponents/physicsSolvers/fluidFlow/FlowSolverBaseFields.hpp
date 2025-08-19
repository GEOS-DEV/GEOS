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

DECLARE_FIELD( pressure_k,
               "pressure_k",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Pressure at the previous sequential iteration" );

DECLARE_FIELD( initialPressure,
               "initialPressure",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Initial pressure" );

DECLARE_FIELD( deltaPressure,
               "deltaPressure",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Delta pressure: current pressure - initial pressure" );

DECLARE_FIELD( facePressure,
               "facePressure",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Face pressure" );

DECLARE_FIELD( facePressure_n,
               "facePressure_n",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Face pressure at the previous converged time step" );

DECLARE_FIELD( pressureGradient,
               "pressureGradient",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Pressure gradient" );

DECLARE_FIELD( temperature,
               "temperature",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Temperature" );

DECLARE_FIELD( temperature_n,
               "temperature_n",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Temperature at the previous converged time step" );

DECLARE_FIELD( temperature_k,
               "temperature_k",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Temperature at the previous sequential iteration" );

DECLARE_FIELD( initialTemperature,
               "initialTemperature",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Initial temperature" );

DECLARE_FIELD( faceTemperature,
               "faceTemperature",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Face temperature" );

DECLARE_FIELD( netToGross,
               "netToGross",
               array1d< real64 >,
               1,
               NOPLOT,
               NO_WRITE,
               "Net to gross" );

DECLARE_FIELD( deltaVolume,
               "deltaVolume",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Accumulated volume updates" );

DECLARE_FIELD( aperture0,
               "aperture_n",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Initial aperture" );

DECLARE_FIELD( hydraulicAperture,
               "hydraulicAperture",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Hydraulic aperture" );

DECLARE_FIELD( gravityCoefficient,
               "gravityCoefficient",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Gravity coefficient (dot product of gravity acceleration by gravity vector)" );

DECLARE_FIELD( mimGravityCoefficient,
               "mimGravityCoefficient",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Mimetic gravity coefficient" );

DECLARE_FIELD( macroElementIndex,
               "macroElementIndex",
               array1d< integer >,
               -1,
               LEVEL_1,
               WRITE_AND_READ,
               "Index of the macroelement for a given element" );

DECLARE_FIELD( bcPressure,
               "bcPressure",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Boundary condition pressure" );

DECLARE_FIELD( bcTemperature,
               "bcTemperature",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Temperature boundary condition" );

// ========================================
// New Flux and Phase Property Fields
// ========================================

// Total flux fields
DECLARE_FIELD( totalMassFlux,
               "totalMassFlux",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Total mass flux across all phases" );

DECLARE_FIELD( totalCapillaryFlux,
               "totalCapillaryFlux",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Total capillary flux across all phases" );

// Phase-specific flux fields
DECLARE_FIELD( phaseMassicFlux,
               "phaseMassicFlux",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Massic flux for each phase" );

DECLARE_FIELD( phaseCapillaryFlux,
               "phaseCapillaryFlux",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Capillary flux for each phase" );

DECLARE_FIELD( phaseBuoyancyFlux,
               "phaseBuoyancyFlux",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Buoyancy flux for each phase" );

// Fractional flow and mobility fields
DECLARE_FIELD( fractionalFlow,
               "fractionalFlow",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Fractional flow for each phase" );

DECLARE_FIELD( totalMobility,
               "totalMobility",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Total mobility across all phases" );

DECLARE_FIELD( phaseMobility,
               "phaseMobility",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Mobility for each phase" );

// Phase saturation fields (for multiphase flow)
DECLARE_FIELD( phaseSaturation,
               "phaseSaturation",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Saturation for each phase" );

DECLARE_FIELD( phaseSaturation_n,
               "phaseSaturation_n",
               array2d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase saturation at previous time step" );

DECLARE_FIELD( phaseSaturation_k,
               "phaseSaturation_k",
               array2d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Phase saturation at previous sequential iteration" );

// Mixture property fields
DECLARE_FIELD( mixtureDensity,
               "mixtureDensity",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Mixture density" );

DECLARE_FIELD( mixtureEnthalpy,
               "mixtureEnthalpy",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Mixture enthalpy" );

DECLARE_FIELD( mixtureViscosity,
               "mixtureViscosity",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Mixture viscosity" );

// Phase density fields
DECLARE_FIELD( phaseDensity,
               "phaseDensity",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Density for each phase" );

DECLARE_FIELD( phaseViscosity,
               "phaseViscosity",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Viscosity for each phase" );

DECLARE_FIELD( phaseEnthalpy,
               "phaseEnthalpy",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Enthalpy for each phase" );

// Capillary pressure fields
DECLARE_FIELD( capillaryPressure,
               "capillaryPressure",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Capillary pressure between phases" );

// Transport equation residuals
DECLARE_FIELD( transportResidual,
               "transportResidual",
               array2d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Transport equation residuals" );

// Flux derivatives for Jacobian assembly
DECLARE_FIELD( dTotalMassFlux_dPressure,
               "dTotalMassFlux_dPressure",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of total mass flux w.r.t. pressure" );

DECLARE_FIELD( dTotalMassFlux_dSaturation,
               "dTotalMassFlux_dSaturation",
               array2d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of total mass flux w.r.t. saturation" );

DECLARE_FIELD( dFractionalFlow_dSaturation,
               "dFractionalFlow_dSaturation",
               array3d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of fractional flow w.r.t. saturation" );

// ...existing code...
}

}

}

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERBASEFIELDS_HPP_
