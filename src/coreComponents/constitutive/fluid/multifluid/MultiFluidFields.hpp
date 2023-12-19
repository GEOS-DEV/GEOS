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
 * @file MultiFluidFields.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDFIELDS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDFIELDS_HPP_

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace multifluid
{

using array2dLayout = Array< real64, 2 >;
using array3dLayout = Array< real64, 3 >;
using array4dLayout = Array< real64, 4 >;
using array5dLayout = Array< real64, 5 >;

DECLARE_FIELD( phaseFraction,
               "phaseFraction",
               array3dLayout,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Phase fraction" );

DECLARE_FIELD( dPhaseFraction,
               "dPhaseFraction",
               array4dLayout,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase fraction with respect to pressure, temperature, and global component fractions" );

DECLARE_FIELD( phaseDensity,
               "phaseDensity",
               array3dLayout,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Phase density" );

DECLARE_FIELD( phaseDensity_n,
               "phaseDensity_n",
               array3dLayout,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase density at the previous converged time step" );

DECLARE_FIELD( dPhaseDensity,
               "dPhaseDensity",
               array4dLayout,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase density with respect to pressure, temperature, and global component fractions" );

DECLARE_FIELD( phaseMassDensity,
               "phaseMassDensity",
               array3dLayout,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Phase mass density" );

DECLARE_FIELD( dPhaseMassDensity,
               "dPhaseMassDensity",
               array4dLayout,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase mass density with respect to pressure, temperature, and global component fractions" );

DECLARE_FIELD( phaseViscosity,
               "phaseViscosity",
               array3dLayout,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Phase viscosity" );

DECLARE_FIELD( dPhaseViscosity,
               "dPhaseViscosity",
               array4dLayout,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase viscosity with respect to pressure, temperature, and global component fractions" );

DECLARE_FIELD( phaseEnthalpy,
               "phaseEnthalpy",
               array3dLayout,
               0,
               NOPLOT, // default behavior overridden by thermal models
               NO_WRITE,
               "Phase enthalpy" );

DECLARE_FIELD( phaseEnthalpy_n,
               "phaseEnthalpy_n",
               array3dLayout,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase enthalpy at the previous converged time step" );

DECLARE_FIELD( dPhaseEnthalpy,
               "dPhaseEnthalpy",
               array4dLayout,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase enthalpy with respect to pressure, temperature, and global component fractions" );

DECLARE_FIELD( phaseInternalEnergy,
               "phaseInternalEnergy",
               array3dLayout,
               0,
               NOPLOT, // default behavior overridden by thermal models
               NO_WRITE,
               "Phase internal energy" );

DECLARE_FIELD( phaseInternalEnergy_n,
               "phaseInternalEnergy_n",
               array3dLayout,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase internal energy at the previous converged time step" );

DECLARE_FIELD( dPhaseInternalEnergy,
               "dPhaseInternalEnergy",
               array4dLayout,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase internal energy with respect to pressure, temperature, and global component fractions" );

DECLARE_FIELD( phaseCompFraction,
               "phaseCompFraction",
               array4dLayout,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Phase component fraction" );

DECLARE_FIELD( phaseCompFraction_n,
               "phaseCompFraction_n",
               array4dLayout,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase component fraction at the previous converged time step" );

DECLARE_FIELD( dPhaseCompFraction,
               "dPhaseCompFraction",
               array5dLayout,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase component fraction with respect to pressure, temperature, and global component fractions" );

DECLARE_FIELD( totalDensity,
               "totalDensity",
               array2dLayout,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Total density" );

DECLARE_FIELD( totalDensity_n,
               "totalDensity_n",
               array2dLayout,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Total density at the previous converged time step" );

DECLARE_FIELD( dTotalDensity,
               "dTotalDensity",
               array3dLayout,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of total density with respect to pressure, temperature, and global component fractions" );

}

}

}

#endif // GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDFIELDS_HPP_
