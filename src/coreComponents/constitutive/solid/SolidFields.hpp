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
 * @file SolidFields.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_SOLIDFIELDS_HPP
#define GEOS_CONSTITUTIVE_SOLID_SOLIDFIELDS_HPP

#include "mesh/MeshFields.hpp"
#include "common/DataLayouts.hpp"

namespace geos
{

namespace fields
{

namespace solid
{

using array3dLayoutStress = array3d< real64, geos::solid::STRESS_PERMUTATION >;

DECLARE_FIELD( stress,
               "stress",
               array3dLayoutStress,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Current material stress" );

DECLARE_FIELD( oldStress,
               "oldStress",
               array3dLayoutStress,
               0,
               NOPLOT,
               NO_WRITE,
               "Previous material stress" );

DECLARE_FIELD( density,
               "density",
               array2d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Material density" );

DECLARE_FIELD( bulkModulus,
               "bulkModulus",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Bulk modulus" );

DECLARE_FIELD( shearModulus,
               "shearModulus",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Shear modulus" );

DECLARE_FIELD( biotCoefficient,
               "biotCoefficient",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Biot coefficient" );

DECLARE_FIELD( thermalExpansionCoefficient,
               "thermalExpansionCoefficient",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Thermal expansion coefficient" );

DECLARE_FIELD( yieldStress,
               "yieldStress",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Yield stress" );

DECLARE_FIELD( damage,
               "damage",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Material damage" );

DECLARE_FIELD( oldDamage,
               "oldDamage",
               array2d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Material old damage" );

DECLARE_FIELD( jacobian,
               "jacobian",
               array2d< real64 >,
               1.0,
               NOPLOT,
               NO_WRITE,
               "Jacobian of the deformation" );

DECLARE_FIELD( lengthScale,
               "lengthScale",
               array1d< real64 >,
               DBL_MIN,
               NOPLOT,
               NO_WRITE,
               "Length scale" );

DECLARE_FIELD( damageGrad,
               "damageGrad",
               array3d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Material damage gradient" );

DECLARE_FIELD( strainEnergyDensity,
               "strainEnergyDensity",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Material strain energy density" );

DECLARE_FIELD( volStrain,
               "volStrain",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Volumetric strain" );

DECLARE_FIELD( extDrivingForce,
               "extDrivingForce",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "External driving force" );

DECLARE_FIELD( criticalFractureEnergy,
               "criticalFractureEnergy",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Critical fracture energy" );

DECLARE_FIELD( tensileStrength,
               "tensileStrength",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Tensile strength" );

DECLARE_FIELD( compressStrength,
               "compressStrength",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Compressive strength" );

DECLARE_FIELD( deltaCoefficient,
               "deltaCoefficient",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Delta coefficient" );

DECLARE_FIELD( recompressionIndex,
               "recompressionIndex",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Recompression index" );

DECLARE_FIELD( virginCompressionIndex,
               "virginCompressionIndex",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Virgin compression index" );

DECLARE_FIELD( cslSlope,
               "cslSlope",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Critical state line slope" );

DECLARE_FIELD( shapeParameter,
               "shapeParameter",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Shape parameter" );

DECLARE_FIELD( preConsolidationPressure,
               "preConsolidationPressure",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Current preconsolidation pressure" );

DECLARE_FIELD( oldPreConsolidationPressure,
               "oldPreConsolidationPressure",
               array2d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Old preconsolidation pressure" );

DECLARE_FIELD( friction,
               "friction",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Friction angle" );

DECLARE_FIELD( initialFriction,
               "initialFriction",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Initial friction angle" );

DECLARE_FIELD( residualFriction,
               "residualFriction",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Final friction angle" );

DECLARE_FIELD( dilation,
               "dilation",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Dilation angle" );

DECLARE_FIELD( hardening,
               "hardening",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Hardening rate" );

DECLARE_FIELD( cohesion,
               "cohesion",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Cohesion" );

DECLARE_FIELD( oldCohesion,
               "oldCohesion",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Old cohesion" );

DECLARE_FIELD( pressureIntercept,
               "pressureIntercept",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Pressure intercept" );

DECLARE_FIELD( state,
               "state",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Material state" );

DECLARE_FIELD( oldState,
               "oldState",
               array2d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Old material state" );

// TODO merge in one
DECLARE_FIELD( c11,
               "c11",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "C11 component of Voigt stiffness tensor" );

DECLARE_FIELD( c12,
               "c12",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "C12 component of Voigt stiffness tensor" );

DECLARE_FIELD( c13,
               "c13",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "C13 component of Voigt stiffness tensor" );

DECLARE_FIELD( c22,
               "c22",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "C22 component of Voigt stiffness tensor" );

DECLARE_FIELD( c23,
               "c23",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "C23 component of Voigt stiffness tensor" );

DECLARE_FIELD( c33,
               "c33",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "C33 component of Voigt stiffness tensor" );

DECLARE_FIELD( c44,
               "c44",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "C44 component of Voigt stiffness tensor" );

DECLARE_FIELD( c55,
               "c55",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "C55 component of Voigt stiffness tensor" );

DECLARE_FIELD( c66,
               "c66",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "C66 component of Voigt stiffness tensor" );
//

DECLARE_FIELD( internalEnergy,
               "internalEnergy",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Internal energy of the solid per unit volume [J/m^3]" );

DECLARE_FIELD( oldInternalEnergy,
               "oldInternalEnergy",
               array2d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Internal energy of the solid per unit volume at the previous time-step [J/m^3]" );

DECLARE_FIELD( dInternalEnergy_dTemperature,
               "dInternalEnergy_dTemperature",
               array2d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of the solid internal energy w.r.t. temperature [J/(m^3.K)]" );

}

}

}

#endif // GEOS_CONSTITUTIVE_SOLID_SOLIDFIELDS_HPP_
