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
 * @file MPMSolverBaseFields.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_MPMSOLVERBASEFIELDS_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_MPMSOLVERBASEFIELDS_HPP_

#include "mesh/MeshFields.hpp"

namespace geos
{
/**
 * A scope for field traits.
 */
namespace fields
{

namespace mpm
{

DECLARE_FIELD( isBad,
               "isBad",
               array1d< int >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "An array that remembers particles that should be deleted at the end of the time step." );

DECLARE_FIELD( particleCrystalHealFlag,
               "particleCrystalHealFlag",
               array1d< int >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that remembers particles that are undergoing crystal healing for the mpm event." );

DECLARE_FIELD( particleMass,
               "particleMass",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that stores particle masses." );

DECLARE_FIELD( particleWavespeed,
               "particleWavespeed",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that stores particle wavespeeds." );

DECLARE_FIELD( particleHeatCapacity,
               "particleHeatCapacity",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that stores particle temperature." );

DECLARE_FIELD( particleTemperature,
               "particleTemperature",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that stores particle temperature." );

DECLARE_FIELD( particleInternalEnergy,
               "particleInternalEnergy",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that stores particle internal energy." );

DECLARE_FIELD( particleKineticEnergy,
               "particleKineticEnergy",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that stores particle kinetic energy." );

DECLARE_FIELD( particleArtificialViscosity,
               "particleArtificialViscosity",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that stores particle internal energy." );

DECLARE_FIELD( particleSPHJacobian,
               "particleSPHJacobian",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that stores particle SPH computed jacobian." );

DECLARE_FIELD( particleInitialVolume,
               "particleInitialVolume",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "ParticleInitialVolume" );

DECLARE_FIELD( particleInitialRVectors,
               "particleInitialRVectors",
               array3d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "ParticleInitialRVectors" );

DECLARE_FIELD( particleDeformationGradient,
               "particleDeformationGradient",
               array3d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleDeformationGradient" );

DECLARE_FIELD( particleFDot,
               "particleFDot",
               array3d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Material time derivative of the particle deformation gradient." );

DECLARE_FIELD( particleVelocityGradient,
               "particleVelocityGradient",
               array3d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "ParticleVelocityGradient" );

DECLARE_FIELD( particleStress,
               "particleStress",
               array2d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that holds particle stresses in Voigt notation." );

DECLARE_FIELD( particleBodyForce,
               "particleBodyForce",
               array2d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that holds particle body forces." );

DECLARE_FIELD( particlePlasticStrain,
               "particlePlasticStrain",
               array2d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that holds particle plastic strain in Voigt notation." );

DECLARE_FIELD( particleDensity,
               "particleDensity",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that holds particle densities." );

DECLARE_FIELD( particleDamageGradient,
               "particleDamageGradient",
               array2d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that holds particle damage gradients as calculated with an SPH kernel." );

DECLARE_FIELD( particleSurfaceFlag,
               "particleSurfaceFlag",
               array1d< int >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that holds particle surface flags." );

DECLARE_FIELD( particleSphF,
               "particleSphF",
               array3d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleSphF" );

DECLARE_FIELD( particleOverlap,
               "particleOverlap",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleOverlap" );

DECLARE_FIELD( particleReferencePosition,
               "particleReferencePosition",
               array2d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleReferencePosition" );

DECLARE_FIELD( particlePorosity,
               "particlePorosity",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticlePorosity" );

DECLARE_FIELD( particleReferencePorosity,
               "particleReferencePorosity",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleReferencePorosity" );

DECLARE_FIELD( particleCohesiveNormalForce,
               "particleCohesiveNormalForce",
               array2d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleCohesiveNormalForce" );

DECLARE_FIELD( particleCohesiveShearForce,
               "particleCohesiveShearForce",
               array2d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleCohesiveShearForce" );

DECLARE_FIELD( particleCohesiveZoneFlag, 
               "particleCohesiveZoneFlag", 
               array1d< int >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleCohesiveZoneFlag" );

DECLARE_FIELD( particleInitialMappedNodes, 
               "particleInitialMappedNodes", 
               array2d< globalIndex >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleInitialMappedNodes" );

DECLARE_FIELD( particleInitialShapeFunctionValues, 
               "particleInitialShapeFunctionValues", 
               array2d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleInitialShapeFunctionValues" );

DECLARE_FIELD( particleInitialShapeFunctionGradientValues, 
               "particleInitialShapeFunctionGradientValues", 
               array3d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleInitialShapeFunctionGradientValues" );

DECLARE_FIELD( particleCohesiveTraction,
               "particleCohesiveTraction",
               array2d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleCohesiveTraction" );

DECLARE_FIELD( particleAreaFromGrid, 
               "particleAreaFromGrid",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleAreaFromGrid" );

// DECLARE_FIELD( particleProjectedArea,
//                "particelProjectedArea",
//                array1d< real64 >,
//                0,
//                LEVEL_1,
//                WRITE_AND_READ,
//                "ParticleProjectedArea" );-

}

}

}

#endif // GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_MPMSOLVERBASEFIELDS_HPP_
