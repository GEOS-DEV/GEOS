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

DECLARE_FIELD( particleMass,
               "particleMass",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that stores particle masses." );

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

DECLARE_FIELD( particleDensity,
               "particleDensity",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that holds particle densities." );

// This lives on particle subregions directly for now
// DECLARE_FIELD( particleDamage,
//                "particleDamage",
//                array1d< real64 >,
//                0,
//                LEVEL_1,
//                WRITE_AND_READ,
//                "An array that holds particle damage values." );

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
}

}

}

#endif // GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_MPMSOLVERBASEFIELDS_HPP_
