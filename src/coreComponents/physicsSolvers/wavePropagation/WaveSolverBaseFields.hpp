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
 * @file AcousticFirstOrderWaveEquationSEM.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION__HPP_WAVESOLVERBASEFIELDS
#define GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION__HPP_WAVESOLVERBASEFIELDS

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geosx
{

namespace fields
{

namespace wavesolverfields
{

DECLARE_FIELD( Pressure_np1,
               "pressure_np1",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Scalar pressure at time n+1." );

DECLARE_FIELD( Velocity_x,
               "velocity_x",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Velocity in the x-direction." );

DECLARE_FIELD( Velocity_y,
               "velocity_y",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Velocity in the y-direction." );

DECLARE_FIELD( Velocity_z,
               "velocity_z",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Velocity in the z-direction." );

DECLARE_FIELD( ForcingRHS,
               "rhs",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "RHS" );

DECLARE_FIELD( MassVector,
               "massVector",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Mass Matrix." );

DECLARE_FIELD( DampingVector,
               "dampingVector",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Damping Matrix." );

DECLARE_FIELD( MediumVelocity,
               "mediumVelocity",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Medium velocity of the cell" );

DECLARE_FIELD( MediumDensity,
               "mediumDensity",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Medium density of the cell" );

DECLARE_FIELD( FreeSurfaceFaceIndicator,
               "freeSurfaceFaceIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Free surface indicator, 1 if a face is on free surface 0 otherwise." );

DECLARE_FIELD( FreeSurfaceNodeIndicator,
               "freeSurfaceNodeIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Free surface indicator, 1 if a node is on free surface 0 otherwise." );

DECLARE_FIELD( Displacementx_np1,
               "displacementx_np1",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "x-component of displacement at time n+1." );

DECLARE_FIELD( Displacementy_np1,
               "displacementy_np1",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "y-component of displacement at time n+1." );

DECLARE_FIELD( Displacementz_np1,
               "displacementz_np1",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "z-component of displacement at time n+1." );

DECLARE_FIELD( Stresstensorxx,
               "stresstensorxx",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "xx-components of the stress tensor." );

DECLARE_FIELD( Stresstensoryy,
               "stresstensoryy",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "yy-components of the stress tensor." );

DECLARE_FIELD( Stresstensorzz,
               "stresstensorzz",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "zz-components of the stress tensor." );

DECLARE_FIELD( Stresstensorxy,
               "stresstensorxy",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "xy-components of the stress tensor (symetric of yx-component)." );

DECLARE_FIELD( Stresstensorxz,
               "stresstensorxz",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "xz-components of the stress tensor (symetric of zx-component)." );

DECLARE_FIELD( Stresstensoryz,
               "stresstensoryz",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "yz-components of the stress tensor (symetric of zy-component)." );

DECLARE_FIELD( DampingVectorx,
               "dampingVectorx",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal Damping Matrix in x-direction." );

DECLARE_FIELD( DampingVectory,
               "dampingVectory",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal Damping Matrix in y-direction." );

DECLARE_FIELD( DampingVectorz,
               "dampingVectorz",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal Damping Matrix in z-direction." );

DECLARE_FIELD( MediumVelocityVp,
               "mediumVelocityVp",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "P-waves speed in the cell" );

DECLARE_FIELD( MediumVelocityVs,
               "mediumVelocityVs",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "S-waves speed in the cell" );

DECLARE_FIELD( Lambda,
               "lambda",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "First Lame parameter: lambda" );

DECLARE_FIELD( Mu,
               "mu",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Second Lame parameter: mu" );


}

}

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_AcousticFirstOrderWaveEquationSEM_HPP_ */
