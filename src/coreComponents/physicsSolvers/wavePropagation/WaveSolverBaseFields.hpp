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

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION__HPP_WAVESOLVERBASEFIELDS
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION__HPP_WAVESOLVERBASEFIELDS

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace matricialFields
{ 
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

DECLARE_FIELD( StiffnessVector,
               "stiffnessVector",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Stiffness vector contains R_h*Pressure_n." );

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

DECLARE_FIELD( StiffnessVectorx,
               "stiffnessVectorx",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "x-component of stiffness vector." );

DECLARE_FIELD( StiffnessVectory,
               "stiffnessVectory",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "y-component of stiffness vector." );

DECLARE_FIELD( StiffnessVectorz,
               "stiffnessVectorz",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "z-component of stiffness vector." );

}

namespace geophysicalFields
{
 
 DECLARE_FIELD( Density,
               "Density",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Medium density of the cell" );

  DECLARE_FIELD( Pwavespeed,
               "Pwavespeed",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "P-waves speed in the cell" );

DECLARE_FIELD( Swavespeed,
               "Swavespeed",
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
}

namespace acousticSecondOrderSemFields
{
  DECLARE_FIELD( Pressure_nm1,
               "pressure_nm1",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar pressure at time n-1." );

  DECLARE_FIELD( Pressure_n,
               "pressure_n",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar pressure at time n." );

  DECLARE_FIELD( Pressure_np1,
               "pressure_np1",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Scalar pressure at time n+1." );

  DECLARE_FIELD( ForcingRHS,
               "rhs",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "RHS" );

  DECLARE_FIELD( AuxiliaryVar1PML,
               "auxiliaryVar1PML",
               array2d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "PML vectorial auxiliary variable 1." );

DECLARE_FIELD( AuxiliaryVar2PML,
               "auxiliaryVar2PML",
               array2d< real32 >,
               0,
               NOPLOT,
               NO_WRITE,
               "PML vectorial auxiliary variable 2." );

DECLARE_FIELD( AuxiliaryVar3PML,
               "auxiliaryVar3PML",
               array1d< real32 >,
               0,
               NOPLOT,
               NO_WRITE,
               "PML scalar auxiliary variable 3." );

DECLARE_FIELD( AuxiliaryVar4PML,
               "auxiliaryVar4PML",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "PML scalar auxiliary variable 4." );

DECLARE_FIELD( PressureDoubleDerivative,
               "pressureDoubleDerivative",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Double derivative of the pressure for each node to compute the gradient" );

DECLARE_FIELD( PartialGradient,
               "partialGradient",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Partiel gradient computed during backward propagation" );
}

namespace acousticFirstOrderSemFields
{
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

   DECLARE_FIELD( Pressure_np1,
               "pressure_np1",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Scalar pressure at time n+1." );
}

namespace elasticSecondOrderSemFields
{
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

DECLARE_FIELD( Displacementx_nm1,
               "displacementx_nm1",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "x-component of displacement at time n-1." );

DECLARE_FIELD( Displacementy_nm1,
               "displacementy_nm1",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "y-component of displacement at time n-1." );

DECLARE_FIELD( Displacementz_nm1,
               "displacementz_nm1",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "z-component of displacement at time n-1." );

DECLARE_FIELD( Displacementx_n,
               "displacementx_n",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "x-component of displacement at time n." );

DECLARE_FIELD( Displacementy_n,
               "displacementy_n",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "y-component of displacement at time n." );

DECLARE_FIELD( Displacementz_n,
               "displacementz_n",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "z-component of displacement at time n." );

DECLARE_FIELD( ForcingRHSx,
               "rhsx",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "RHS for x-direction" );

DECLARE_FIELD( ForcingRHSy,
               "rhsy",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "RHS for y-direction" );

DECLARE_FIELD( ForcingRHSz,
               "rhsz",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "RHS for z-direction" );

}

namespace elasticFirstOrderSemFields
{
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
        
}

}

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION__HPP_WAVESOLVERBASEFIELDS */
