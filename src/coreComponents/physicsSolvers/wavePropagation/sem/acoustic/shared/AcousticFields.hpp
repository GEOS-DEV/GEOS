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
 * @file AcousticFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICFIELDS_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"


namespace geos
{

namespace fields
{

namespace acousticfields
{

DECLARE_FIELD_WITH_NAMESPACE( Pressure_nm1,
                              "pressure_nm1",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Scalar pressure at time n-1.",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( Pressure_n,
                              "pressure_n",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Scalar pressure at time n.",
                              "acousticfields" );


DECLARE_FIELD_WITH_NAMESPACE( Pressure_np1,
                              "pressure_np1",
                              array1d< real32 >,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Scalar pressure at time n+1.",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( DivPsi,
                              "divpsi",
                              array2d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "memory variable for acoustic attenuation.",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( DivPsi_p,
                              "divpsi_p",
                              array2d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "p-type memory variable for acoustic VTI attenuation.",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( DivPsi_q,
                              "divpsi_q",
                              array2d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "q-type memory variable for acoustic VTI attenuation.",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( PressureForward,
                              "pressureForward",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Pressure field from forward pass on each node to compute the gradient",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( Velocity_x,
                              "velocity_x",
                              array2d< real32 >,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Velocity in the x-direction.",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( Velocity_y,
                              "velocity_y",
                              array2d< real32 >,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Velocity in the y-direction.",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( Velocity_z,
                              "velocity_z",
                              array2d< real32 >,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Velocity in the z-direction.",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( PartialGradient,
                              "partialGradient",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Partial gradient or imaging condition computed during backward propagation",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( PartialGradient2,
                              "partialGradient2",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Partial gradient for density/velocity computed during backward propagation",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( ForcingRHS,
                              "rhs",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "RHS",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( AcousticMassVector,
                              "acousticMassVector",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Diagonal of the Mass Matrix.",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( StiffnessVector,
                              "stiffnessVector",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Stiffness vector contains R_h*Pressure_n.",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( StiffnessVectorA,
                              "stiffnessVectorA",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Acoustic attenuation stiffness vector.",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( DampingVector,
                              "dampingVector",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Diagonal of the Damping Matrix.",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( AcousticVelocity,
                              "acousticVelocity",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Medium velocity of the cell",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( AcousticDensity,
                              "acousticDensity",
                              array1d< real32 >,
                              1,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Medium density of the cell",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( AcousticQualityFactor,
                              "acousticQualityFactor",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Quality factor for acoustic wave attenuation in the cell",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( AcousticFreeSurfaceFaceIndicator,
                              "acousticFreeSurfaceFaceIndicator",
                              array1d< localIndex >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Free surface indicator, 1 if a face is on free surface 0 otherwise.",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( AcousticFreeSurfaceNodeIndicator,
                              "acousticFreeSurfaceNodeIndicator",
                              array1d< localIndex >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Free surface indicator, 1 if a node is on free surface 0 otherwise.",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( AuxiliaryVar1PML,
                              "auxiliaryVar1PML",
                              array2d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "PML vectorial auxiliary variable 1.",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( AuxiliaryVar2PML,
                              "auxiliaryVar2PML",
                              array2d< real32 >,
                              0,
                              NOPLOT,
                              NO_WRITE,
                              "PML vectorial auxiliary variable 2.",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( AuxiliaryVar3PML,
                              "auxiliaryVar3PML",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              NO_WRITE,
                              "PML scalar auxiliary variable 3.",
                              "acousticfields" );

DECLARE_FIELD_WITH_NAMESPACE( AuxiliaryVar4PML,
                              "auxiliaryVar4PML",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "PML scalar auxiliary variable 4.",
                              "acousticfields" );

}

}

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION__HPP_ACOUSTICFIELDS */
