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
 * @file AcousticVTIFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICVTIFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICVTIFIELDS_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"


namespace geos
{

namespace fields
{

namespace acousticvtifields
{

DECLARE_FIELD( AcousticVelocity,
               "acousticVelocity",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Medium velocity of the cell" );

DECLARE_FIELD( AcousticDensity,
               "acousticDensity",
               array1d< real32 >,
               1,
               NOPLOT,
               WRITE_AND_READ,
               "Medium density of the cell" );

DECLARE_FIELD( AcousticDelta,
               "acousticDelta",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Delta thomsen anisotropy parameter" );

DECLARE_FIELD( AcousticEpsilon,
               "acousticEpsilon",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Epsilon thomsen anisotropy parameter" );

DECLARE_FIELD( AcousticSigma,
               "acousticSigma",
               array1d< real32 >,
               0.75,
               NOPLOT,
               WRITE_AND_READ,
               "Sigma quantity in VTI/TTI Fletcher's equations" );

DECLARE_FIELD( AcousticMassVector,
               "acousticMassVector",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Mass Matrix." );

DECLARE_FIELD( StiffnessVector_p,
               "stiffnessVector_p",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Stiffness vector contains R_h*Pressure_n." );

DECLARE_FIELD( StiffnessVector_q,
               "stiffnessVector_q",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Stiffness vector contains R_h*Pressure_n." );

DECLARE_FIELD( Pressure_p_nm1,
               "pressure_p_nm1",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar pressure at time n-1." );

DECLARE_FIELD( Pressure_p_n,
               "pressure_p_n",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar pressure at time n." );

DECLARE_FIELD( Pressure_p_np1,
               "pressure_p_np1",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Scalar pressure at time n+1." );

DECLARE_FIELD( Pressure_q_nm1,
               "pressure_q_nm1",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar auxiliary pressure q at time n-1." );

DECLARE_FIELD( Pressure_q_n,
               "pressure_q_n",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar auxiliary pressure q at time n." );

DECLARE_FIELD( Pressure_q_np1,
               "pressure_q_np1",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Scalar auxiliary pressure q at time n+1." );

DECLARE_FIELD( DampingVector_p,
               "dampingVector_p",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Damping Matrix for p terms in p equation." );

DECLARE_FIELD( DampingVector_pq,
               "dampingVector_pq",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Damping Matrix for q terms in p equation." );

DECLARE_FIELD( DampingVector_q,
               "dampingVector_q",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Damping Matrix for q terms in q equation." );

DECLARE_FIELD( DampingVector_qp,
               "dampingVector_qp",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Damping Matrix for p terms in q equation." );

DECLARE_FIELD( AcousticFreeSurfaceFaceIndicator,
               "acousticFreeSurfaceFaceIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Free surface indicator, 1 if a face is on free surface 0 otherwise." );

DECLARE_FIELD( AcousticFreeSurfaceNodeIndicator,
               "acousticFreeSurfaceNodeIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Free surface indicator, 1 if a node is on free surface 0 otherwise." );

DECLARE_FIELD( AcousticLateralSurfaceFaceIndicator,
               "lateralSurfaceFaceIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Free surface indicator, 1 if a face is on a lateral surface 0 otherwise." );

DECLARE_FIELD( AcousticLateralSurfaceNodeIndicator,
               "lateralSurfaceNodeIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Lateral surface indicator, 1 if a face is on a lateral surface 0 otherwise." );

DECLARE_FIELD( AcousticBottomSurfaceFaceIndicator,
               "bottomSurfaceFaceIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Bottom surface indicator, 1 if a face is on the bottom surface 0 otherwise." );

DECLARE_FIELD( AcousticBottomSurfaceNodeIndicator,
               "bottomSurfaceNodeIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Bottom surface indicator, 1 if a face is on the bottom surface 0 otherwise." );

DECLARE_FIELD( ForcingRHS,
               "rhs",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "RHS" );
}

}

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION__HPP_ACOUSTICVTIFIELDS */
