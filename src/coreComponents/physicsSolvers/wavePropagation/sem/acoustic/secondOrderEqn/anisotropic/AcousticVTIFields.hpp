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

DECLARE_FIELD_WITH_NAMESPACE( Delta,
                              "delta",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Delta thomsen anisotropy parameter",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( Epsilon,
                              "epsilon",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Epsilon thomsen anisotropy parameter",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( F,
                              "f",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "f quantity in VTI/TTI Fletcher's equations",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( StiffnessVector_p,
                              "stiffnessVector_p",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Stiffness vector contains R_h*Pressure_n.",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( StiffnessVector_q,
                              "stiffnessVector_q",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Stiffness vector contains R_h*Pressure_n.",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( StiffnessVectorA_p,
                              "stiffnessVectorA_p",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "P-type acoustic attenuation stiffness vector.",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( StiffnessVectorA_q,
                              "stiffnessVectorA_q",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Q-type acoustic attenuation stiffness vector.",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( Pressure_p_nm1,
                              "pressure_p_nm1",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Scalar pressure at time n-1.",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( Pressure_p_n,
                              "pressure_p_n",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Scalar pressure at time n.",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( Pressure_p_np1,
                              "pressure_p_np1",
                              array1d< real32 >,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Scalar pressure at time n+1.",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( Pressure_q_nm1,
                              "pressure_q_nm1",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Scalar auxiliary pressure q at time n-1.",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( Pressure_q_n,
                              "pressure_q_n",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Scalar auxiliary pressure q at time n.",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( Pressure_q_np1,
                              "pressure_q_np1",
                              array1d< real32 >,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Scalar auxiliary pressure q at time n+1.",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( DivPsi_p,
                              "divpsi_p",
                              array2d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "p-type memory variable for acoustic VTI attenuation.",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( DivPsi_q,
                              "divpsi_q",
                              array2d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "q-type memory variable for acoustic VTI attenuation.",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( DampingVector_p,
                              "dampingVector_p",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Diagonal of the Damping Matrix for p terms in p equation.",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( DampingVector_pq,
                              "dampingVector_pq",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Diagonal of the Damping Matrix for q terms in p equation.",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( DampingVector_q,
                              "dampingVector_q",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Diagonal of the Damping Matrix for q terms in q equation.",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( DampingVector_qp,
                              "dampingVector_qp",
                              array1d< real32 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Diagonal of the Damping Matrix for p terms in q equation.",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( LateralSurfaceFaceIndicator,
                              "lateralSurfaceFaceIndicator",
                              array1d< localIndex >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Free surface indicator, 1 if a face is on a lateral surface 0 otherwise.",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( LateralSurfaceNodeIndicator,
                              "lateralSurfaceNodeIndicator",
                              array1d< localIndex >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Lateral surface indicator, 1 if a face is on a lateral surface 0 otherwise.",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( BottomSurfaceFaceIndicator,
                              "bottomSurfaceFaceIndicator",
                              array1d< localIndex >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Bottom surface indicator, 1 if a face is on the bottom surface 0 otherwise.",
                              "acousticvtifields" );

DECLARE_FIELD_WITH_NAMESPACE( BottomSurfaceNodeIndicator,
                              "bottomSurfaceNodeIndicator",
                              array1d< localIndex >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Bottom surface indicator, 1 if a face is on the bottom surface 0 otherwise.",
                              "acousticvtifields" );
}

}

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION__HPP_ACOUSTICVTIFIELDS */
