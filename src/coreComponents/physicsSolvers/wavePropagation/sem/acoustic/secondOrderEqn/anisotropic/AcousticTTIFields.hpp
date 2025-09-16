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
 * @file AcousticTTIFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICTTIFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICTTIFIELDS_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace acousticttifields
{
DECLARE_FIELD( AcousticDipX,
               "acousticDipX",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Inline Dip for for pseudo-acoustic TTI" );

DECLARE_FIELD( AcousticDipY,
               "acousticDipY",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Crossline Dip for pseudo-acoustic TTI" );


DECLARE_FIELD( AcousticDofTilt,
               "acousticDofTilt",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Tilt anisotropy Dof parameter" );

DECLARE_FIELD( AcousticDofAzimuth,
               "acousticDofAzimuth",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Azimuth anisotropy Dof parameter" );
}

}

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION__HPP_ACOUSTICTTIFIELDS */
