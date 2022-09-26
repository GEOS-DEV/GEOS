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
 * @file SolidExtrinsicData.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_SOLIDEXTRINSICDATA_HPP
#define GEOSX_CONSTITUTIVE_SOLID_SOLIDEXTRINSICDATA_HPP

#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{

namespace extrinsicMeshData
{

namespace solid
{

EXTRINSIC_MESH_DATA_TRAIT( bulkModulus,
                           "bulkModulus",
                           array1d< real64 >,
                           -1,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Elastic Bulk Modulus Field" );

EXTRINSIC_MESH_DATA_TRAIT( shearModulus,
                           "shearModulus",
                           array1d< real64 >,
                           -1,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Elastic Shear Modulus Field" );



} // end namespace solid

} // end namespace extrinsicMeshData

} // end namespace geosx



#endif //GEOSX_CONSTITUTIVE_SOLID_SOLIDEXTRINSICDATA_HPP
