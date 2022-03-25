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
 * @file BiotExtrinsicData.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_POROSITY_BIOTEXTRINSICDATA_HPP
#define GEOSX_CONSTITUTIVE_SOLID_POROSITY_BIOTEXTRINSICDATA_HPP

#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{

namespace extrinsicMeshData
{

namespace biot
{

EXTRINSIC_MESH_DATA_TRAIT(biotCoefficient,
                          "biotCoefficient",
                          array1d< real64 >,
                          -1,
                          LEVEL_0,
                          WRITE_AND_READ,
                          "Biot Coefficient" );




} // end namespace biot

} // end namespace extrinsicMeshData

} // end namespace geosx



#endif //GEOSX_CONSTITUTIVE_SOLID_POROSITY_BIOTEXTRINSICDATA_HPP