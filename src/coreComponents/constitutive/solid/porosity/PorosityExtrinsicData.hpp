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
 * @file PorosityExtrinsicData.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_POROSITY_POROSITYEXTRINSICDATA_HPP
#define GEOSX_CONSTITUTIVE_SOLID_POROSITY_POROSITYEXTRINSICDATA_HPP

#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{

namespace extrinsicMeshData
{

namespace porosity
{

EXTRINSIC_MESH_DATA_TRAIT( porosity,
                           "porosity",
                           array2d< real64 >,
                           1.0, // important for newly created face elements
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Rock porosity" );

EXTRINSIC_MESH_DATA_TRAIT( porosity_n,
                           "porosity_n",
                           array2d< real64 >,
                           1.0, // important for newly created face elements
                           NOPLOT,
                           WRITE_AND_READ,
                           "Rock porosity at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( dPorosity_dPressure,
                           "dPorosity_dPressure",
                           array2d< real64 >,
                           0.0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of rock porosity with respect to pressure" );

EXTRINSIC_MESH_DATA_TRAIT( initialPorosity,
                           "initialPorosity",
                           array2d< real64 >,
                           0.0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Initial porosity" );

EXTRINSIC_MESH_DATA_TRAIT( referencePorosity,
                           "referencePorosity",
                           array1d< real64 >,
                           1.0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Reference porosity" );

EXTRINSIC_MESH_DATA_TRAIT( grainBulkModulus,
                           "grainBulkModulus",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Grain bulk modulus" );

EXTRINSIC_MESH_DATA_TRAIT( biotCoefficient,
                           "biotCoefficient",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Biot coefficient" );

}

}

}

#endif // GEOSX_CONSTITUTIVE_SOLID_POROSITY_POROSITYEXTRINSICDATA_HPP_
