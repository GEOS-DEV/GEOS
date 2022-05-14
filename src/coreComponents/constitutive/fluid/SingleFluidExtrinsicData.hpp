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
 * @file SingleFluidExtrinsicData.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_SINGLEFLUIDEXTRINSICDATA_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_SINGLEFLUIDEXTRINSICDATA_HPP_

#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{

namespace extrinsicMeshData
{

namespace singlefluid
{

EXTRINSIC_MESH_DATA_TRAIT( density,
                           "density",
                           array2d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Density" );

EXTRINSIC_MESH_DATA_TRAIT( dDensity_dPressure,
                           "dDensity_dPressure",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of density with respect to pressure" );

EXTRINSIC_MESH_DATA_TRAIT( initialDensity,
                           "initialDensity",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Initial density" );

EXTRINSIC_MESH_DATA_TRAIT( density_n,
                           "density_n",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Density at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( viscosity,
                           "viscosity",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Viscosity" );

EXTRINSIC_MESH_DATA_TRAIT( dViscosity_dPressure,
                           "dViscosity_dPressure",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of viscosity with respect to pressure" );



}

}

}

#endif // GEOSX_CONSTITUTIVE_FLUID_SINGLEFLUIDEXTRINSICDATA_HPP_
