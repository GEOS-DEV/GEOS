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
 * @file SlurryFluidExtrinsicData.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_SLURRYFLUIDEXTRINSICDATA_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_SLURRYFLUIDEXTRINSICDATA_HPP_

#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{

namespace extrinsicMeshData
{

namespace slurryfluid
{

EXTRINSIC_MESH_DATA_TRAIT( density,
                           "density",
                           array2d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Density" );

EXTRINSIC_MESH_DATA_TRAIT( dDensity_dPressure,
                           "dDens_dPres",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of density with respect to pressure" );

EXTRINSIC_MESH_DATA_TRAIT( dDensity_dProppantConcentration,
                           "dDens_dProppantConc",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of density with respect to proppant concentration" );

EXTRINSIC_MESH_DATA_TRAIT( dDensity_dComponentConcentration,
                           "dDens_dCompConc",
                           array3d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of density with respect to component concentration" );

EXTRINSIC_MESH_DATA_TRAIT( componentDensity,
                           "componentDensity",
                           array3d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Component density" );

EXTRINSIC_MESH_DATA_TRAIT( dComponentDensity_dPressure,
                           "dCompDens_dPres",
                           array3d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of component density with respect to pressure" );

EXTRINSIC_MESH_DATA_TRAIT( dComponentDensity_dComponentConcentration,
                           "dCompDens_dCompConc",
                           array4d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of component density with respect to component concentration" );

EXTRINSIC_MESH_DATA_TRAIT( fluidDensity,
                           "FluidDensity",
                           array2d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Fluid density" );

EXTRINSIC_MESH_DATA_TRAIT( dFluidDensity_dPressure,
                           "dFluidDens_dPres",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of fluid density with respect to pressure" );

EXTRINSIC_MESH_DATA_TRAIT( dFluidDensity_dComponentConcentration,
                           "dFluidDens_dCompConc",
                           array3d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of fluid density with respect to component concentration" );

EXTRINSIC_MESH_DATA_TRAIT( viscosity,
                           "viscosity",
                           array2d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Viscosity" );

EXTRINSIC_MESH_DATA_TRAIT( dViscosity_dPressure,
                           "dVisc_dPres",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of viscosity with respect to pressure" );

EXTRINSIC_MESH_DATA_TRAIT( dViscosity_dProppantConcentration,
                           "dVisc_dProppantConc",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of viscosity with respect to proppant concentration" );

EXTRINSIC_MESH_DATA_TRAIT( dViscosity_dComponentConcentration,
                           "dVisc_dCompConc",
                           array3d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of viscosity with respect to component concentration" );

EXTRINSIC_MESH_DATA_TRAIT( fluidViscosity,
                           "FluidViscosity",
                           array2d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Fluid viscosity" );

EXTRINSIC_MESH_DATA_TRAIT( dFluidViscosity_dPressure,
                           "dFluidVisc_dPres",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of fluid viscosity with respect to pressure" );

EXTRINSIC_MESH_DATA_TRAIT( dFluidViscosity_dComponentConcentration,
                           "dFluidVisc_dCompConc",
                           array3d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of fluid viscosity with respect to component concentration" );


}

}

}

#endif // GEOSX_CONSTITUTIVE_FLUID_SLURRYFLUIDEXTRINSICDATA_HPP_
