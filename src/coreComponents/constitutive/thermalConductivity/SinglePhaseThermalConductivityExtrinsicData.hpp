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
 * @file SinglePhaseThermalConductivityExtrinsicData.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SINGLEPHASE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYEXTRINSICDATA_HPP_
#define GEOSX_CONSTITUTIVE_SINGLEPHASE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYEXTRINSICDATA_HPP_

#include "constitutive/relativePermeability/layouts.hpp"
#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{

namespace extrinsicMeshData
{

namespace thermalconductivity
{

EXTRINSIC_MESH_DATA_TRAIT( dEffectiveConductivity_dPorosity,
                           "dEffectiveConductivity_dPorosity",
                           array3d< real64 >,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of effective conductivity with respect to porosity" );

}

}

}

#endif // GEOSX_CONSTITUTIVE_SINGLEPHASE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYEXTRINSICDATA_HPP_
