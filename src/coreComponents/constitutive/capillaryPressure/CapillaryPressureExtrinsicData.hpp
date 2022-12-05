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
 * @file CapillaryPressureExtrinsicData.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSUREEXTRINSICDATA_HPP_
#define GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSUREEXTRINSICDATA_HPP_

#include "constitutive/capillaryPressure/layouts.hpp"
#include "mesh/ExtrinsicMeshData.hpp"

#include "dataRepository/BufferOps.hpp"

namespace geosx
{

namespace extrinsicMeshData
{

namespace cappres
{

using array3dLayoutCapPressure = array3d< real64, constitutive::cappres::LAYOUT_CAPPRES >;
using array4dLayoutCapPressure_dS = array4d< real64, constitutive::cappres::LAYOUT_CAPPRES_DS >;

  enum ModeIndexType : integer
  {
    DRAINAGE = 0,//to be used in array of Kernels
    IMBIBITION = 1,
    DRAINAGE_TO_IMBIBITION = 2,
    IMBIBITION_TO_DRAINAGE = 3
    };
//  };

//class ExpF{};
//inline std::enable_if< bufferOps::can_memcpy< ModeIndexType >, bool>::type
//inline std::enable_if< geosx::bufferOps::can_memcpy< ModeIndexType >, localIndex >::type test(){ return localIndex(0); };

//inline
//std::ostream & operator<<( std::ostream & os, ModeIndexType const & mode )
//{ return os << static_cast< int >( mode ); }


EXTRINSIC_MESH_DATA_TRAIT( phaseTrappedVolFraction,
                           "phaseTrappedVolFraction",
                           array3dLayoutCapPressure,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Phase trapped volume fraction" );

EXTRINSIC_MESH_DATA_TRAIT( phaseCapPressure,
                           "phaseCapPressure",
                           array3dLayoutCapPressure,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Phase capillary pressure" );

EXTRINSIC_MESH_DATA_TRAIT( dPhaseCapPressure_dPhaseVolFraction,
                           "dPhaseCapPressure_dPhaseVolFraction",
                           array4dLayoutCapPressure_dS,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of phase capillary pressure with respect to phase volume fraction" );

EXTRINSIC_MESH_DATA_TRAIT( jFuncMultiplier,
                           "jFuncMultiplier",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Multiplier for the Leverett J-function" );

EXTRINSIC_MESH_DATA_TRAIT( mode,
                           "mode",
                           array1d< integer >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Imbibition/drainage status mode" );
}

}

}

#endif // GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSUREEXTRINSICDATA_HPP_
