/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MeshData.hpp
 */
#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_MESHDATA_HPP
#define GEOSX_LINEARALGEBRA_MULTISCALE_MESHDATA_HPP

#include "mesh/MeshFields.hpp"

namespace geos
{
namespace fields
{
namespace multiscale
{

// Doxygen throws spurious warnings even though these decls are identical to ones in MeshFields
///@cond DO_NOT_DOCUMENT

DECLARE_FIELD( OrigElementRegion,
               "origElementRegion",
               array1d< localIndex >,
               -1,
               LEVEL_0,
               NO_WRITE,
               "Region index of the original geosx element" );

DECLARE_FIELD( OrigElementSubRegion,
               "origElementSubRegion",
               array1d< localIndex >,
               -1,
               LEVEL_0,
               NO_WRITE,
               "Subregion index of the original geosx element" );

DECLARE_FIELD( OrigElementIndex,
               "origElementIndex",
               array1d< localIndex >,
               -1,
               LEVEL_0,
               NO_WRITE,
               "Index of the original geosx element" );

DECLARE_FIELD( OrigNodeIndex,
               "origNodeIndex",
               array1d< localIndex >,
               -1,
               LEVEL_0,
               NO_WRITE,
               "Index of the original geosx node" );

DECLARE_FIELD( CoarseCellLocalIndex,
               "coarseCellLocalIndex",
               array1d< localIndex >,
               -1,
               LEVEL_0,
               NO_WRITE,
               "Local index of the coarse scale cell" );

DECLARE_FIELD( CoarseCellGlobalIndex,
               "coarseCellGlobalIndex",
               array1d< globalIndex >,
               -1,
               LEVEL_0,
               NO_WRITE,
               "Global index of the coarse scale cell" );

DECLARE_FIELD( CoarseNodeLocalIndex,
               "coarseNodeLocalIndex",
               array1d< localIndex >,
               -1,
               LEVEL_0,
               NO_WRITE,
               "Local index of the coarse scale node" );

DECLARE_FIELD( CoarseNodeGlobalIndex,
               "coarseNodeGlobalIndex",
               array1d< globalIndex >,
               -1,
               LEVEL_0,
               NO_WRITE,
               "Global index of the coarse scale node" );

DECLARE_FIELD( FineNodeLocalIndex,
               "fineNodeLocalIndex",
               array1d< localIndex >,
               -1,
               LEVEL_0,
               NO_WRITE,
               "Local index of the fine scale node corresponding to a coarse node" );

///@endcond DO_NOT_DOCUMENT

} // namespace multiscale
} // namespace fields
} // namespace geos

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_MESHDATA_HPP
