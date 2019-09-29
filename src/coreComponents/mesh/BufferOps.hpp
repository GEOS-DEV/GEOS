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


#ifndef MESH_BUFFEROPS_H_
#define MESH_BUFFEROPS_H_

#include "common/DataTypes.hpp"

class OrderedVariableToManyElementRelation;
class FixedToManyElementRelation;
class ElementRegionManager;

namespace geosx
{
namespace bufferOps
{

template< bool DO_PACKING >
localIndex Pack( char*& buffer,
                 OrderedVariableToManyElementRelation const & var,
                 arrayView1d<localIndex const> const & packList,
                 ElementRegionManager const * const elementRegionManager );

template< bool DO_PACKING >
localIndex Pack( char*& buffer,
                 FixedToManyElementRelation const & var,
                 arrayView1d<localIndex const> const & packList,
                 ElementRegionManager const * const elementRegionManager );



localIndex Unpack( char const * & buffer,
                   OrderedVariableToManyElementRelation & var,
                   arrayView1d<localIndex const> const & packList,
                   ElementRegionManager const * const elementRegionManager,
                   bool const clearFlag );

localIndex Unpack( char const * & buffer,
                   FixedToManyElementRelation & var,
                   arrayView1d<localIndex const> const & packList,
                   ElementRegionManager const * const elementRegionManager,
                   bool const clearFlag );
}
}
#endif
