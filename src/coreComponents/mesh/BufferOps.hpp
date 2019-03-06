/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
