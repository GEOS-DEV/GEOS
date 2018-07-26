/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

/*
 * @file StencilCollection.cpp
 *
 */

#include "StencilCollection.hpp"

namespace geosx
{

void StencilCollection::set(localIndex index,
                            CellDescriptor *connCells,
                            array<CellDescriptor> cells,
                            array<real64> weights)
{
  assert(index < numConnections());
  CellConnection & conn = m_faceConnectors[index];
  conn.connectedCellIndices[0] = connCells[0];
  conn.connectedCellIndices[1] = connCells[1];
  conn.stencilCellIndices = cells;
  conn.stencilWeights = weights;
}

}