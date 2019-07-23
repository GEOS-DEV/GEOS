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

/**
 * @file CellElementStencilMPFA.cpp
 */


#include "CellElementStencilMPFA.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geosx
{

CellElementStencilMPFA::CellElementStencilMPFA():
  m_elementRegionIndices(),
  m_elementSubRegionIndices(),
  m_elementIndices(),
  m_weights(),
  m_connectorIndices()
{
  // TODO Auto-generated constructor stub

}

//CellElementStencilMPFA::~CellElementStencilMPFA()
//{
//  // TODO Auto-generated destructor stub
//}

void CellElementStencilMPFA::reserve( localIndex const size )
{
  m_elementRegionIndices.reserve( size * 9 );
  m_elementSubRegionIndices.reserve( size * 9 );
  m_elementIndices.reserve( size * 9 );
  m_weights.reserve( size * 9 );
}


void CellElementStencilMPFA::add( localIndex const numPts,
                                  localIndex const * const elementRegionIndices,
                                  localIndex const * const elementSubRegionIndices,
                                  localIndex const * const elementIndices,
                                  real64 const * const weights,
                                  localIndex const connectorIndex )
{
  GEOS_ERROR_IF( numPts >= MAX_STENCIL_SIZE, "Maximum stencil size exceeded" );

  m_elementRegionIndices.appendArray( elementRegionIndices, numPts );
  m_elementSubRegionIndices.appendArray( elementSubRegionIndices, numPts );
  m_elementIndices.appendArray( elementIndices, numPts );
  m_weights.appendArray( weights, numPts );

  m_connectorIndices[connectorIndex] = m_elementRegionIndices.size()-1;
}

bool CellElementStencilMPFA::zero( localIndex const connectorIndex )
{
  return
  executeOnMapValue( m_connectorIndices, connectorIndex, [&]( localIndex const connectionListIndex )
  {
    for (localIndex i = 0; i < m_weights.sizeOfArray(connectionListIndex) ; ++i)
    {
      m_weights[connectionListIndex][i] = 0;
    }
  });
}
} /* namespace geosx */
