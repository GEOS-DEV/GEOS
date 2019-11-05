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
 * @file FaceElementStencil.cpp
 */

#include "FaceElementStencil.hpp"

namespace geosx
{

FaceElementStencil::FaceElementStencil():
  StencilBase<FaceElementStencil_Traits,FaceElementStencil>()
{}


void FaceElementStencil::add( localIndex const numPts,
                              localIndex  const * const elementRegionIndices,
                              localIndex  const * const elementSubRegionIndices,
                              localIndex  const * const elementIndices,
                              real64 const * const weights,
                              real64 const * const weightedElementCenterToConnectorCenter,
                              localIndex const connectorIndex )
{
  GEOS_ERROR_IF( numPts >= MAX_STENCIL_SIZE, "Maximum stencil size exceeded" );

  localIndex constexpr maxElems = FaceElementStencil::MAX_STENCIL_SIZE;
  stackArray1d<real64, FaceElementStencil::MAX_STENCIL_SIZE> zeroVector;
//  localIndex const maxElems = (numPts - 1 ) * (numPts - 1 );
//  zeroVector.resize(maxElems);
  zeroVector = 0;

  typename decltype( m_stencilIndices )::iterator iter = m_stencilIndices.find(connectorIndex);
  if( iter==m_stencilIndices.end() )
  {
    m_elementRegionIndices.appendArray( elementRegionIndices, numPts );
    m_elementSubRegionIndices.appendArray( elementSubRegionIndices, numPts );
    m_elementIndices.appendArray( elementIndices, numPts );
    m_weights.appendArray( weights, numPts );
    m_weightedElementCenterToConnectorCenter.appendArray( weightedElementCenterToConnectorCenter, numPts );
    m_stencilIndices[connectorIndex] = m_weights.size() - 1;
    m_buffer.appendArray( zeroVector.data(), maxElems );
  }
  else
  {
    localIndex const stencilIndex = iter->second;
    m_elementRegionIndices.clearArray( stencilIndex );
    m_elementSubRegionIndices.clearArray( stencilIndex );
    m_elementIndices.clearArray( stencilIndex );
    m_weights.clearArray( stencilIndex );
    m_weightedElementCenterToConnectorCenter.clearArray( stencilIndex );
    m_buffer.clearArray( stencilIndex );

    m_elementRegionIndices.appendToArray( stencilIndex, elementRegionIndices, numPts );
    m_elementSubRegionIndices.appendToArray( stencilIndex, elementSubRegionIndices, numPts );
    m_elementIndices.appendToArray( stencilIndex, elementIndices, numPts );
    m_weights.appendToArray( stencilIndex, weights, numPts );
    m_weightedElementCenterToConnectorCenter.appendToArray( stencilIndex, weightedElementCenterToConnectorCenter, numPts );
    m_buffer.appendToArray( stencilIndex, zeroVector.data(), maxElems );
  }
}


} /* namespace geosx */
