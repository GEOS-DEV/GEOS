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
                              localIndex const connectorIndex )
{
  GEOSX_ERROR_IF( numPts >= MAX_STENCIL_SIZE, "Maximum stencil size exceeded" );

  typename decltype( m_connectorIndices )::iterator iter = m_connectorIndices.find(connectorIndex);
  if( iter==m_connectorIndices.end() )
  {
    m_elementRegionIndices.appendArray( elementRegionIndices, numPts );
    m_elementSubRegionIndices.appendArray( elementSubRegionIndices, numPts );
    m_elementIndices.appendArray( elementIndices, numPts );
    m_weights.appendArray( weights, numPts );

    m_connectorIndices[connectorIndex] = m_weights.size() - 1;
  }
  else
  {
    localIndex const stencilIndex = iter->second;
    m_elementRegionIndices.clearArray( stencilIndex );
    m_elementSubRegionIndices.clearArray( stencilIndex );
    m_elementIndices.clearArray( stencilIndex );
    m_weights.clearArray( stencilIndex );

    m_elementRegionIndices.appendToArray( stencilIndex, elementRegionIndices, numPts );
    m_elementSubRegionIndices.appendToArray( stencilIndex, elementSubRegionIndices, numPts );
    m_elementIndices.appendToArray( stencilIndex, elementIndices, numPts );
    m_weights.appendToArray( stencilIndex, weights, numPts );
  }
}


} /* namespace geosx */
