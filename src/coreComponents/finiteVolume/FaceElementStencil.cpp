/*
 * FaceElementStencil.cpp
 *
 *  Created on: Jul 1, 2019
 *      Author: settgast
 */

#include "FaceElementStencil.hpp"

namespace geosx
{

FaceElementStencil::FaceElementStencil():
  m_elementRegionIndices(),
  m_elementSubRegionIndices(),
  m_elementIndices(),
  m_weights(),
  m_connectorIndices()
{}

//FaceElementStencil::~FaceElementStencil()
//{}

void FaceElementStencil::add( localIndex const numPts,
                              localIndex  const * const elementRegionIndices,
                              localIndex  const * const elementSubRegionIndices,
                              localIndex  const * const elementIndices,
                              real64 const * const weights,
                              localIndex const connectorIndex )
{
  GEOS_ERROR_IF( numPts >= MAX_STENCIL_SIZE, "Maximum stencil size exceeded" );

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
    m_elementRegionIndices.clearArray( iter->second );
    m_elementSubRegionIndices.clearArray( iter->second );
    m_elementIndices.clearArray( iter->second );
    m_weights.clearArray( iter->second );

    m_elementRegionIndices.appendToArray( iter->second, elementRegionIndices, numPts );
    m_elementSubRegionIndices.appendToArray( iter->second, elementSubRegionIndices, numPts );
    m_elementIndices.appendToArray( iter->second, elementIndices, numPts );
    m_weights.appendToArray( iter->second, weights, numPts );
  }


}


} /* namespace geosx */
