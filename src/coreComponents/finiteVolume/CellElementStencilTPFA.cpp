/*
 * CellElementStencil.cpp
 *
 *  Created on: Jul 10, 2019
 *      Author: settgast
 */

#include "CellElementStencilTPFA.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geosx
{

CellElementStencilTPFA::CellElementStencilTPFA():
  m_elementRegionIndices(),
  m_elementSubRegionIndices(),
  m_elementIndices(),
  m_weights(),
  m_connectorIndices()
{
  // TODO Auto-generated constructor stub

}

//CellElementStencilTPFA::~CellElementStencilTPFA()
//{
//  // TODO Auto-generated destructor stub
//}

void CellElementStencilTPFA::reserve( localIndex const size )
{
  m_elementRegionIndices.reserve( size * 2 );
  m_elementSubRegionIndices.reserve( size * 2 );
  m_elementIndices.reserve( size * 2 );
  m_weights.reserve( size * 2 );
}


void CellElementStencilTPFA::add( localIndex const numPts,
                                  localIndex const * const elementRegionIndices,
                                  localIndex const * const elementSubRegionIndices,
                                  localIndex const * const elementIndices,
                                  real64 const * const weights,
                                  localIndex const connectorIndex )
{
  GEOS_ERROR_IF( numPts!=2, "number of cells in TPFA stencil should be 2");

  localIndex const oldSize = m_elementRegionIndices.size(0);
  localIndex const newSize = oldSize + 1;
  m_elementRegionIndices.resize( newSize, numPts );
  m_elementSubRegionIndices.resize( newSize, numPts );
  m_elementIndices.resize( newSize, numPts );
  m_weights.resize( newSize, numPts );

  for( localIndex a=0 ; a<numPts ; ++a )
  {
    m_elementRegionIndices(oldSize,a) = elementRegionIndices[a];
    m_elementSubRegionIndices(oldSize,a) = elementSubRegionIndices[a];
    m_elementIndices(oldSize,a) = elementIndices[a];
    m_weights(oldSize,a) = weights[a];
  }
  m_connectorIndices[connectorIndex] = oldSize;
}

bool CellElementStencilTPFA::zero( localIndex const connectorIndex )
{
  return
  executeOnMapValue( m_connectorIndices, connectorIndex, [&]( localIndex const connectionListIndex )
  {
    for (localIndex i = 0; i < 2; ++i)
    {
      m_weights[connectionListIndex][i] = 0;
    }
  });
}
} /* namespace geosx */
