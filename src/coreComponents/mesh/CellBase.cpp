/*
 * CellBase.cpp
 *
 *  Created on: Jan 14, 2019
 *      Author: settgast
 */

#include "CellBase.hpp"

namespace geosx
{

CellBase::CellBase( string const & name, ManagedGroup * const parent ):
  ObjectManagerBase(name,parent),
  m_numNodesPerElement(),
  m_numEdgesPerElement(),
  m_numFacesPerElement(),
  m_elementCenter(),
  m_elementVolume()
{
}

CellBase::~CellBase()
{
}

} /* namespace geosx */
