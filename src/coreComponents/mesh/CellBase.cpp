/*
 * CellBase.cpp
 *
 *  Created on: Jan 14, 2019
 *      Author: settgast
 */

#include "CellBase.hpp"

namespace geosx
{
using namespace dataRepository;

CellBase::CellBase( string const & name, ManagedGroup * const parent ):
  ObjectManagerBase(name,parent),
  m_constitutiveModels(groupKeyStruct::constitutiveModelsString,this),
  m_numNodesPerElement(),
  m_numEdgesPerElement(),
  m_numFacesPerElement(),
  m_elementCenter(),
  m_elementVolume()
{
  RegisterGroup( groupKeyStruct::constitutiveModelsString, &m_constitutiveModels, 0 );
}

CellBase::~CellBase()
{
}




} /* namespace geosx */
