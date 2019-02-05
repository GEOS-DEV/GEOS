/*
 * CellBase.cpp
 *
 *  Created on: Jan 14, 2019
 *      Author: settgast
 */

#include "ElementSubRegionBase.hpp"

namespace geosx
{
using namespace dataRepository;

ElementSubRegionBase::ElementSubRegionBase( string const & name, ManagedGroup * const parent ):
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

ElementSubRegionBase::~ElementSubRegionBase()
{
}




} /* namespace geosx */
