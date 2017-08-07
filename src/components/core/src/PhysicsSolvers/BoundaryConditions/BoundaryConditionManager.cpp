/*
 * BoundaryConditionManager.cpp
 *
 *  Created on: May 26, 2017
 *      Author: rrsettgast
 */

#include "BoundaryConditionManager.hpp"
#include "BoundaryConditionBase.hpp"

namespace geosx
{
using namespace dataRepository;
BoundaryConditionManager::BoundaryConditionManager( string const & name, ManagedGroup * const parent ):
ManagedGroup(name,parent)
{
  // TODO Auto-generated constructor stub

}


BoundaryConditionManager & BoundaryConditionManager::get()
{
  static BoundaryConditionManager bcman("bcMan",nullptr);
  return bcman;
}


BoundaryConditionManager::~BoundaryConditionManager()
{
  // TODO Auto-generated destructor stub
}

void BoundaryConditionManager::ReadXMLsub( xmlWrapper::xmlNode const & targetNode )
{
  for (xmlWrapper::xmlNode childNode=targetNode.first_child(); childNode; childNode=childNode.next_sibling())
  {
    string const typeName = childNode.name();
    string const name = childNode.attribute("name").value();
    std::unique_ptr<BoundaryConditionBase> bc = BoundaryConditionBase::CatalogInterface::Factory( typeName, name, this );
    this->RegisterGroup(name, std::move(bc) );
  }
}

} /* namespace geosx */
