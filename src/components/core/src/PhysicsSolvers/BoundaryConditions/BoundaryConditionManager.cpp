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
    bc->SetDocumentationNodes(nullptr);
    bc->ReadXML(childNode);
    this->RegisterGroup(name, std::move(bc) );
  }
}

void BoundaryConditionManager::ApplyBoundaryCondition( dataRepository::ManagedGroup & object,
                                                       std::string const & fieldName,
                                                       real64 const time )
{
  dataRepository::ManagedGroup const & sets = object.GetGroup(dataRepository::keys::sets);

  // iterate over all boundary conditions.
  forSubGroups<BoundaryConditionBase>( [&]( BoundaryConditionBase & bc ) -> void
  {
    if( time >= bc.GetStartTime() && time < bc.GetEndTime() && ( bc.GetFieldName()==fieldName) )
    {
      string_array setNames = bc.GetSetNames();
      for( auto & setName : setNames )
      {
        dataRepository::ViewWrapper<lSet> const * const setWrapper = sets.getWrapperPtr<lSet>(setName);
        if( setWrapper != nullptr )
        {
          lSet const & set = setWrapper->reference();
          bc.ApplyBounaryConditionDefaultMethod(set,time, object, fieldName);
        }
      }
    }
  });
}
} /* namespace geosx */
