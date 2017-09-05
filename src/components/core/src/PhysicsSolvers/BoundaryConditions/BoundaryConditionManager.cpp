/*
 * BoundaryConditionManager.cpp
 *
 *  Created on: May 26, 2017
 *      Author: rrsettgast
 */

#include "BoundaryConditionManager.hpp"
#include "BoundaryConditionBase.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "managers/CellBlockSubRegion.hpp"

#include "managers/ElementRegionManager.hpp"
#include "managers/ElementRegion.hpp"
#include "finiteElement/FiniteElementManager.hpp"
#include "finiteElement/ElementLibrary/FiniteElement.h"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;
BoundaryConditionManager::BoundaryConditionManager( string const & name, ManagedGroup * const parent ):
ManagedGroup(name,parent)
{
  // TODO Auto-generated constructor stub

}


BoundaryConditionManager * BoundaryConditionManager::get()
{
  static BoundaryConditionManager bcman("bcMan",nullptr);
  return &bcman;
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

void BoundaryConditionManager::ApplyBoundaryCondition( dataRepository::ManagedGroup * object,
                                                       std::string const & fieldName,
                                                       real64 const time )
{
  dataRepository::ManagedGroup const * sets = object->GetGroup(dataRepository::keys::sets);

  // iterate over all boundary conditions.
  forSubGroups<BoundaryConditionBase>( [&]( BoundaryConditionBase * bc ) -> void
  {
    if( time >= bc->GetStartTime() && time < bc->GetEndTime() && ( bc->GetFieldName()==fieldName) )
    {
      string_array setNames = bc->GetSetNames();
      for( auto & setName : setNames )
      {
        dataRepository::ViewWrapper<lSet> const * const setWrapper = sets->getWrapper<lSet>(setName);
        if( setWrapper != nullptr )
        {
          lSet const & set = setWrapper->reference();
          bc->ApplyBounaryConditionDefaultMethod(set,time, object, fieldName);
        }
      }
    }
  });
}

void BoundaryConditionManager::ApplyInitialConditions( ManagedGroup * domain ) const
{

  forSubGroups<BoundaryConditionBase>( [&] ( BoundaryConditionBase const * bc )-> void
  {
    if( bc->initialCondition() )
    {

      if( bc->GetElementRegion().empty() )
      {
        string_array objectPath = stringutilities::Tokenize( bc->GetFieldName(), "/");
        int32 const pathLength = objectPath.size();
        ManagedGroup * currentGroup = domain;
        for( int32 a=0 ; a<(pathLength-1) ; ++a )
        {
          currentGroup = currentGroup->GetGroup(objectPath[a]);
        }
        string const fieldName = objectPath[pathLength-1];

        dataRepository::ManagedGroup const * setGroup = currentGroup->GetGroup(dataRepository::keys::sets);
        string_array setNames = bc->GetSetNames();
        for( auto & setName : setNames )
        {
          dataRepository::ViewWrapper<lSet> const * const setWrapper = setGroup->getWrapper<lSet>(setName);
          if( setWrapper != nullptr )
          {
            lSet const & set = setWrapper->reference();
            bc->ApplyBounaryConditionDefaultMethod( set, 0.0, currentGroup, fieldName );
          }
        }
      }
      else
      {
        ConstitutiveManager * constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);
        ConstitutiveManager::constitutiveMaps const & constitutiveMaps = constitutiveManager->GetMaps(0);

        lSet targetSet;



        string_array targetPath = stringutilities::Tokenize( bc->GetFieldName(), "/");
        int32 const targetPathLength = targetPath.size();
        ManagedGroup * targetGroup = domain;
        for( int32 a=0 ; a<(targetPathLength-1) ; ++a )
        {
          targetGroup = targetGroup->GetGroup(targetPath[a]);
        }
        string const materialName = targetPath[targetPathLength-3];
        string const fieldName = targetPath[targetPathLength-1];


        FiniteElementManager const * numericalMethodManager = domain->getParent()->GetGroup<FiniteElementManager>(keys::finiteElementManager);

        // Get element Region
        string const elementRegionName = bc->GetElementRegion();
        ManagedGroup * ElementRegionManager = domain->GetGroup(keys::FEM_Elements);
        ManagedGroup * ElementRegions = ElementRegionManager->GetGroup(keys::elementRegions);
        ElementRegion * elementRegion = ElementRegions->GetGroup<ElementRegion>(elementRegionName);
        // ManagedGroup * elementSubRegions = elementRegion->GetGroup(dataRepository::keys::cellBlockSubRegions);

        auto const & numMethodName = elementRegion->getData<string>(keys::numericalMethod);
        FiniteElementSpace const * feSpace = numericalMethodManager->GetGroup<FiniteElementSpace>(numMethodName);

        string_array setNames = bc->GetSetNames();

        elementRegion->forCellBlocks( [&] ( CellBlockSubRegion * subRegion ) -> void
        {
          view_rtype_const< Array2dT<mapPair> > constitutiveMap = subRegion->getData< Array2dT<mapPair> >(keys::constitutiveMap);
          ManagedGroup const * sets = subRegion->GetGroup(keys::sets);

          for( auto & setName : setNames )
          {
            dataRepository::ViewWrapper<lSet> const * const setWrapper = sets->getWrapper<lSet>(setName);
            if( setWrapper != nullptr )
            {
              lSet const & set = setWrapper->reference();
              int32 materialIndex = constitutiveMaps.second.at(materialName);
              for( auto const & k : set )
              {
                if( constitutiveMap[k]->first == materialIndex )
                {
                  for( auto q=0 ; q < feSpace->m_finiteElement->n_quadrature_points() ; ++q )
                  {
                    targetSet.insert(constitutiveMap[k][q].second);
                  }
                }
              }
            }
          }

        });

        // ManagedGroup const * elementSets = elementRegion->GetGroup(dataRepository::keys::sets);


        bc->ApplyBounaryConditionDefaultMethod( targetSet, 0.0, targetGroup, fieldName );
      }
    }
  });
}

} /* namespace geosx */
