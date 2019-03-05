/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file ElementManagerT.h
 * @author Randolph Settgast
 * @date created on Sep 14, 2010
 */

#ifndef ZONEMANAGER_H
#define ZONEMANAGER_H

//#include "Common.h"
//#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "CellBlock.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "CellElementSubRegion.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "dataRepository/ReferenceWrapper.hpp"
//#include "legacy/ArrayT/bufvector.h"
#include "ElementRegion.hpp"
#include "fileIO/schema/SchemaUtilities.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const elementRegions = "elementRegions";
//string const elementRegionManager="ElementRegions";
}
}

/**
 * Class to manage the data stored at the element level.
 */
class ElementRegionManager : public ObjectManagerBase
{
public:

  constexpr static int maxNumNodesPerElem = 8;

  template< typename VIEWTYPE >
  using ElementViewAccessor = array1d< array1d< VIEWTYPE > > ;

  template< typename VIEWTYPE >
  using ElementReferenceAccessor = array1d< array1d< ReferenceWrapper<VIEWTYPE> > > ;

  /**
   * The MaterialViewAccessor at the ElementRegionManager level is an 3d array that contains a
   * ReferenceWrapper around the VIEWTYPE. The dimensions are denoted as follows:
   * var[elementRegionIndex][elementSubRegionIndex][materialIndexInRegion]
   */
  template< typename VIEWTYPE >
  using MaterialViewAccessor = array1d< array1d< array1d < VIEWTYPE > > > ;

  template< typename CONSTITUTIVE_TYPE >
  using ConstitutiveRelationAccessor = array1d< array1d< array1d< CONSTITUTIVE_TYPE * > > >;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  static const string CatalogName()
  { return "ZoneManager"; }

  virtual const string getCatalogName() const override final
  { return ElementRegionManager::CatalogName(); }



  ///@}

  ElementRegionManager( string const &, ManagedGroup * const parent );
  virtual ~ElementRegionManager() override;

  localIndex getNumberOfElements() const;

//  void Initialize(  ){}

  void GenerateMesh( ManagedGroup const * const cellBlockManager );

  void GenerateFractureMesh( FaceManager const * const faceManager );

  virtual ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;
//  virtual void ReadXMLsub( xmlWrapper::xmlNode const & targetNode ) override;


  /**
   * This function is used to expand any objects in the data structure.
   * Currently, there is only one type of element region
   */
  virtual void ExpandObjectCatalogs() override;

  /**
   * This function is used to inform the schema generator of any
   * deviations between the xml and GEOS data structures.
   */
  virtual void SetSchemaDeviations(xmlWrapper::xmlNode schemaRoot,
                                   xmlWrapper::xmlNode schemaParent,
                                   integer documentationType) override;

  using ManagedGroup::resize;

  void resize( integer_array const & numElements,
               string_array const & regionNames,
               string_array const & elementTypes );

//  CellBlock & CreateRegion( string const & regionName,
//                               string const & elementType,
//                               integer const & numElements );


  subGroupMap const & GetRegions() const
  {
    return this->GetGroup(dataRepository::keys::elementRegions)->GetSubGroups();
  }
  subGroupMap & GetRegions()
  {
    return this->GetGroup(dataRepository::keys::elementRegions)->GetSubGroups();
  }

  ElementRegion const * GetRegion( string const & regionName ) const
  {
    return this->GetGroup(dataRepository::keys::elementRegions)->GetGroup<ElementRegion>(regionName);
  }
  ElementRegion * GetRegion( string const & regionName )
  {
    return this->GetGroup(dataRepository::keys::elementRegions)->GetGroup<ElementRegion>(regionName);
  }

  ElementRegion const * GetRegion( localIndex const & index ) const
  {
    return this->GetGroup(dataRepository::keys::elementRegions)->GetGroup<ElementRegion>(index);
  }
  ElementRegion * GetRegion( localIndex const & index )
  {
    return this->GetGroup(dataRepository::keys::elementRegions)->GetGroup<ElementRegion>(index);
  }

  localIndex numRegions() const
  {
    return this->GetGroup(dataRepository::keys::elementRegions)->GetSubGroups().size();
  }

  localIndex numCellBlocks() const;


  template< typename LAMBDA >
  void forElementRegions( LAMBDA lambda )
  {
    ManagedGroup * elementRegions = this->GetGroup(dataRepository::keys::elementRegions);
    elementRegions->forSubGroups<ElementRegion>( lambda );
  }

  template< typename LAMBDA >
  void forElementRegions( LAMBDA lambda ) const
  {
    ManagedGroup const * elementRegions = this->GetGroup(dataRepository::keys::elementRegions);
    elementRegions->forSubGroups<ElementRegion>( lambda );
  }


  template< typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda )
  {
    forElementSubRegions<CellElementSubRegion,FaceElementSubRegion>( std::forward<LAMBDA>(lambda) );
  }


  template< typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda ) const
  {
    forElementSubRegions<CellElementSubRegion,FaceElementSubRegion>( std::forward<LAMBDA>(lambda) );
  }

  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda )
  {
    ManagedGroup * elementRegions = this->GetGroup(dataRepository::keys::elementRegions);

    for( auto & region : elementRegions->GetSubGroups() )
    {
      ElementRegion * const elemRegion = region.second->group_cast<ElementRegion *>();
      elemRegion->forElementSubRegions<SUBREGIONTYPE,SUBREGIONTYPES...>( std::forward<LAMBDA>(lambda) );
    }
  }

  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda ) const
  {
    ManagedGroup const * elementRegions = this->GetGroup(dataRepository::keys::elementRegions);

    for( auto & region : elementRegions->GetSubGroups() )
    {
      ElementRegion const * const elemRegion = region.second->group_cast<ElementRegion const *>();
      elemRegion->forElementSubRegions<SUBREGIONTYPE,SUBREGIONTYPES...>( std::forward<LAMBDA>(lambda) );
    }
  }


  template< typename LAMBDA >
  void forElementSubRegionsComplete( LAMBDA lambda ) const
  {
    forElementSubRegionsComplete<CellElementSubRegion,FaceElementSubRegion>( std::forward<LAMBDA>(lambda) );
  }
  template< typename LAMBDA >
  void forElementSubRegionsComplete( LAMBDA lambda )
  {
    forElementSubRegionsComplete<CellElementSubRegion,FaceElementSubRegion>( std::forward<LAMBDA>(lambda) );
  }


  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegionsComplete( LAMBDA lambda )
  {
    for( localIndex er=0 ; er<this->numRegions() ; ++er )
    {
      ElementRegion * const elementRegion = this->GetRegion(er);

      for( localIndex esr=0 ;  esr<elementRegion->numSubRegions() ; ++esr )
      {
        ElementSubRegionBase * const subRegion = elementRegion->GetSubRegion(esr);

        bool validCast =
        ElementRegion::applyLambdaToCellBlocks<SUBREGIONTYPE,SUBREGIONTYPES...>( subRegion, [&]( auto * const castedSubRegion )
        {
          lambda( er, esr, elementRegion, castedSubRegion );
        });
      }
    }
  }

  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegionsComplete( LAMBDA lambda ) const
  {
    for( localIndex er=0 ; er<this->numRegions() ; ++er )
    {
      ElementRegion const * const elementRegion = this->GetRegion(er);

      for( localIndex esr=0 ;  esr<elementRegion->numSubRegions() ; ++esr )
      {
        ElementSubRegionBase const * const subRegion = elementRegion->GetSubRegion(esr);

        ElementRegion::applyLambdaToCellBlocks<SUBREGIONTYPE,SUBREGIONTYPES...>( subRegion, [&]( auto const * const castedSubRegion )
        {
          lambda( er, esr, elementRegion, castedSubRegion );
        });
      }
    }
  }


  // template< typename VIEWTYPE >
  // ElementViewAccessor<VIEWTYPE> ConstructViewAccessor( string const & name,
  //                                                      string const & neighborName = string() );

  template< typename VIEWTYPE, typename LHS=VIEWTYPE >
  ElementViewAccessor< LHS > ConstructViewAccessor( string const & name,
                                                    string const & neighborName = string() ) const;
  template< typename VIEWTYPE, typename LHS=VIEWTYPE >
  ElementViewAccessor< LHS > ConstructViewAccessor( string const & name,
                                                    string const & neighborName = string() );

  template< typename VIEWTYPE >
  ElementViewAccessor< ReferenceWrapper<VIEWTYPE> >
  ConstructReferenceAccessor( string const & viewName, string const & neighborName = string() ) const;

  template< typename VIEWTYPE >
  ElementViewAccessor< ReferenceWrapper<VIEWTYPE> >
  ConstructReferenceAccessor( string const & viewName, string const & neighborName = string() );

  template< typename VIEWTYPE, typename LHS=VIEWTYPE >
  MaterialViewAccessor< LHS >
  ConstructFullMaterialViewAccessor( string const & name,
                                 constitutive::ConstitutiveManager const * const cm ) const;

//  template< typename VIEWTYPE, typename LHS=VIEWTYPE >
//  MaterialViewAccessor< LHS >
//  ConstructFullMaterialViewAccessor( string const & name,
//                                     constitutive::ConstitutiveManager const * const cm ) const;


  template< typename CONSTITUTIVE_TYPE >
  ConstitutiveRelationAccessor< CONSTITUTIVE_TYPE >
  ConstructFullConstitutiveAccessor( constitutive::ConstitutiveManager const * const cm );

  using ManagedGroup::PackSize;
  using ManagedGroup::Pack;
  using ObjectManagerBase::PackGlobalMapsSize;
  using ObjectManagerBase::PackGlobalMaps;
  using ObjectManagerBase::UnpackGlobalMaps;
  using ObjectManagerBase::PackUpDownMapsSize;
  using ObjectManagerBase::PackUpDownMaps;
  using ObjectManagerBase::UnpackUpDownMaps;



  int PackSize( string_array const & wrapperNames,
                ElementViewAccessor<arrayView1d<localIndex>> const & packList ) const;

  int Pack( buffer_unit_type * & buffer,
            string_array const & wrapperNames,
            ElementViewAccessor<arrayView1d<localIndex>> const & packList ) const;

  using ObjectManagerBase::Unpack;
  int Unpack( buffer_unit_type const * & buffer,
              ElementViewAccessor<arrayView1d<localIndex>> & packList );

  int Unpack( buffer_unit_type const * & buffer,
              ElementReferenceAccessor<array1d<localIndex>> & packList );



  int PackGlobalMapsSize( ElementViewAccessor<arrayView1d<localIndex>> const & packList ) const;

  int PackGlobalMaps( buffer_unit_type * & buffer,
                              ElementViewAccessor<arrayView1d<localIndex>> const & packList ) const;


  int UnpackGlobalMaps( buffer_unit_type const * & buffer,
                                ElementViewAccessor<ReferenceWrapper<localIndex_array>> & packList );

  int PackUpDownMapsSize( ElementViewAccessor<arrayView1d<localIndex>> const & packList ) const;
  int PackUpDownMapsSize( ElementReferenceAccessor<array1d<localIndex>> const & packList ) const;

  int PackUpDownMaps( buffer_unit_type * & buffer,
                      ElementViewAccessor<arrayView1d<localIndex>> const & packList ) const;
  int PackUpDownMaps( buffer_unit_type * & buffer,
                      ElementReferenceAccessor<array1d<localIndex>> const & packList ) const;


  int UnpackUpDownMaps( buffer_unit_type const * & buffer,
                        ElementReferenceAccessor<localIndex_array> & packList,
                        bool const overwriteMap );




private:
  template< bool DOPACK >
  int PackPrivate( buffer_unit_type * & buffer,
                   string_array const & wrapperNames,
                   ElementViewAccessor<arrayView1d<localIndex>> const & viewAccessor ) const;

  template< bool DOPACK >
  int PackGlobalMapsPrivate( buffer_unit_type * & buffer,
                             ElementViewAccessor<arrayView1d<localIndex>> const & viewAccessor ) const;

  template< bool DOPACK, typename T >
  int
  PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                         T const & packList ) const;

  template< typename T >
  int UnpackPrivate( buffer_unit_type const * & buffer,
                     T & packList );

  ElementRegionManager( const ElementRegionManager& );
  ElementRegionManager& operator=( const ElementRegionManager&);
};


template< typename VIEWTYPE, typename LHS >
ElementRegionManager::ElementViewAccessor<LHS>
ElementRegionManager::ConstructViewAccessor( string const & viewName, string const & neighborName ) const
{
  ElementViewAccessor<LHS> viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg=0 ; kReg<numRegions() ; ++kReg  )
  {
    ElementRegion const * const elemRegion = GetRegion(kReg);
    viewAccessor[kReg].resize( elemRegion->numSubRegions() );

    for( typename dataRepository::indexType kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      ManagedGroup const * group = elemRegion->GetSubRegion(kSubReg);

      if( !neighborName.empty() )
      {
        group = group->GetGroup(ObjectManagerBase::groupKeyStruct::neighborDataString)->GetGroup(neighborName);
      }

      if ( group->hasView(viewName) )
      {
        viewAccessor[kReg][kSubReg] = group->getReference<VIEWTYPE>(viewName);
      }
    }
  }
  return viewAccessor;
}


template< typename VIEWTYPE, typename LHS >
ElementRegionManager::ElementViewAccessor<LHS>
ElementRegionManager::
ConstructViewAccessor( string const & viewName, string const & neighborName )
{
  ElementViewAccessor<LHS> viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg=0 ; kReg<numRegions() ; ++kReg  )
  {
    ElementRegion * const elemRegion = GetRegion(kReg);
    viewAccessor[kReg].resize( elemRegion->numSubRegions() );

    for( typename dataRepository::indexType kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      ManagedGroup * group = elemRegion->GetSubRegion(kSubReg);

      if( !neighborName.empty() )
      {
        group = group->GetGroup(ObjectManagerBase::groupKeyStruct::neighborDataString)->GetGroup(neighborName);
      }

      if ( group->hasView(viewName) )
      {
        viewAccessor[kReg][kSubReg] = group->getReference<VIEWTYPE>(viewName);
      }
    }
  }
  return viewAccessor;
}

template< typename VIEWTYPE >
ElementRegionManager::ElementViewAccessor<ReferenceWrapper<VIEWTYPE>>
ElementRegionManager::
ConstructReferenceAccessor( string const & viewName, string const & neighborName ) const
{
  ElementViewAccessor<ReferenceWrapper<VIEWTYPE>> viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg=0 ; kReg<numRegions() ; ++kReg  )
  {
    ElementRegion const * const elemRegion = GetRegion(kReg);
    viewAccessor[kReg].resize( elemRegion->numSubRegions() );

    for( typename dataRepository::indexType kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      ManagedGroup const * group = elemRegion->GetSubRegion(kSubReg);

      if( !neighborName.empty() )
      {
        group = group->GetGroup(ObjectManagerBase::groupKeyStruct::neighborDataString)->GetGroup(neighborName);
      }

      if ( group->hasView(viewName) )
      {
        viewAccessor[kReg][kSubReg].set(group->getReference<VIEWTYPE>(viewName));
      }
    }
  }
  return viewAccessor;
}

template< typename VIEWTYPE >
ElementRegionManager::ElementViewAccessor<ReferenceWrapper<VIEWTYPE>>
ElementRegionManager::
ConstructReferenceAccessor( string const & viewName, string const & neighborName )
{
  ElementViewAccessor<ReferenceWrapper<VIEWTYPE>> viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg=0 ; kReg<numRegions() ; ++kReg  )
  {
    ElementRegion * const elemRegion = GetRegion(kReg);
    viewAccessor[kReg].resize( elemRegion->numSubRegions() );

    for( typename dataRepository::indexType kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      ManagedGroup * group = elemRegion->GetSubRegion(kSubReg);

      if( !neighborName.empty() )
      {
        group = group->GetGroup(ObjectManagerBase::groupKeyStruct::neighborDataString)->GetGroup(neighborName);
      }

      if ( group->hasView(viewName) )
      {
        viewAccessor[kReg][kSubReg].set(group->getReference<VIEWTYPE>(viewName));
      }
    }
  }
  return viewAccessor;
}

//template< typename VIEWTYPE, typename LHS >
//ElementRegionManager::MaterialViewAccessor<LHS>
//ElementRegionManager::
//ConstructMaterialViewAccessor( string const & viewName,
//                               constitutive::ConstitutiveManager const * const cm  ) const
//{
//  MaterialViewAccessor<LHS> accessor;
//  accessor.resize( numRegions() );
//  for( localIndex kReg=0 ; kReg<numRegions() ; ++kReg  )
//  {
//    ElementRegion const * const elemRegion = GetRegion(kReg);
//    accessor[kReg].resize( elemRegion->numSubRegions() );
//
//    for( localIndex kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
//    {
//      ElementSubRegionBase const * const subRegion = elemRegion->GetSubRegion(kSubReg);
//      dataRepository::ManagedGroup const * const constitutiveGroup = subRegion->GetConstitutiveModels();
//      accessor[kReg][kSubReg].resize( constitutiveGroup->numSubGroups() );
//
//      for( localIndex matIndex=0 ; matIndex<constitutiveGroup->numSubGroups() ; ++matIndex )
//      {
//        dataRepository::ManagedGroup const * const
//        constitutiveRelation = constitutiveGroup->GetGroup(matIndex);
//        if( constitutiveRelation != nullptr )
//        {
//          dataRepository::ViewWrapper<VIEWTYPE> const * const
//          wrapper = constitutiveRelation->getWrapper<VIEWTYPE>(viewName);
//
//          if( wrapper != nullptr )
//          {
//            accessor[kReg][kSubReg][matIndex] = wrapper->reference();
//          }
//        }
//      }
//    }
//  }
//  return accessor;
//}


template< typename VIEWTYPE, typename LHS >
ElementRegionManager::MaterialViewAccessor<LHS>
ElementRegionManager::
ConstructFullMaterialViewAccessor( string const & viewName,
                                   constitutive::ConstitutiveManager const * const cm  ) const
{
  MaterialViewAccessor<LHS> accessor;
  accessor.resize( numRegions() );
  for( localIndex kReg=0 ; kReg<numRegions() ; ++kReg  )
  {
    ElementRegion const * const elemRegion = GetRegion(kReg);
    accessor[kReg].resize( elemRegion->numSubRegions() );

    for( localIndex kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      ElementSubRegionBase const * const subRegion = elemRegion->GetSubRegion(kSubReg);
      dataRepository::ManagedGroup const * const constitutiveGroup = subRegion->GetConstitutiveModels();

      accessor[kReg][kSubReg].resize( cm->numSubGroups() );

      for( localIndex matIndex=0 ; matIndex<cm->numSubGroups() ; ++matIndex )
      {
        string constitutiveName = cm->GetGroup(matIndex)->getName();
        dataRepository::ManagedGroup const * const constitutiveRelation = constitutiveGroup->GetGroup(constitutiveName);
        if( constitutiveRelation != nullptr )
        {
          dataRepository::ViewWrapper<VIEWTYPE> const * const
          wrapper = constitutiveRelation->getWrapper<VIEWTYPE>(viewName);

          if( wrapper != nullptr )
          {
            accessor[kReg][kSubReg][matIndex] = wrapper->reference();
          }
        }
      }
    }
  }
  return accessor;
}

template< typename CONSTITUTIVE_TYPE >
ElementRegionManager::ConstitutiveRelationAccessor<CONSTITUTIVE_TYPE>
ElementRegionManager::ConstructFullConstitutiveAccessor( constitutive::ConstitutiveManager const * const cm )
{
  ConstitutiveRelationAccessor<CONSTITUTIVE_TYPE> accessor;
  accessor.resize( numRegions() );
  for( localIndex kReg=0 ; kReg<numRegions() ; ++kReg  )
  {
    ElementRegion * const elemRegion = GetRegion(kReg);
    accessor[kReg].resize( elemRegion->numSubRegions() );

    for( localIndex kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      ElementSubRegionBase * const subRegion = elemRegion->GetSubRegion(kSubReg);
      dataRepository::ManagedGroup * const
      constitutiveGroup = subRegion->GetConstitutiveModels();
      accessor[kReg][kSubReg].resize( cm->numSubGroups() );

      for( localIndex matIndex=0 ; matIndex<cm->numSubGroups() ; ++matIndex )
      {
        string const constitutiveName = cm->GetGroup(matIndex)->getName();

        CONSTITUTIVE_TYPE * const
        constitutiveRelation = constitutiveGroup->GetGroup<CONSTITUTIVE_TYPE>(constitutiveName);
        if( constitutiveRelation != nullptr )
        {
          accessor[kReg][kSubReg][matIndex] = constitutiveRelation;
        }
      }
    }
  }
  return accessor;
}

}
#endif /* ZONEMANAGER_H */
