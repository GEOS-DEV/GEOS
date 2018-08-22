/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
#include "CellBlockSubRegion.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "dataRepository/ReferenceWrapper.hpp"
//#include "legacy/ArrayT/bufvector.h"
#include "ElementRegion.hpp"

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

  template< typename VIEWTYPE >
  using ElementViewAccessor = array1d< array1d< ReferenceWrapper< VIEWTYPE > > > ;

  /**
   * The MaterialViewAccessor at the ElementRegionManager level is an 3d array that contains a
   * ReferenceWrapper around the VIEWTYPE. The dimensions are denoted as follows:
   * var[elementRegionIndex][elementSubRegionIndex][materialIndexInRegion]
   */
  template< typename VIEWTYPE >
  using MaterialViewAccessor = array1d< array1d< array1d < ReferenceWrapper< VIEWTYPE > > > > ;

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

  void InitializePreSubGroups( ManagedGroup * const ) override final;
  void InitializePostSubGroups( ManagedGroup * const ) override final;

  // virtual void CreateChild( string const & childKey, string const & childName ) override;
  virtual void ReadXMLsub( xmlWrapper::xmlNode const & targetNode ) override;


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
  void forCellBlocks( LAMBDA lambda )
  {
    ManagedGroup * elementRegions = this->GetGroup(dataRepository::keys::elementRegions);

    for( auto & region : elementRegions->GetSubGroups() )
    {
      ManagedGroup * cellBlockSubRegions = region.second->GetGroup(dataRepository::keys::cellBlockSubRegions);
      for( auto & iterCellBlocks : cellBlockSubRegions->GetSubGroups() )
      {
        CellBlockSubRegion * cellBlock = cellBlockSubRegions->GetGroup<CellBlockSubRegion>(iterCellBlocks.first);
        lambda( cellBlock );
      }
    }
  }

  template< typename LAMBDA >
  void forCellBlocks( LAMBDA lambda ) const
  {
    ManagedGroup const * elementRegions = this->GetGroup(dataRepository::keys::elementRegions);

    for( auto const & region : elementRegions->GetSubGroups() )
    {
      ManagedGroup const * cellBlockSubRegions = region.second->GetGroup(dataRepository::keys::cellBlockSubRegions);
      for( auto const & iterCellBlocks : cellBlockSubRegions->GetSubGroups() )
      {
        CellBlockSubRegion const * cellBlock = cellBlockSubRegions->GetGroup<CellBlockSubRegion>(iterCellBlocks.first);
        lambda( cellBlock );
      }
    }
  }

  template< typename LAMBDA >
  void forCellBlocksComplete( LAMBDA lambda )
  {
    for( localIndex er=0 ; er<this->numRegions() ; ++er )
    {
      ElementRegion * elementRegion = this->GetRegion(er);

      for( localIndex esr=0 ;  esr<elementRegion->numSubRegions() ; ++esr )
      {
        CellBlockSubRegion * cellBlock = elementRegion->GetSubRegion(esr);
        lambda( er, esr, elementRegion, cellBlock );
      }
    }
  }

  template< typename LAMBDA >
  void forCellBlocksComplete( LAMBDA lambda ) const
  {
    for( localIndex er=0 ; er<this->numRegions() ; ++er )
    {
      ElementRegion const * elementRegion = this->GetRegion(er);

      for( localIndex esr=0 ;  esr<elementRegion->numSubRegions() ; ++esr )
      {
        CellBlockSubRegion const * cellBlock = elementRegion->GetSubRegion(esr);
        lambda( er, esr, elementRegion, cellBlock );
      }
    }
  }


  template< typename VIEWTYPE >
  ElementViewAccessor<VIEWTYPE> ConstructViewAccessor( string const & name,
                                                       string const & neighborName = string() );

  template< typename VIEWTYPE >
  ElementViewAccessor<VIEWTYPE const> ConstructViewAccessor( string const & name,
                                                             string const & neighborName = string() ) const;


  template< typename VIEWTYPE >
  MaterialViewAccessor<VIEWTYPE const>
  ConstructMaterialViewAccessor( string const & name,
                                 constitutive::ConstitutiveManager const * const cm ) const;

  template< typename VIEWTYPE >
  MaterialViewAccessor<VIEWTYPE>
  ConstructMaterialViewAccessor( string const & name,
                                 constitutive::ConstitutiveManager const * const cm  );

  template< typename CONSTITUTIVE_TYPE >
  ConstitutiveRelationAccessor<CONSTITUTIVE_TYPE>
  ConstructConstitutiveAccessor( constitutive::ConstitutiveManager const * const cm );

  using ManagedGroup::PackSize;
  using ManagedGroup::Pack;
  using ObjectManagerBase::PackGlobalMapsSize;
  using ObjectManagerBase::PackGlobalMaps;
  using ObjectManagerBase::UnpackGlobalMaps;
  using ObjectManagerBase::PackUpDownMapsSize;
  using ObjectManagerBase::PackUpDownMaps;
  using ObjectManagerBase::UnpackUpDownMaps;



  int PackSize( array1d<string> const & wrapperNames,
                ElementViewAccessor<localIndex_array> const & packList ) const;

  int Pack( buffer_unit_type * & buffer,
            array1d<string> const & wrapperNames,
            ElementViewAccessor<localIndex_array> const & packList ) const;

  using ObjectManagerBase::Unpack;
  int Unpack( buffer_unit_type const * & buffer,
              ElementViewAccessor<localIndex_array> & packList );



  int PackGlobalMapsSize( ElementViewAccessor<localIndex_array> const & packList ) const;

  int PackGlobalMaps( buffer_unit_type * & buffer,
                              ElementViewAccessor<localIndex_array> const & packList ) const;


  int UnpackGlobalMaps( buffer_unit_type const * & buffer,
                                ElementViewAccessor<localIndex_array> & packList );

  int PackUpDownMapsSize( ElementViewAccessor<localIndex_array> const & packList ) const;

  int PackUpDownMaps( buffer_unit_type * & buffer,
                      ElementViewAccessor<localIndex_array> const & packList ) const;


  int UnpackUpDownMaps( buffer_unit_type const * & buffer,
                        ElementViewAccessor<localIndex_array> const & packList );




private:
  template< bool DOPACK >
  int PackPrivate( buffer_unit_type * & buffer,
                   array1d<string> const & wrapperNames,
                   ElementViewAccessor<localIndex_array> const & viewAccessor ) const;

  template< bool DOPACK >
  int PackGlobalMapsPrivate( buffer_unit_type * & buffer,
                             ElementViewAccessor<localIndex_array> const & viewAccessor ) const;

  template< bool DOPACK >
  int
  PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                         ElementViewAccessor<localIndex_array> const & packList ) const;

  ElementRegionManager( const ElementRegionManager& );
  ElementRegionManager& operator=( const ElementRegionManager&);
};



template< typename VIEWTYPE >
ElementRegionManager::ElementViewAccessor<VIEWTYPE const>
ElementRegionManager::
ConstructViewAccessor( string const & viewName,
                       string const & neighborName ) const
{
  ElementViewAccessor<VIEWTYPE const> viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg=0 ; kReg<numRegions() ; ++kReg  )
  {
    ElementRegion const * const elemRegion = GetRegion(kReg);
    viewAccessor[kReg].resize( elemRegion->numSubRegions() );

    for( typename dataRepository::indexType kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion(kSubReg);

      if( neighborName.empty() )
      {
        //        viewAccessor[kReg].push_back( subRegion->getReference<VIEWTYPE>(viewName)) ;
        viewAccessor[kReg][kSubReg].set(subRegion->getReference<VIEWTYPE>(viewName));
      }
      else
      {
//        viewAccessor[kReg].push_back( subRegion->
        viewAccessor[kReg][kSubReg].set(subRegion->GetGroup(ObjectManagerBase::groupKeyStruct::neighborDataString)->
                                        GetGroup(neighborName)->getReference<VIEWTYPE>(viewName));
      }
    }
  }
  return viewAccessor;
}


template< typename VIEWTYPE >
ElementRegionManager::ElementViewAccessor<VIEWTYPE>
ElementRegionManager::
ConstructViewAccessor( string const & viewName,
                       string const & neighborName )
{
  ElementViewAccessor<VIEWTYPE> viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg=0 ; kReg<numRegions() ; ++kReg  )
  {
    ElementRegion * const elemRegion = GetRegion(kReg);
    viewAccessor[kReg].resize( elemRegion->numSubRegions() );

    for( typename dataRepository::indexType kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      CellBlockSubRegion * const subRegion = elemRegion->GetSubRegion(kSubReg);

      if( neighborName.empty() )
      {
        viewAccessor[kReg][kSubReg].set(subRegion->getReference<VIEWTYPE>(viewName));
      }
      else
      {
        viewAccessor[kReg][kSubReg].set(subRegion->GetGroup(ObjectManagerBase::groupKeyStruct::neighborDataString)->
                                        GetGroup(neighborName)->getReference<VIEWTYPE>(viewName));
      }
    }
  }
  return viewAccessor;
}




template< typename VIEWTYPE  >
ElementRegionManager::MaterialViewAccessor<VIEWTYPE const>
ElementRegionManager::
ConstructMaterialViewAccessor( string const & viewName,
                               constitutive::ConstitutiveManager const * const cm  ) const
{
  MaterialViewAccessor<VIEWTYPE const> accessor;
  accessor.resize( numRegions() );
  for( localIndex kReg=0 ; kReg<numRegions() ; ++kReg  )
  {
    ElementRegion const * const elemRegion = GetRegion(kReg);
    accessor[kReg].resize( elemRegion->numSubRegions() );

    for( localIndex kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion(kSubReg);
      dataRepository::ManagedGroup const * const
      constitutiveGroup = subRegion->GetConstitutiveModels();
      accessor[kReg][kSubReg].resize( cm->numSubGroups() );

      for( localIndex matIndex=0 ; matIndex<cm->numSubGroups() ; ++matIndex )
      {
        dataRepository::ManagedGroup const * const
        constitutiveRelation = constitutiveGroup->GetGroup(matIndex);
        if( constitutiveRelation != nullptr )
        {
          dataRepository::ViewWrapper<VIEWTYPE> * const
          wrapper = constitutiveRelation->getWrapper<VIEWTYPE>(viewName);

          if( wrapper != nullptr )
          {
            accessor[kReg][kSubReg][matIndex].set(wrapper->reference() );
          }
        }
      }
    }
  }
  return accessor;
}

/**
 * @tparam VIEWTYPE is the type that you would like to return inside the MaterialViewAccessor
 * @tparam REPOTYPE is the type that the data is stored as inside the repository
 * @param viewName The name of the material data field
 * @return a ElementRegionManager::MaterialViewAccessor<VIEWTYPE> that contains the data
 *
 * This function extracts data from the material arrays and provides multi-dimensional array-like
 * access to the data. The resulting return type is dereferenced using square bracket operators with
 * each level being, accessor[elementRegionIndex][elementSubRegionIndex][materialTypeIndex].
 *
 *
 */
template< typename VIEWTYPE >
ElementRegionManager::MaterialViewAccessor<VIEWTYPE>
ElementRegionManager::
ConstructMaterialViewAccessor( string const & viewName,
                               constitutive::ConstitutiveManager const * const cm  )
{
  MaterialViewAccessor<VIEWTYPE> accessor;
  accessor.resize( numRegions() );
  for( localIndex kReg=0 ; kReg<numRegions() ; ++kReg  )
  {
    ElementRegion * const elemRegion = GetRegion(kReg);
    accessor[kReg].resize( elemRegion->numSubRegions() );

    for( localIndex kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      CellBlockSubRegion * const subRegion = elemRegion->GetSubRegion(kSubReg);
      dataRepository::ManagedGroup * const
      constitutiveGroup = subRegion->GetConstitutiveModels();
      accessor[kReg][kSubReg].resize( cm->numSubGroups() );

      for( localIndex matIndex=0 ; matIndex<cm->numSubGroups() ; ++matIndex )
      {
        dataRepository::ManagedGroup * const
        constitutiveRelation = constitutiveGroup->GetGroup(matIndex);
        if( constitutiveRelation != nullptr )
        {
          dataRepository::ViewWrapper<VIEWTYPE> * const
          wrapper = constitutiveRelation->getWrapper<VIEWTYPE>(viewName);

          if( wrapper != nullptr )
          {
            accessor[kReg][kSubReg][matIndex].set(constitutiveRelation->getReference<VIEWTYPE>(viewName));
          }
        }
      }
    }
  }
  return accessor;
}

template< typename CONSTITUTIVE_TYPE >
ElementRegionManager::ConstitutiveRelationAccessor<CONSTITUTIVE_TYPE>
ElementRegionManager::ConstructConstitutiveAccessor( constitutive::ConstitutiveManager const * const cm )
{
  ConstitutiveRelationAccessor<CONSTITUTIVE_TYPE> accessor;
  accessor.resize( numRegions() );
  for( localIndex kReg=0 ; kReg<numRegions() ; ++kReg  )
  {
    ElementRegion * const elemRegion = GetRegion(kReg);
    accessor[kReg].resize( elemRegion->numSubRegions() );

    for( localIndex kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      CellBlockSubRegion * const subRegion = elemRegion->GetSubRegion(kSubReg);
      dataRepository::ManagedGroup * const
      constitutiveGroup = subRegion->GetConstitutiveModels();
      accessor[kReg][kSubReg].resize( cm->numSubGroups() );

      for( localIndex matIndex=0 ; matIndex<cm->numSubGroups() ; ++matIndex )
      {
        CONSTITUTIVE_TYPE * const
        constitutiveRelation = constitutiveGroup->GetGroup<CONSTITUTIVE_TYPE>(matIndex);
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
