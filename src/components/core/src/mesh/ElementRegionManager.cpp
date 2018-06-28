// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/*
 * ElementManagerT.cpp
 *
 *  Created on: Sep 14, 2010
 *      Author: settgast1
 */

#include "ElementRegion.hpp"
#include "ElementRegionManager.hpp"
#include "FaceManager.hpp"
//#include "legacy/IO/BinStream.h"
#include <map>
#include <vector>
//#include "legacy/Constitutive/Material/MaterialFactory.h"
//#include "legacy/ArrayT/ArrayT.h"

namespace geosx
{
using namespace dataRepository;

ElementRegionManager::ElementRegionManager(  string const & name, ManagedGroup * const parent ):
  ObjectManagerBase(name,parent)
{
  this->RegisterGroup<ManagedGroup>(keys::elementRegions);
}

ElementRegionManager::~ElementRegionManager()
{
  // TODO Auto-generated destructor stub
}

localIndex ElementRegionManager::getNumberOfElements() const
{
  localIndex numElem = 0;
  this->forCellBlocks([&]( ManagedGroup const * cellBlock ) -> void
    {
      numElem += cellBlock->size();
    });
  return numElem;
}

localIndex ElementRegionManager::numCellBlocks() const
{
  localIndex numCellBlocks = 0;
  this->forCellBlocks([&]( ManagedGroup const * cellBlock ) -> void
    {
    numCellBlocks += 1;
    });
  return numCellBlocks;
}

void ElementRegionManager::resize( integer_array const & numElements,
                                   string_array const & regionNames,
                                   string_array const & elementTypes )
{
  localIndex const n_regions = regionNames.size();
//  ManagedGroup * elementRegions = this->GetGroup(keys::cellBlocks);
  for( localIndex reg=0 ; reg<n_regions ; ++reg )
  {
    ElementRegion * elemRegion = this->GetRegion( regionNames[reg] );
    elemRegion->resize(numElements[reg]);
  }
}


//CellBlock & ZoneManager::CreateRegion( string const & regionName,
//                                             string const & elementType,
//                                             integer const & numElements )
//{
////  ElementRegion * elemRegion = elementRegions.RegisterGroup( regionNames );
////  elemRegion->resize(numElements);
//}

// void ElementRegionManager::CreateChild( string const & childKey, string const & childName )
// {
// }


void ElementRegionManager::ReadXMLsub( xmlWrapper::xmlNode const & targetNode )
{
  ManagedGroup * elementRegions = this->GetGroup(keys::elementRegions);
  for (xmlWrapper::xmlNode childNode=targetNode.first_child() ; childNode ; childNode=childNode.next_sibling())
  {
    if( childNode.name() == string("ElementRegion") )
    {
      std::string regionName = childNode.attribute("name").value();
      std::cout<<regionName<<std::endl;

      ElementRegion * elemRegion = elementRegions->RegisterGroup<ElementRegion>( regionName );
      elemRegion->SetDocumentationNodes();
      elemRegion->ReadXML(childNode);
    }
  }
}


void ElementRegionManager::InitializePreSubGroups( ManagedGroup * const )
{
//    map<string,integer> constitutiveSizes;
//    ManagedGroup * domain = problemManager.GetGroup(keys::domain);
//    forElementRegions([&]( ElementRegion& elementRegion ) -> void
//    {
//      map<string,integer> sizes = elementRegion.SetConstitutiveMap(
// problemManager );
//      for( auto& entry : sizes )
//      {
//        constitutiveSizes[entry.first] += entry.second;
//      }
//    });
//
//    ManagedGroup * constitutiveManager =
// domain->GetGroup(keys::ConstitutiveManager);
//    for( auto & material : constitutiveManager->GetSubGroups() )
//    {
//      string name = material.first;
//      if( constitutiveSizes.count(name) > 0 )
//      {
//        material.second->resize( constitutiveSizes.at(name) );
//      }
//    }
}

void ElementRegionManager::InitializePostSubGroups( ManagedGroup * const problemManager )
{
  ObjectManagerBase::InitializePostSubGroups(nullptr);

  map<string,localIndex> constitutiveSizes;
  ManagedGroup * domain = problemManager->GetGroup(keys::domain);
  forElementRegions([&]( ElementRegion * elementRegion ) -> void
    {
//      map<string,localIndex> sizes;
      elementRegion->SetConstitutiveMap(problemManager, constitutiveSizes);
//      for( auto& entry : sizes )
//      {
//        constitutiveSizes[entry.first] += entry.second;
//      }
    });

  ManagedGroup * constitutiveManager = domain->GetGroup(keys::ConstitutiveManager);
  for( auto & material : constitutiveManager->GetSubGroups() )
  {
    string name = material.first;
    if( constitutiveSizes.count(name) > 0 )
    {
      material.second->resize(constitutiveSizes.at(name));
    }
  }
}

int ElementRegionManager::PackSize( array<string> const & wrapperNames,
              ElementViewAccessor<localIndex_array> const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackPrivate<false>( junk, wrapperNames, packList );
}

int ElementRegionManager::Pack( buffer_unit_type * & buffer,
          array<string> const & wrapperNames,
          ElementViewAccessor<localIndex_array> const & packList ) const
{
  return PackPrivate<true>( buffer, wrapperNames, packList );
}

template< bool DOPACK >
int
ElementRegionManager::PackPrivate( buffer_unit_type * & buffer,
                                   array<string> const & wrapperNames,
                                   ElementViewAccessor<localIndex_array> const & packList ) const
{
  int packedSize = 0;

//  packedSize += ManagedGroup::Pack( buffer, wrapperNames, {}, 0, 0);


  packedSize += bufferOps::Pack<DOPACK>( buffer, numRegions() );

  for( typename dataRepository::indexType kReg=0 ; kReg<numRegions() ; ++kReg  )
  {
    ElementRegion const * const elemRegion = GetRegion(kReg);
    packedSize += bufferOps::Pack<DOPACK>( buffer, elemRegion->getName() );

    packedSize += bufferOps::Pack<DOPACK>( buffer, elemRegion->numSubRegions() );
    for( typename dataRepository::indexType kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion(kSubReg);
      packedSize += bufferOps::Pack<DOPACK>( buffer, subRegion->getName() );

      localIndex_array const & elemList = packList[kReg][kSubReg];
      if( DOPACK )
      {
        packedSize += subRegion->Pack( buffer, wrapperNames, elemList, 0 );
      }
      else
      {
        packedSize += subRegion->PackSize( wrapperNames, elemList, 0 );
      }
    }
  }

  return packedSize;
}
//template int ElementRegionManager::PackPrivate<true>( buffer_unit_type * &,
//                                                      array<string> const &,
//                                                      ElementViewAccessor<localIndex_array> const & ) const;
//template int ElementRegionManager::PackPrivate<false>( buffer_unit_type * &,
//                                                      array<string> const &,
//                                                      ElementViewAccessor<localIndex_array> const & ) const;

int
ElementRegionManager::Unpack( buffer_unit_type const * & buffer,
                              ElementViewAccessor<localIndex_array> & packList )
{
  int unpackedSize = 0;

  localIndex numRegionsRead;
  unpackedSize += bufferOps::Unpack( buffer, numRegionsRead );

  for( localIndex kReg=0 ; kReg<numRegionsRead ; ++kReg  )
  {
    string regionName;
    unpackedSize += bufferOps::Unpack( buffer, regionName );

    ElementRegion * const elemRegion = GetRegion(regionName);

    localIndex numSubRegionsRead;
    unpackedSize += bufferOps::Unpack( buffer, numSubRegionsRead );
    for( localIndex kSubReg=0 ; kSubReg<numSubRegionsRead ; ++kSubReg  )
    {
      string subRegionName;
      unpackedSize += bufferOps::Unpack( buffer, subRegionName );

      CellBlockSubRegion * const subRegion = elemRegion->GetSubRegion(subRegionName);

      /// THIS IS WRONG??
      localIndex_array & elemList = packList[kReg][kSubReg];

      unpackedSize += subRegion->Unpack( buffer, elemList, 0 );
    }
  }

  return unpackedSize;
}


 int ElementRegionManager::PackGlobalMapsSize( ElementViewAccessor<localIndex_array> const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackGlobalMapsPrivate<false>( junk, packList);
}

 int ElementRegionManager::PackGlobalMaps( buffer_unit_type * & buffer,
                                           ElementViewAccessor<localIndex_array> const & packList ) const
{
  return PackGlobalMapsPrivate<true>( buffer, packList);
}
template< bool DOPACK >
int
ElementRegionManager::PackGlobalMapsPrivate( buffer_unit_type * & buffer,
                                             ElementViewAccessor<localIndex_array> const & packList ) const
{
  int packedSize = 0;

  packedSize += bufferOps::Pack<DOPACK>( buffer, numRegions() );

  for( typename dataRepository::indexType kReg=0 ; kReg<numRegions() ; ++kReg  )
  {
    ElementRegion const * const elemRegion = GetRegion(kReg);
    packedSize += bufferOps::Pack<DOPACK>( buffer, elemRegion->getName() );

    packedSize += bufferOps::Pack<DOPACK>( buffer, elemRegion->numSubRegions() );
    for( typename dataRepository::indexType kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion(kSubReg);
      packedSize += bufferOps::Pack<DOPACK>( buffer, subRegion->getName() );

      localIndex_array const & elemList = packList[kReg][kSubReg];
      if( DOPACK )
      {
        packedSize += subRegion->PackGlobalMaps( buffer, elemList, 0 );
      }
      else
      {
        packedSize += subRegion->PackGlobalMapsSize( elemList, 0 );
      }
    }
  }

  return packedSize;
}




int
ElementRegionManager::UnpackGlobalMaps( buffer_unit_type const * & buffer,
                                        ElementViewAccessor<localIndex_array> & packList )
{
  int unpackedSize = 0;

  localIndex numRegionsRead;
  unpackedSize += bufferOps::Unpack( buffer, numRegionsRead );

  packList.resize(numRegionsRead);
  for( localIndex kReg=0 ; kReg<numRegionsRead ; ++kReg  )
  {
    string regionName;
    unpackedSize += bufferOps::Unpack( buffer, regionName );

    ElementRegion * const elemRegion = GetRegion(regionName);

    localIndex numSubRegionsRead;
    unpackedSize += bufferOps::Unpack( buffer, numSubRegionsRead );
    packList[kReg].resize(numSubRegionsRead);
    for( localIndex kSubReg=0 ; kSubReg<numSubRegionsRead ; ++kSubReg  )
    {
      string subRegionName;
      unpackedSize += bufferOps::Unpack( buffer, subRegionName );

      CellBlockSubRegion * const subRegion = elemRegion->GetSubRegion(subRegionName);

      /// THIS IS WRONG
      localIndex_array & elemList = packList[kReg][kSubReg];

      unpackedSize += subRegion->UnpackGlobalMaps( buffer, elemList, 0 );
    }
  }

  return unpackedSize;
}




int ElementRegionManager::PackUpDownMapsSize( ElementViewAccessor<localIndex_array> const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate<false>( junk, packList);
}

int ElementRegionManager::PackUpDownMaps( buffer_unit_type * & buffer,
                                          ElementViewAccessor<localIndex_array> const & packList ) const
{
  return PackUpDownMapsPrivate<true>( buffer, packList);
}

template< bool DOPACK >
int
ElementRegionManager::PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                             ElementViewAccessor<localIndex_array> const & packList ) const
{
  int packedSize = 0;

  packedSize += bufferOps::Pack<DOPACK>( buffer, numRegions() );

  for( typename dataRepository::indexType kReg=0 ; kReg<numRegions() ; ++kReg  )
  {
    ElementRegion const * const elemRegion = GetRegion(kReg);
    packedSize += bufferOps::Pack<DOPACK>( buffer, elemRegion->getName() );

    packedSize += bufferOps::Pack<DOPACK>( buffer, elemRegion->numSubRegions() );
    for( typename dataRepository::indexType kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion(kSubReg);
      packedSize += bufferOps::Pack<DOPACK>( buffer, subRegion->getName() );

      localIndex_array const & elemList = packList[kReg][kSubReg];
      if( DOPACK )
      {
        packedSize += subRegion->PackUpDownMaps( buffer, elemList );
      }
      else
      {
        packedSize += subRegion->PackUpDownMapsSize( elemList );
      }
    }
  }

  return packedSize;
}
//template int
//ElementRegionManager::
//PackUpDownMapsPrivate<true>( buffer_unit_type * & buffer,
//                             ElementViewAccessor<localIndex_array> const & packList ) const;
//template int
//ElementRegionManager::
//PackUpDownMapsPrivate<false>( buffer_unit_type * & buffer,
//                             ElementViewAccessor<localIndex_array> const & packList ) const;


int
ElementRegionManager::UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                        ElementViewAccessor<localIndex_array> const & packList )
{
  int unpackedSize = 0;

  localIndex numRegionsRead;
  unpackedSize += bufferOps::Unpack( buffer, numRegionsRead );

//  packList.resize(numRegionsRead);
  for( localIndex kReg=0 ; kReg<numRegionsRead ; ++kReg  )
  {
    string regionName;
    unpackedSize += bufferOps::Unpack( buffer, regionName );

    ElementRegion * const elemRegion = GetRegion(regionName);

    localIndex numSubRegionsRead;
    unpackedSize += bufferOps::Unpack( buffer, numSubRegionsRead );
//    packList[kReg].resize(numSubRegionsRead);
    for( localIndex kSubReg=0 ; kSubReg<numSubRegionsRead ; ++kSubReg  )
    {
      string subRegionName;
      unpackedSize += bufferOps::Unpack( buffer, subRegionName );

      CellBlockSubRegion * const subRegion = elemRegion->GetSubRegion(subRegionName);

      /// THIS IS WRONG
      localIndex_array const & elemList = packList[kReg][kSubReg];

      unpackedSize += subRegion->UnpackUpDownMaps( buffer, elemList );
    }
  }

  return unpackedSize;
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, ElementRegionManager, string const &, ManagedGroup * const )
}
