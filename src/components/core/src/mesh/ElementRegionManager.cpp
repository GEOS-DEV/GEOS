//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
// SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S.
// Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National
// Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy,
// completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe
// privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or
// services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or
// favoring by the United States Government or Lawrence Livermore National
// Security,
//     LLC. The views and opinions of authors expressed herein do not
// necessarily state or reflect those of the United States Government or
// Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising
// or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The
// BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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


template< bool DOPACK >
int
ElementRegionManager::PackPrivate( buffer_unit_type * & buffer,
                                   array<string> const & wrapperNames,
                                   ElementViewAccessor<localIndex_array> const & packList ) const
{
  int packedSize = 0;

//  packedSize += ManagedGroup::Pack( buffer, wrapperNames, {}, 0, 0);


  packedSize += CommBufferOps::Pack<DOPACK>( buffer, numRegions() );

  for( typename dataRepository::indexType kReg=0 ; kReg<numRegions() ; ++kReg  )
  {
    ElementRegion const * const elemRegion = GetRegion(kReg);
    packedSize += CommBufferOps::Pack<DOPACK>( buffer, elemRegion->getName() );

    packedSize += CommBufferOps::Pack<DOPACK>( buffer, elemRegion->numSubRegions() );
    for( typename dataRepository::indexType kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion(kSubReg);
      packedSize += CommBufferOps::Pack<DOPACK>( buffer, subRegion->getName() );

      localIndex_array const & elemList = *(packList[kReg][kSubReg]);
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
template int ElementRegionManager::PackPrivate<true>( buffer_unit_type * &,
                                                      array<string> const &,
                                                      ElementViewAccessor<localIndex_array> const & ) const;
template int ElementRegionManager::PackPrivate<false>( buffer_unit_type * &,
                                                      array<string> const &,
                                                      ElementViewAccessor<localIndex_array> const & ) const;

int
ElementRegionManager::Unpack( buffer_unit_type const * & buffer,
                              ElementViewAccessor<localIndex_array> & packList )
{
  int unpackedSize = 0;

  localIndex numRegionsRead;
  unpackedSize += CommBufferOps::Unpack( buffer, numRegionsRead );

  for( localIndex kReg=0 ; kReg<numRegionsRead ; ++kReg  )
  {
    string regionName;
    unpackedSize += CommBufferOps::Unpack( buffer, regionName );

    ElementRegion * const elemRegion = GetRegion(regionName);

    localIndex numSubRegionsRead;
    unpackedSize += CommBufferOps::Unpack( buffer, numSubRegionsRead );
    for( localIndex kSubReg=0 ; kSubReg<numSubRegionsRead ; ++kSubReg  )
    {
      string subRegionName;
      unpackedSize += CommBufferOps::Unpack( buffer, subRegionName );

      CellBlockSubRegion * const subRegion = elemRegion->GetSubRegion(subRegionName);

      /// THIS IS WRONG??
      localIndex_array & elemList = *(packList[kReg][kSubReg]);

      unpackedSize += subRegion->Unpack( buffer, elemList, 0 );
    }
  }

  return unpackedSize;
}


template< bool DOPACK >
int
ElementRegionManager::PackGlobalMapsPrivate( buffer_unit_type * & buffer,
                                             ElementViewAccessor<localIndex_array> const & packList ) const
{
  int packedSize = 0;

  packedSize += CommBufferOps::Pack<DOPACK>( buffer, numRegions() );

  for( typename dataRepository::indexType kReg=0 ; kReg<numRegions() ; ++kReg  )
  {
    ElementRegion const * const elemRegion = GetRegion(kReg);
    packedSize += CommBufferOps::Pack<DOPACK>( buffer, elemRegion->getName() );

    packedSize += CommBufferOps::Pack<DOPACK>( buffer, elemRegion->numSubRegions() );
    for( typename dataRepository::indexType kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion(kSubReg);
      packedSize += CommBufferOps::Pack<DOPACK>( buffer, subRegion->getName() );

      localIndex_array const & elemList = *(packList[kReg][kSubReg]);
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
template int
ElementRegionManager::PackGlobalMapsPrivate<true>( buffer_unit_type * & buffer,
                                                   ElementViewAccessor<localIndex_array> const & packList ) const;
template int
ElementRegionManager::PackGlobalMapsPrivate<false>( buffer_unit_type * & buffer,
                                                   ElementViewAccessor<localIndex_array> const & packList ) const;





int
ElementRegionManager::UnpackGlobalMaps( buffer_unit_type const * & buffer,
                                        ElementViewAccessor<localIndex_array> & packList )
{
  int unpackedSize = 0;

  localIndex numRegionsRead;
  unpackedSize += CommBufferOps::Unpack( buffer, numRegionsRead );

  packList.resize(numRegionsRead);
  for( localIndex kReg=0 ; kReg<numRegionsRead ; ++kReg  )
  {
    string regionName;
    unpackedSize += CommBufferOps::Unpack( buffer, regionName );

    ElementRegion * const elemRegion = GetRegion(regionName);

    localIndex numSubRegionsRead;
    unpackedSize += CommBufferOps::Unpack( buffer, numSubRegionsRead );
    packList[kReg].resize(numSubRegionsRead);
    for( localIndex kSubReg=0 ; kSubReg<numSubRegionsRead ; ++kSubReg  )
    {
      string subRegionName;
      unpackedSize += CommBufferOps::Unpack( buffer, subRegionName );

      CellBlockSubRegion * const subRegion = elemRegion->GetSubRegion(subRegionName);

      /// THIS IS WRONG
      localIndex_array & elemList = *(packList[kReg][kSubReg]);

      unpackedSize += subRegion->UnpackGlobalMaps( buffer, elemList, 0 );
    }
  }

  return unpackedSize;
}





template< bool DOPACK >
int
ElementRegionManager::PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                             ElementViewAccessor<localIndex_array> const & packList ) const
{
  int packedSize = 0;

  packedSize += CommBufferOps::Pack<DOPACK>( buffer, numRegions() );

  for( typename dataRepository::indexType kReg=0 ; kReg<numRegions() ; ++kReg  )
  {
    ElementRegion const * const elemRegion = GetRegion(kReg);
    packedSize += CommBufferOps::Pack<DOPACK>( buffer, elemRegion->getName() );

    packedSize += CommBufferOps::Pack<DOPACK>( buffer, elemRegion->numSubRegions() );
    for( typename dataRepository::indexType kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion(kSubReg);
      packedSize += CommBufferOps::Pack<DOPACK>( buffer, subRegion->getName() );

      localIndex_array const & elemList = *(packList[kReg][kSubReg]);
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
template int
ElementRegionManager::
PackUpDownMapsPrivate<true>( buffer_unit_type * & buffer,
                             ElementViewAccessor<localIndex_array> const & packList ) const;
template int
ElementRegionManager::
PackUpDownMapsPrivate<false>( buffer_unit_type * & buffer,
                             ElementViewAccessor<localIndex_array> const & packList ) const;


int
ElementRegionManager::UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                        ElementViewAccessor<localIndex_array> const & packList )
{
  int unpackedSize = 0;

  localIndex numRegionsRead;
  unpackedSize += CommBufferOps::Unpack( buffer, numRegionsRead );

//  packList.resize(numRegionsRead);
  for( localIndex kReg=0 ; kReg<numRegionsRead ; ++kReg  )
  {
    string regionName;
    unpackedSize += CommBufferOps::Unpack( buffer, regionName );

    ElementRegion * const elemRegion = GetRegion(regionName);

    localIndex numSubRegionsRead;
    unpackedSize += CommBufferOps::Unpack( buffer, numSubRegionsRead );
//    packList[kReg].resize(numSubRegionsRead);
    for( localIndex kSubReg=0 ; kSubReg<numSubRegionsRead ; ++kSubReg  )
    {
      string subRegionName;
      unpackedSize += CommBufferOps::Unpack( buffer, subRegionName );

      CellBlockSubRegion * const subRegion = elemRegion->GetSubRegion(subRegionName);

      /// THIS IS WRONG
      localIndex_array & elemList = *(packList[kReg][kSubReg]);

      unpackedSize += subRegion->UnpackUpDownMaps( buffer, elemList );
    }
  }

  return unpackedSize;
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, ElementRegionManager, string const &, ManagedGroup * const )
}
