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

/*
 * ObjectManagerBase.cpp
 *
 *  Created on: Sep 15, 2016
 *      Author: settgast1
 */

#include "ObjectManagerBase.hpp"
#include "common/TimingMacros.hpp"

namespace geosx
{
using namespace dataRepository;

ObjectManagerBase::ObjectManagerBase( std::string const & name,
                                      ManagedGroup * const parent ):
  ManagedGroup(name,parent),
  m_sets(keys::sets,this),
  m_localToGlobalMap(),
  m_globalToLocalMap()
{

  RegisterViewWrapper(viewKeyStruct::localToGlobalMapString, &m_localToGlobalMap, false );
  RegisterViewWrapper(viewKeyStruct::globalToLocalMapString, &m_globalToLocalMap, false );

  RegisterGroup( keys::sets, &m_sets, false );
  RegisterViewWrapper(viewKeyStruct::isExternalString, &m_isExternal, false );
  RegisterViewWrapper(viewKeyStruct::ghostRankString, &m_ghostRank, false )->
    setPlotLevel(PlotLevel::LEVEL_0);

  this->RegisterGroup(m_ObjectManagerBaseGroupKeys.neighborData);

  m_sets.RegisterViewWrapper<set<localIndex>>( this->m_ObjectManagerBaseViewKeys.externalSet );
}
//ObjectManagerBase::ObjectManagerBase( std::string const & name,
//                                      ManagedGroup * const parent,
//                                      cxx_utilities::DocumentationNode *
// docNode ):
//    ManagedGroup(name,parent,docNode),
//    m_localToGlobalMap( RegisterViewWrapper< globalIndex_array
// >("localToGlobal")->reference() )
//{
//
//
//  this->RegisterGroup<ManagedGroup>("Sets");
//  this->RegisterViewWrapper< array1d<integer> >("isExternal");
//}


ObjectManagerBase::~ObjectManagerBase()
{}



ObjectManagerBase::CatalogInterface::CatalogType& ObjectManagerBase::GetCatalog()
{
  static ObjectManagerBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void ObjectManagerBase::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->AllocateChildNode( m_ObjectManagerBaseViewKeys.ghostRank.Key(),
                              m_ObjectManagerBaseViewKeys.ghostRank.Key(),
                              -1,
                              "integer_array",
                              "integer_array",
                              "Array that indicates whether or not an index is a ghost. ",
                              "Array that indicates whether or not an index is a ghost. "
                              "If it is not a ghost the value will be -1. If it "
                              "is a ghost, then the value will be the owning rank.",
                              "-1",
                              "",
                              1,
                              0,
                              0 );

  docNode->AllocateChildNode( m_ObjectManagerBaseViewKeys.domainBoundaryIndicator.Key(),
                              m_ObjectManagerBaseViewKeys.domainBoundaryIndicator.Key(),
                              -1,
                              "integer_array",
                              "integer_array",
                              "List containing the element regions of the faces",
                              "List containing the element regions of the faces",
                              "0",
                              "",
                              1,
                              0,
                              2 );

//  docNode->AllocateChildNode( viewKeys.globalToLocalMap.Key(),
//                              viewKeys.globalToLocalMap.Key(),
//                              -1,
//                              "localIndex_map",
//                              "localIndex_map",
//                              "Map from globalIndex to localIndex. ",
//                              "Map from globalIndex to localIndex. ",
//                              "",
//                              "",
//                              1,
//                              0,
//                              0 );

//
//  docNode->AllocateChildNode( viewKeys.localToGlobalMap.Key(),
//                              viewKeys.localToGlobalMap.Key(),
//                              -1,
//                              "globalIndex_array",
//                              "globalIndex_array",
//                              "Array that maps from localIndex to globalIndex. ",
//                              "Array that maps from localIndex to globalIndex. ",
//                              "",
//                              "",
//                              1,
//                              0,
//                              0 );

}

void ObjectManagerBase::InitializePostSubGroups( ManagedGroup * const )
{
  m_ghostRank=-2;
}


void ObjectManagerBase::CreateSet( const std::string& newSetName )
{
  ManagedGroup * sets = GetGroup(std::string("Sets"));
  sets->RegisterViewWrapper<set<localIndex>>(newSetName);
}

void ObjectManagerBase::ConstructSetFromSetAndMap( const set<localIndex>& inputSet,
                                                   const array2d<localIndex>& map,
                                                   const std::string& setName )
{
  ManagedGroup * sets = GetGroup(std::string("Sets"));
  set<localIndex>& newset = sets->getReference<set<localIndex>>(setName);
  newset.clear();

  localIndex mapSize = map.size(1);
  for( localIndex ka=0 ; ka<size() ; ++ka )
  {
    localIndex addToSet = 0;
    for( localIndex a=0 ; a<mapSize ; ++a )
    {
      if( inputSet.count( map[ka][a] ) == 1 )
      {
        ++addToSet;
      }
    }
    if( addToSet == mapSize )
    {
      newset.insert( ka );
    }
  }
}

void ObjectManagerBase::ConstructSetFromSetAndMap( const set<localIndex>& inputSet,
                                                   const array1d<localIndex_array>& map,
                                                   const std::string& setName )
{
  ManagedGroup * sets = GetGroup(std::string("Sets"));
  set<localIndex>& newset = sets->getReference<set<localIndex>>(setName);
  newset.clear();

  for( localIndex ka=0 ; ka<size() ; ++ka )
  {
    localIndex addToSet = 0;
    localIndex mapSize = map[ka].size();
    for( localIndex a=0 ; a<mapSize ; ++a )
    {
      if( inputSet.count( map[ka][a] ) == 1 )
      {
        ++addToSet;
      }
    }
    if( addToSet == mapSize )
    {
      newset.insert( ka );
    }
  }
}

void ObjectManagerBase::ConstructLocalListOfBoundaryObjects( localIndex_array& objectList ) const
{
  const array1d<integer>& isDomainBoundary = this->getReference<integer_array>(m_ObjectManagerBaseViewKeys.domainBoundaryIndicator);
  for( localIndex k=0 ; k<size() ; ++k )
  {
    if( isDomainBoundary[k] == 1 )
    {
      objectList.push_back( k );
    }
  }
}

void ObjectManagerBase::ConstructGlobalListOfBoundaryObjects( globalIndex_array& objectList ) const
{
  const array1d<integer>& isDomainBoundary = this->getReference<integer_array>(m_ObjectManagerBaseViewKeys.domainBoundaryIndicator);
  for( localIndex k=0 ; k<size() ; ++k )
  {
    if( isDomainBoundary[k] == 1 )
    {
      objectList.push_back( this->m_localToGlobalMap[k] );
    }
  }
  std::sort( objectList.begin(), objectList.end() );
}

void ObjectManagerBase::ConstructGlobalToLocalMap()
{
  m_globalToLocalMap.clear();
  for( localIndex k=0 ; k<size() ; ++k )
  {
    m_globalToLocalMap[m_localToGlobalMap[k]] = k;
  }

}






localIndex ObjectManagerBase::PackSize( string_array const & wrapperNames,
                                        arrayView1d<localIndex const> const & packList,
                                        integer const recursive ) const
{
  localIndex packedSize = 0;
  buffer_unit_type * junk;
  packedSize += this->PackPrivate<false>( junk,
                                          wrapperNames,
                                          packList,
                                          recursive );

  return packedSize;
}




localIndex ObjectManagerBase::Pack( buffer_unit_type * & buffer,
                                    string_array const & wrapperNames,
                                    arrayView1d<localIndex const> const & packList,
                                    integer const recursive ) const
{
  localIndex packedSize = 0;

  packedSize += this->PackPrivate<true>( buffer, wrapperNames, packList, recursive );

  return packedSize;
}

template< bool DOPACK >
localIndex ObjectManagerBase::PackPrivate( buffer_unit_type * & buffer,
                                           string_array const & wrapperNames,
                                           arrayView1d<localIndex const> const & packList,
                                           integer const recursive ) const
{
  localIndex packedSize = 0;
  packedSize += bufferOps::Pack<DOPACK>( buffer, this->getName() );

  int rank=0;
  MPI_Comm_rank(MPI_COMM_GEOSX, &rank );
  packedSize += bufferOps::Pack<DOPACK>( buffer, rank );


  localIndex const numPackedIndices = packList.size()==0 ? this->size() : packList.size();
  packedSize += bufferOps::Pack<DOPACK>( buffer, numPackedIndices );

  packedSize += bufferOps::Pack<DOPACK>( buffer, string("Wrappers") );


  string_array wrapperNamesForPacking;
  if( wrapperNames.size()==0 )
  {
    set<localIndex> exclusionList;
    ViewPackingExclusionList(exclusionList);
    wrapperNamesForPacking.resize( this->wrappers().size() );
    localIndex count = 0;
    for( localIndex k=0 ; k<this->wrappers().size() ; ++k )
    {
      if( exclusionList.count(k) == 0)
      {
        wrapperNamesForPacking[count++] = wrappers().values()[k].first;
      }
    }
    wrapperNamesForPacking.resize(count);
  }
  else
  {
    wrapperNamesForPacking = wrapperNames;
  }

  packedSize += bufferOps::Pack<DOPACK>( buffer, wrapperNamesForPacking.size() );
  for( auto const & wrapperName : wrapperNamesForPacking )
  {
    dataRepository::ViewWrapperBase const * const wrapper = this->getWrapperBase(wrapperName);
    packedSize += bufferOps::Pack<DOPACK>( buffer, wrapperName );
    if( packList.empty() )
    {
      if(DOPACK)
      {
        packedSize += wrapper->Pack( buffer );
      }
      else
      {
        packedSize += wrapper->PackSize();
      }
    }
    else
    {
      if(DOPACK)
      {
        packedSize += wrapper->Pack( buffer, packList );
      }
      else
      {
        packedSize += wrapper->PackSize( packList );
      }
    }
  }


  if( recursive > 0 )
  {
    packedSize += bufferOps::Pack<DOPACK>( buffer, string("SubGroups") );
    packedSize += bufferOps::Pack<DOPACK>( buffer, this->GetSubGroups().size() );
    for( auto const & keyGroupPair : this->GetSubGroups() )
    {
      packedSize += bufferOps::Pack<DOPACK>( buffer, keyGroupPair.first );
      packedSize += keyGroupPair.second->Pack( buffer, wrapperNames, packList, recursive );
    }
  }

  return packedSize;
}
//template int ObjectManagerBase::PackPrivate<true>( buffer_unit_type * & ,string_array const & , localIndex_array const &, integer const ) const;
//template int ObjectManagerBase::PackPrivate<false>( buffer_unit_type * & ,string_array const & , localIndex_array const &, integer const ) const;



localIndex ObjectManagerBase::Unpack( buffer_unit_type const *& buffer,
                                      arrayView1d<localIndex> & packList,
                                      integer const recursive )
{
  localIndex unpackedSize = 0;
  string groupName;
  unpackedSize += bufferOps::Unpack( buffer, groupName );
  GEOS_ERROR_IF( groupName != this->getName(), "ObjectManagerBase::Unpack(): group names do not match");

  int rank=0;
  MPI_Comm_rank(MPI_COMM_GEOSX, &rank );
  int sendingRank;
  unpackedSize += bufferOps::Unpack( buffer, sendingRank );

  localIndex numUnpackedIndices;
  unpackedSize += bufferOps::Unpack( buffer, numUnpackedIndices );

//  integer readIncludeGlobalIndices;
//  unpackedSize += CommBufferOps::Unpack( buffer, readIncludeGlobalIndices);

  string wrappersLabel;
  unpackedSize += bufferOps::Unpack( buffer, wrappersLabel);
  GEOS_ERROR_IF( wrappersLabel != "Wrappers", "ObjectManagerBase::Unpack(): wrapper label incorrect");


//  if( readIncludeGlobalIndices )
//  {
//    localIndex_array unpackedLocalIndices;
//    unpackedLocalIndices.resize(numUnpackedIndices);
//
//  }

  localIndex numWrappers;
  unpackedSize += bufferOps::Unpack( buffer, numWrappers);
  for( localIndex a=0 ; a<numWrappers ; ++a )
  {
    string wrapperName;
    unpackedSize += bufferOps::Unpack( buffer, wrapperName );
    ViewWrapperBase * const wrapper = this->getWrapperBase(wrapperName);
    unpackedSize += wrapper->Unpack(buffer,packList);
  }


  if( recursive > 0 )
  {
    string subGroups;
    unpackedSize += bufferOps::Unpack( buffer, subGroups );
    GEOS_ERROR_IF( subGroups != "SubGroups", "ManagedGroup::Unpack(): group names do not match");

    decltype( this->GetSubGroups().size()) numSubGroups;
    unpackedSize += bufferOps::Unpack( buffer, numSubGroups );
    GEOS_ERROR_IF( numSubGroups != this->GetSubGroups().size(), "ManagedGroup::Unpack(): incorrect number of subGroups");

    for( auto const & index : this->GetSubGroups() )
    {
      string subGroupName;
      unpackedSize += bufferOps::Unpack( buffer, subGroupName );
      unpackedSize += this->GetGroup(subGroupName)->Unpack(buffer,packList,recursive);
    }
  }

  return unpackedSize;
}


localIndex ObjectManagerBase::PackGlobalMapsSize( arrayView1d<localIndex> const & packList,
                                integer const recursive ) const
{
  buffer_unit_type * junk = nullptr;
  return PackGlobalMapsPrivate<false>( junk, packList, recursive);
}

localIndex ObjectManagerBase::PackGlobalMaps( buffer_unit_type * & buffer,
                            arrayView1d<localIndex> const & packList,
                            integer const recursive ) const
{
  return PackGlobalMapsPrivate<true>( buffer, packList, recursive);
}

template< bool DOPACK >
localIndex ObjectManagerBase::PackGlobalMapsPrivate( buffer_unit_type * & buffer,
                                                     arrayView1d<localIndex const> const & packList,
                                                     integer const recursive ) const
{
  localIndex packedSize = bufferOps::Pack<DOPACK>( buffer, this->getName() );

  // this doesn't link without the string()...no idea why.
  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::localToGlobalMapString) );

  int rank=0;
  MPI_Comm_rank(MPI_COMM_GEOSX, &rank );
  packedSize += bufferOps::Pack<DOPACK>( buffer, rank );

  localIndex const numPackedIndices = packList.size()==0 ? this->size() : packList.size();
  packedSize += bufferOps::Pack<DOPACK>( buffer, numPackedIndices );

  globalIndex_array globalIndices;
  globalIndices.resize(numPackedIndices);
  if( packList.size()==0 )
  {
    for( localIndex a=0 ; a<numPackedIndices ; ++a )
    {
      globalIndices[a] = this->m_localToGlobalMap[a];
    }
  }
  else
  {
    for( localIndex a=0 ; a<numPackedIndices ; ++a )
    {
      globalIndices[a] = this->m_localToGlobalMap[packList[a]];
    }
  }
  packedSize += bufferOps::Pack<DOPACK>( buffer, globalIndices );


//  packedSize += CommBufferOps::Pack<DOPACK>( buffer, recursive );
  if( recursive > 0 )
  {
    packedSize += bufferOps::Pack<DOPACK>( buffer, string("SubGroups") );
    packedSize += bufferOps::Pack<DOPACK>( buffer, this->GetSubGroups().size() );
    for( auto const & keyGroupPair : this->GetSubGroups() )
    {
      packedSize += bufferOps::Pack<DOPACK>( buffer, keyGroupPair.first );
      ObjectManagerBase const * const subObjectManager = ManagedGroup::group_cast<ObjectManagerBase const *>(keyGroupPair.second);
      if( subObjectManager )
      {
        packedSize += subObjectManager->PackGlobalMapsPrivate<DOPACK>( buffer, packList, recursive );
      }
    }
  }

  return packedSize;
}
//template int ObjectManagerBase::PackGlobalMapsPrivate<true>( buffer_unit_type * &, localIndex_array const &, integer const ) const;
//template int ObjectManagerBase::PackGlobalMapsPrivate<false>( buffer_unit_type * &, localIndex_array const &, integer const ) const;

localIndex ObjectManagerBase::UnpackGlobalMaps( buffer_unit_type const *& buffer,
                                         localIndex_array & packList,
                                         integer const recursive )
{
  localIndex unpackedSize = 0;
  string groupName;
  unpackedSize += bufferOps::Unpack( buffer, groupName );
  GEOS_ERROR_IF( groupName != this->getName(), "ObjectManagerBase::Unpack(): group names do not match");

  string localToGlobalString;
  unpackedSize += bufferOps::Unpack( buffer, localToGlobalString);
  GEOS_ERROR_IF( localToGlobalString != viewKeyStruct::localToGlobalMapString, "ObjectManagerBase::Unpack(): label incorrect");

  int rank=0;
  MPI_Comm_rank(MPI_COMM_GEOSX, &rank );
  int sendingRank;
  unpackedSize += bufferOps::Unpack( buffer, sendingRank );

  localIndex numUnpackedIndices;
  unpackedSize += bufferOps::Unpack( buffer, numUnpackedIndices );

  localIndex_array unpackedLocalIndices;
  unpackedLocalIndices.resize(numUnpackedIndices);

  globalIndex_array globalIndices;
  unpackedSize += bufferOps::Unpack( buffer, globalIndices );
  localIndex numNewIndices = 0;
  globalIndex_array newGlobalIndices;
  localIndex const oldSize = this->size();
  for( localIndex a=0 ; a<numUnpackedIndices ; ++a )
  {
    // check to see if the object already exists by checking for the global
    // index in m_globalToLocalMap. If it doesn't, then add the object
    map<globalIndex,localIndex>::iterator iterG2L = m_globalToLocalMap.find(globalIndices[a]);
    if( iterG2L == m_globalToLocalMap.end() )
    {
      // object does not exist on this domain
      const localIndex newLocalIndex = oldSize + numNewIndices;

      // add the global index of the new object to the globalToLocal map
      m_globalToLocalMap[globalIndices[a]] = newLocalIndex;

      unpackedLocalIndices(a) = newLocalIndex;

      newGlobalIndices.push_back( globalIndices[a] );

      ++numNewIndices;

      GEOS_ERROR_IF( packList.size() != 0,
                   "ObjectManagerBase::Unpack(): packList specified, "
                   "but a new globalIndex is unpacked");
    }
    else
    {
      // object already exists on this domain
      // get the local index of the node
      localIndex b = iterG2L->second;
      unpackedLocalIndices(a) = b;
      if( ( sendingRank < rank && m_ghostRank[b] <= -1) || ( sendingRank < m_ghostRank[b] ) )
      {
        m_ghostRank[b] = sendingRank;
      }
    }
  }
  newGlobalIndices.resize(numNewIndices);

  // figure out new size of object container, and resize it
  const localIndex newSize = oldSize + numNewIndices;
  this->resize( newSize );

  // add the new indices to the maps.
  for( int a=0 ; a<numNewIndices ; ++a )
  {
    localIndex const b = oldSize + a;
    m_localToGlobalMap[b] = newGlobalIndices(a);
    m_ghostRank[b] = sendingRank;
  }



//    for( localIndex a=0 ; a<numUnpackedIndices ; ++a )
//    {
//      localIndex const b = m_globalToLocalMap.at(globalIndices[a]);
//      if( m_ghostRank[b]==-1 )
//      {
//        if( sendingRank < rank )
//        {
//          m_ghostRank[b] = sendingRank;
//        }
//      }
//      else
//      {
//        if( sendingRank < m_ghostRank[b] )
//        {
//          m_ghostRank[b] = sendingRank;
//        }
//
//      }
//    }

    packList = unpackedLocalIndices;


  if( recursive > 0 )
  {
    string subGroups;
    unpackedSize += bufferOps::Unpack( buffer, subGroups );
    GEOS_ERROR_IF( subGroups != "SubGroups", "ManagedGroup::Unpack(): group names do not match");

    decltype( this->GetSubGroups().size()) numSubGroups;
    unpackedSize += bufferOps::Unpack( buffer, numSubGroups );
    GEOS_ERROR_IF( numSubGroups != this->GetSubGroups().size(), "ManagedGroup::Unpack(): incorrect number of subGroups");

    for( auto const & index : this->GetSubGroups() )
    {
      string subGroupName;
      unpackedSize += bufferOps::Unpack( buffer, subGroupName );
      unpackedSize += this->GetGroup<ObjectManagerBase>(subGroupName)->UnpackGlobalMaps(buffer,packList,recursive);
    }
  }

  return unpackedSize;
}



void ObjectManagerBase::ViewPackingExclusionList( set<localIndex> & exclusionList ) const
{
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::localToGlobalMapString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::globalToLocalMapString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::ghostRankString));
}


localIndex ObjectManagerBase::GetNumberOfGhosts() const
{
  localIndex rval = 0;
  for( localIndex i=0 ; i<size() ; ++i )
  {
    if( m_ghostRank[i] > -1 )
    {
      ++rval;
    }
  }
  return rval;
//  return std::count_if( m_ghostRank.begin(), m_ghostRank.end(), [](integer i)->localIndex {return i>-1;} );
}

localIndex ObjectManagerBase::GetNumberOfLocalIndices() const
{
  localIndex rval = 0;
  for( localIndex i=0 ; i<size() ; ++i )
  {
    if( m_ghostRank[i] <= -1 )
    {
      ++rval;
    }
  }
  return rval;
  //return std::count_if( m_ghostRank.begin(), m_ghostRank.end(), [](integer i)->localIndex {return i==-1;} );
}

void ObjectManagerBase::SetReceiveLists(  )
{

  map<int,localIndex_array>  receiveIndices;
  for( localIndex a=0 ; a<size() ; ++a )
  {
    if( m_ghostRank[a] > -1 )
    {
      receiveIndices[m_ghostRank[a]].push_back(a);
    }
  }

  for( map<int,localIndex_array>::const_iterator iter=receiveIndices.begin() ; iter!=receiveIndices.end() ; ++iter )
  {
    ManagedGroup * const neighborData = GetGroup(m_ObjectManagerBaseGroupKeys.neighborData)->GetGroup( std::to_string( iter->first ) );

    localIndex_array & nodeAdjacencyList = neighborData->getReference<localIndex_array>( m_ObjectManagerBaseViewKeys.ghostsToReceive );
    nodeAdjacencyList = iter->second;
  }

}

integer ObjectManagerBase::SplitObject( localIndex const indexToSplit,
                                        int const rank,
                                        localIndex & newIndex )
{

  // if the object index has a zero sized childIndices entry, then this object can be split into two
  // new objects

  if( size()+1 > capacity() )
  {
    reserve( static_cast<localIndex>( size() * m_overAllocationFactor ) );
  }

  // the new indices are tacked on to the end of the arrays
  newIndex = size() ;
  this->resize( newIndex + 1 );

  // copy the fields
  CopyObject( indexToSplit, newIndex );

  localIndex_array * const
  parentIndex = this->getPointer<localIndex_array>( m_ObjectManagerBaseViewKeys.parentIndex );
  if( parentIndex != nullptr )
  {
    (*parentIndex)[newIndex] = indexToSplit;
  }

  localIndex_array * const
  childIndex = this->getPointer<localIndex_array>( m_ObjectManagerBaseViewKeys.childIndex );
  if( childIndex != nullptr )
  {
    (*childIndex)[indexToSplit] = newIndex;
  }

  m_localToGlobalMap[newIndex] = -1;

  if( m_isExternal[indexToSplit]==1 )
  {
    m_isExternal[indexToSplit] = 1;
    m_isExternal[newIndex]     = 1;
  }
  else
  {
    m_isExternal[indexToSplit] = 2;
    m_isExternal[newIndex]     = 2;
  }

  return 1;

}

void ObjectManagerBase::CopyObject( const localIndex source, const localIndex destination )
{
  for( auto & wrapper : wrappers() )
  {
    wrapper.second->copy( source, destination );
  }

  for( localIndex i=0 ; i<m_sets.wrappers().size() ; ++i )
  {
    set<localIndex> & targetSet = m_sets.getReference< set<localIndex> >(i);
    if( targetSet.count(source) > 0 )
    {
      targetSet.insert(destination);
    }
  }
}



} /* namespace geosx */
