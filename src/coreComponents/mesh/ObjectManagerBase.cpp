/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ObjectManagerBase.cpp
 */

#include "ObjectManagerBase.hpp"

#include "common/TimingMacros.hpp"
#include "mesh/ExtrinsicMeshData.hpp"
#include "common/MpiWrapper.hpp"

namespace geosx
{
using namespace dataRepository;

ObjectManagerBase::ObjectManagerBase( string const & name,
                                      Group * const parent ):
  Group( name, parent ),
  m_sets( groupKeyStruct::setsString(), this ),
  m_neighborGroup( groupKeyStruct::neighborDataString(), this ),
  m_localToGlobalMap(),
  m_globalToLocalMap(),
  m_isExternal(),
  m_domainBoundaryIndicator(),
  m_ghostRank(),
  m_neighborData()
{
  registerGroup( groupKeyStruct::setsString(), &m_sets );
  registerGroup( groupKeyStruct::neighborDataString(), &m_neighborGroup );

  registerWrapper( viewKeyStruct::localToGlobalMapString(), &m_localToGlobalMap ).
    setApplyDefaultValue( -1 ).
    setDescription( "Array that contains a map from localIndex to globalIndex." );

  registerWrapper( viewKeyStruct::globalToLocalMapString(), &m_globalToLocalMap );

  registerWrapper( viewKeyStruct::isExternalString(), &m_isExternal );

  registerWrapper( viewKeyStruct::ghostRankString(), &m_ghostRank ).
    setApplyDefaultValue( -2 ).
    setPlotLevel( PlotLevel::LEVEL_0 );

  registerWrapper< array1d< integer > >( viewKeyStruct::domainBoundaryIndicatorString(), &m_domainBoundaryIndicator );

  m_sets.registerWrapper< SortedArray< localIndex > >( this->m_ObjectManagerBaseViewKeys.externalSet );

  excludeWrappersFromPacking( { viewKeyStruct::localToGlobalMapString(),
                                viewKeyStruct::globalToLocalMapString(),
                                viewKeyStruct::ghostRankString(),
                                extrinsicMeshData::ParentIndex::key(),
                                extrinsicMeshData::ChildIndex::key() } );
}

ObjectManagerBase::~ObjectManagerBase()
{}



ObjectManagerBase::CatalogInterface::CatalogType & ObjectManagerBase::getCatalog()
{
  static ObjectManagerBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

SortedArray< localIndex > & ObjectManagerBase::createSet( const string & newSetName )
{
  return m_sets.registerWrapper< SortedArray< localIndex > >( newSetName ).reference();
}

void ObjectManagerBase::constructSetFromSetAndMap( SortedArrayView< localIndex const > const & inputSet,
                                                   const array2d< localIndex > & map,
                                                   const string & setName )
{
  SortedArray< localIndex > & newset = m_sets.getReference< SortedArray< localIndex > >( setName );
  newset.clear();

  localIndex const numObjects = size();
  GEOSX_ERROR_IF( map.size( 0 ) != numObjects, "Size mismatch. " << map.size( 0 ) << " != " << numObjects );

  if( setName == "all" )
  {
    newset.reserve( numObjects );

    for( localIndex ka=0; ka<numObjects; ++ka )
    {
      newset.insert( ka );
    }
  }
  else
  {
    localIndex const mapSize = map.size( 1 );
    for( localIndex ka=0; ka<numObjects; ++ka )
    {
      if( std::all_of( &map( ka, 0 ), &map( ka, 0 ) + mapSize, [&]( localIndex const i ) { return inputSet.contains( i ); } ) )
      {
        newset.insert( ka );
      }
    }
  }
}

void ObjectManagerBase::constructSetFromSetAndMap( SortedArrayView< localIndex const > const & inputSet,
                                                   const array1d< localIndex_array > & map,
                                                   const string & setName )
{
  SortedArray< localIndex > & newset = m_sets.getReference< SortedArray< localIndex > >( setName );
  newset.clear();

  localIndex const numObjects = size();
  GEOSX_ERROR_IF( map.size() != numObjects, "Size mismatch. " << map.size() << " != " << numObjects );

  if( setName == "all" )
  {
    newset.reserve( numObjects );

    for( localIndex ka=0; ka<numObjects; ++ka )
    {
      newset.insert( ka );
    }
  }
  else
  {
    for( localIndex ka=0; ka<numObjects; ++ka )
    {
      if( std::all_of( map[ka].begin(), map[ka].end(), [&]( localIndex const i ) { return inputSet.contains( i ); } ) )
      {
        newset.insert( ka );
      }
    }
  }
}

void ObjectManagerBase::constructSetFromSetAndMap( SortedArrayView< localIndex const > const & inputSet,
                                                   ArrayOfArraysView< localIndex const > const & map,
                                                   const string & setName )
{
  SortedArray< localIndex > & newset = m_sets.getReference< SortedArray< localIndex > >( setName );
  newset.clear();

  localIndex const numObjects = size();
  GEOSX_ERROR_IF( map.size() != numObjects, "Size mismatch. " << map.size() << " != " << numObjects );

  if( setName == "all" )
  {
    newset.reserve( numObjects );

    for( localIndex ka=0; ka<numObjects; ++ka )
    {
      newset.insert( ka );
    }
  }
  else
  {
    for( localIndex ka=0; ka<numObjects; ++ka )
    {
      localIndex const * const values = map[ka];
      localIndex const numValues = map.sizeOfArray( ka );
      if( std::all_of( values, values + numValues, [&]( localIndex const i ) { return inputSet.contains( i ); } ) )
      {
        newset.insert( ka );
      }
    }
  }
}

void ObjectManagerBase::constructLocalListOfBoundaryObjects( localIndex_array & objectList ) const
{
  arrayView1d< integer const > const & isDomainBoundary = this->getDomainBoundaryIndicator();
  for( localIndex k=0; k<size(); ++k )
  {
    if( isDomainBoundary[k] == 1 )
    {
      objectList.emplace_back( k );
    }
  }
}

void ObjectManagerBase::constructGlobalListOfBoundaryObjects( globalIndex_array & objectList ) const
{
  arrayView1d< integer const > const & isDomainBoundary = this->getDomainBoundaryIndicator();
  for( localIndex k=0; k<size(); ++k )
  {
    if( isDomainBoundary[k] == 1 )
    {
      objectList.emplace_back( this->m_localToGlobalMap[k] );
    }
  }
  std::sort( objectList.begin(), objectList.end() );
}

void ObjectManagerBase::constructGlobalToLocalMap()
{
  GEOSX_MARK_FUNCTION;

  m_globalToLocalMap.clear();
  localIndex const N = size();
  m_globalToLocalMap.reserve( N );
  for( localIndex k = 0; k < N; ++k )
  {
    updateGlobalToLocalMap( k );
  }
}

localIndex ObjectManagerBase::packSize( string_array const & wrapperNames,
                                        arrayView1d< localIndex const > const & packList,
                                        integer const recursive,
                                        bool onDevice,
                                        parallelDeviceEvents & events ) const
{
  buffer_unit_type * junk;
  return this->packImpl< false >( junk,
                                  wrapperNames,
                                  packList,
                                  recursive,
                                  onDevice,
                                  events );
}

localIndex ObjectManagerBase::pack( buffer_unit_type * & buffer,
                                    string_array const & wrapperNames,
                                    arrayView1d< localIndex const > const & packList,
                                    integer const recursive,
                                    bool onDevice,
                                    parallelDeviceEvents & events ) const
{
  return this->packImpl< true >( buffer, wrapperNames, packList, recursive, onDevice, events );
}

template< bool DO_PACKING >
localIndex ObjectManagerBase::packImpl( buffer_unit_type * & buffer,
                                        string_array const & wrapperNames,
                                        arrayView1d< localIndex const > const & packList,
                                        integer const recursive,
                                        bool onDevice,
                                        parallelDeviceEvents & events ) const
{
  localIndex packedSize = 0;
  packedSize += bufferOps::Pack< DO_PACKING >( buffer, this->getName() );

  int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer, rank );

  localIndex const numPackedIndices = packList.size();
  packedSize += bufferOps::Pack< DO_PACKING >( buffer, numPackedIndices );
  if( numPackedIndices > 0 )
  {
    std::set< string > input;
    std::copy( wrapperNames.begin(), wrapperNames.end(), std::inserter( input, input.end() ) );

    std::set< string > const exclusion = getPackingExclusionList();
    std::set< string > const available = mapKeys< std::set >( wrappers() );

    // Only taking the wrappers that are requested and available.
    std::set< string > tmp0;
    std::set_intersection( input.cbegin(), input.cend(), available.cbegin(), available.cend(), std::inserter( tmp0, tmp0.end() ) );

    // Checking that all the requested wrappers are available.
    std::set< string > tmp1;
    std::set_difference( input.cbegin(), input.cend(), available.cbegin(), available.cend(), std::inserter( tmp1, tmp1.end() ) );
    if( !tmp1.empty() )
    {
      GEOSX_LOG( "Wrappers \"" << stringutilities::join( tmp1, ", " ) << "\" were requested from \"" << getName() << "\" but are not available." );
    }

    // Not taking the wrappers that are excluded.
    std::set< string > tmp2;
    std::set_difference( tmp0.cbegin(), tmp0.cend(), exclusion.cbegin(), exclusion.cend(), std::inserter( tmp2, tmp2.end() ) );

    // Now we build the final list.
    // No packing by index is allowed if the registered wrapper does not share the size of the owning group.
    // Hence, the sufficient (but not necessary...) condition on `wrapper.sizedFromParent()`.
    std::vector< string > tmp3;
    auto predicate = [this]( string const & wrapperName ) -> bool
    {
      return bool( this->getWrapperBase( wrapperName ).sizedFromParent() );
    };
    std::copy_if( tmp2.cbegin(), tmp2.cend(), std::back_inserter( tmp3 ), predicate );

    // Extracting the wrappers
    std::vector< WrapperBase const * > wrappers;
    auto transformer = [this]( string const & wrapperName ) -> WrapperBase const *
    {
      return &this->getWrapperBase( wrapperName );
    };
    std::transform( tmp3.cbegin(), tmp3.cend(), std::back_inserter( wrappers ), transformer );

    // Additional refactoring should be done by using `Group::packImpl` that duplicates the following pack code.
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( "Wrappers" ) );
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, LvArray::integerConversion< localIndex >( wrappers.size() ) );
    for( WrapperBase const * wrapper: wrappers )
    {
      packedSize += bufferOps::Pack< DO_PACKING >( buffer, wrapper->getName() );
      packedSize += wrapper->packByIndex< DO_PACKING >( buffer, packList, true, onDevice, events );
    }
  }

  if( recursive > 0 )
  {
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( "SubGroups" ) );
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, this->getSubGroups().size() );
    for( auto const & keyGroupPair : this->getSubGroups() )
    {
      packedSize += bufferOps::Pack< DO_PACKING >( buffer, keyGroupPair.first );
      packedSize += keyGroupPair.second->pack( buffer, wrapperNames, packList, recursive, onDevice, events );
    }
  }

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, this->getName() );

  return packedSize;
}


localIndex ObjectManagerBase::unpack( buffer_unit_type const * & buffer,
                                      arrayView1d< localIndex > & packList,
                                      integer const recursive,
                                      bool onDevice,
                                      parallelDeviceEvents & events )
{
  localIndex unpackedSize = 0;
  string groupName;
  unpackedSize += bufferOps::Unpack( buffer, groupName );
  GEOSX_ERROR_IF_NE( groupName, this->getName() );

  int sendingRank;
  unpackedSize += bufferOps::Unpack( buffer, sendingRank );

  localIndex numUnpackedIndices;
  unpackedSize += bufferOps::Unpack( buffer, numUnpackedIndices );
  if( numUnpackedIndices > 0 )
  {

    string wrappersLabel;
    unpackedSize += bufferOps::Unpack( buffer, wrappersLabel );
    GEOSX_ERROR_IF_NE( wrappersLabel, "Wrappers" );

    localIndex numWrappers;
    unpackedSize += bufferOps::Unpack( buffer, numWrappers );
    for( localIndex a=0; a<numWrappers; ++a )
    {
      string wrapperName;
      unpackedSize += bufferOps::Unpack( buffer, wrapperName );
      unpackedSize += this->getWrapperBase( wrapperName ).unpackByIndex( buffer, packList, true, onDevice, events );
    }
  }

  if( recursive > 0 )
  {
    string subGroups;
    unpackedSize += bufferOps::Unpack( buffer, subGroups );
    GEOSX_ERROR_IF_NE( subGroups, "SubGroups" );

    decltype( this->getSubGroups().size()) numSubGroups;
    unpackedSize += bufferOps::Unpack( buffer, numSubGroups );
    GEOSX_ERROR_IF_NE( numSubGroups, this->getSubGroups().size() );

    for( localIndex i = 0; i < this->numSubGroups(); ++i )
    {
      string subGroupName;
      unpackedSize += bufferOps::Unpack( buffer, subGroupName );
      unpackedSize += this->getGroup( subGroupName ).unpack( buffer, packList, recursive, onDevice, events );
    }
  }

  unpackedSize += bufferOps::Unpack( buffer, groupName );
  GEOSX_ERROR_IF_NE( groupName, this->getName() );

  return unpackedSize;
}

template< bool DO_PACKING >
localIndex ObjectManagerBase::packParentChildMapsImpl( buffer_unit_type * & buffer,
                                                       arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = 0;

  if( this->hasExtrinsicData< extrinsicMeshData::ParentIndex >() )
  {
    arrayView1d< localIndex const > const parentIndex = this->getExtrinsicData< extrinsicMeshData::ParentIndex >();
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( extrinsicMeshData::ParentIndex::key() ) );
    packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                                 parentIndex,
                                                 packList,
                                                 this->m_localToGlobalMap,
                                                 this->m_localToGlobalMap );
  }

  if( this->hasExtrinsicData< extrinsicMeshData::ChildIndex >() )
  {
    arrayView1d< localIndex const > const & childIndex = this->getExtrinsicData< extrinsicMeshData::ChildIndex >();
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( extrinsicMeshData::ChildIndex::key() ) );
    packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                                 childIndex,
                                                 packList,
                                                 this->m_localToGlobalMap,
                                                 this->m_localToGlobalMap );
  }

  return packedSize;
}

template
localIndex ObjectManagerBase::packParentChildMapsImpl< true >( buffer_unit_type * & buffer,
                                                               arrayView1d< localIndex const > const & packList ) const;
template
localIndex ObjectManagerBase::packParentChildMapsImpl< false >( buffer_unit_type * & buffer,
                                                                arrayView1d< localIndex const > const & packList ) const;


localIndex ObjectManagerBase::unpackParentChildMaps( buffer_unit_type const * & buffer,
                                                     localIndex_array & packList )
{
  localIndex unpackedSize = 0;

  if( this->hasExtrinsicData< extrinsicMeshData::ParentIndex >() )
  {
    arrayView1d< localIndex > const & parentIndex = this->getExtrinsicData< extrinsicMeshData::ParentIndex >();
    string shouldBeParentIndexString;
    unpackedSize += bufferOps::Unpack( buffer, shouldBeParentIndexString );
    GEOSX_ERROR_IF( shouldBeParentIndexString != extrinsicMeshData::ParentIndex::key(),
                    "value read from buffer is:" << shouldBeParentIndexString << ". It should be " << extrinsicMeshData::ParentIndex::key() );
    unpackedSize += bufferOps::Unpack( buffer,
                                       parentIndex,
                                       packList,
                                       this->m_globalToLocalMap,
                                       this->m_globalToLocalMap );
  }

  if( this->hasExtrinsicData< extrinsicMeshData::ChildIndex >() )
  {
    arrayView1d< localIndex > const & childIndex = this->getExtrinsicData< extrinsicMeshData::ChildIndex >();
    string shouldBeChildIndexString;
    unpackedSize += bufferOps::Unpack( buffer, shouldBeChildIndexString );
    GEOSX_ERROR_IF( shouldBeChildIndexString != extrinsicMeshData::ChildIndex::key(),
                    "value read from buffer is:" << shouldBeChildIndexString << ". It should be " << extrinsicMeshData::ChildIndex::key() );
    unpackedSize += bufferOps::Unpack( buffer,
                                       childIndex,
                                       packList,
                                       this->m_globalToLocalMap,
                                       this->m_globalToLocalMap );
  }

  return unpackedSize;
}



template< bool DO_PACKING >
localIndex ObjectManagerBase::packSets( buffer_unit_type * & buffer,
                                        arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = 0;
  packedSize += bufferOps::Pack< DO_PACKING >( buffer, m_sets.getName() );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, m_sets.wrappers().size() );
  for( auto const & wrapperIter : m_sets.wrappers() )
  {
    string const & setName = wrapperIter.first;
    SortedArrayView< localIndex const > const & currentSet = m_sets.getReference< SortedArray< localIndex > >( setName );
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, setName );
    SortedArray< globalIndex > emptySet;
    packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                                 currentSet,
                                                 packList,
                                                 emptySet.toViewConst(),
                                                 m_localToGlobalMap );
  }
  return packedSize;
}
//template<>
//localIndex ObjectManagerBase::packSets<true>( buffer_unit_type * &,
//                                              arrayView1d<localIndex const> const & );
//template<>
//localIndex ObjectManagerBase::packSets<false>( buffer_unit_type * &,
//                                               arrayView1d<localIndex const> const & );


localIndex ObjectManagerBase::unpackSets( buffer_unit_type const * & buffer )
{
  localIndex unpackedSize = 0;
  string name;
  unpackedSize += bufferOps::Unpack( buffer, name );
  GEOSX_ERROR_IF( name != m_sets.getName(), "ObjectManagerBase::UnpackSets(): group names do not match" );

  localIndex numUnpackedSets;
  unpackedSize += bufferOps::Unpack( buffer, numUnpackedSets );
  for( localIndex a=0; a<numUnpackedSets; ++a )
  {
    string setName;
    unpackedSize += bufferOps::Unpack( buffer, setName );
    SortedArray< localIndex > & targetSet = m_sets.getReference< SortedArray< localIndex > >( setName );

    SortedArray< globalIndex > junk;
    unpackedSize += bufferOps::Unpack( buffer,
                                       targetSet,
                                       junk,
                                       this->m_globalToLocalMap,
                                       false );
  }


  return unpackedSize;
}


localIndex ObjectManagerBase::packGlobalMapsSize( arrayView1d< localIndex const > const & packList,
                                                  integer const recursive ) const
{
  buffer_unit_type * junk = nullptr;
  return packGlobalMapsImpl< false >( junk, packList, recursive );
}

localIndex ObjectManagerBase::packGlobalMaps( buffer_unit_type * & buffer,
                                              arrayView1d< localIndex const > const & packList,
                                              integer const recursive ) const
{
  return packGlobalMapsImpl< true >( buffer, packList, recursive );
}

template< bool DO_PACKING >
localIndex ObjectManagerBase::packGlobalMapsImpl( buffer_unit_type * & buffer,
                                                  arrayView1d< localIndex const > const & packList,
                                                  integer const recursive ) const
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );

  localIndex packedSize = bufferOps::Pack< DO_PACKING >( buffer, this->getName() );

  // this doesn't link without the string()...no idea why.
  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::localToGlobalMapString() ) );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, rank );

  localIndex const numPackedIndices = packList.size();
  packedSize += bufferOps::Pack< DO_PACKING >( buffer, numPackedIndices );

  if( numPackedIndices > 0 )
  {
    globalIndex_array globalIndices;
    globalIndices.resize( numPackedIndices );
    for( localIndex a=0; a<numPackedIndices; ++a )
    {
      globalIndices[a] = this->m_localToGlobalMap[packList[a]];
    }
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, globalIndices );
  }

  // FIXME is this the responsibility of this instance to do this?
  if( this->hasExtrinsicData< extrinsicMeshData::ParentIndex >() )
  {
    arrayView1d< localIndex const > const & parentIndex = this->getExtrinsicData< extrinsicMeshData::ParentIndex >();
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( extrinsicMeshData::ParentIndex::key() ) );
    packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                                 parentIndex,
                                                 packList,
                                                 this->m_localToGlobalMap,
                                                 this->m_localToGlobalMap );
  }



  if( recursive > 0 )
  {
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( "SubGroups" ) );
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, this->getSubGroups().size() );
    for( auto const & keyGroupPair : this->getSubGroups() )
    {
      packedSize += bufferOps::Pack< DO_PACKING >( buffer, keyGroupPair.first );
      ObjectManagerBase const * const subObjectManager = dynamicCast< ObjectManagerBase const * >( keyGroupPair.second );
      if( subObjectManager )
      {
        packedSize += subObjectManager->packGlobalMapsImpl< DO_PACKING >( buffer, packList, recursive );
      }
    }
  }

  packedSize += packSets< DO_PACKING >( buffer, packList );

  return packedSize;
}

localIndex ObjectManagerBase::unpackGlobalMaps( buffer_unit_type const * & buffer,
                                                localIndex_array & packList,
                                                integer const recursive )
{
  GEOSX_MARK_FUNCTION;
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );

  localIndex unpackedSize = 0;
  string groupName;
  unpackedSize += bufferOps::Unpack( buffer, groupName );
  string msg = "ObjectManagerBase::unpack(): group names do not match as they are groupName = " + groupName + " and this->getName= " + this->getName();
  GEOSX_ERROR_IF( groupName != this->getName(), msg );

  string localToGlobalString;
  unpackedSize += bufferOps::Unpack( buffer, localToGlobalString );
  GEOSX_ERROR_IF( localToGlobalString != viewKeyStruct::localToGlobalMapString(), "ObjectManagerBase::unpack(): label incorrect" );

  int sendingRank;
  unpackedSize += bufferOps::Unpack( buffer, sendingRank );

  localIndex numUnpackedIndices;
  unpackedSize += bufferOps::Unpack( buffer, numUnpackedIndices );

  if( numUnpackedIndices > 0 )
  {
    localIndex_array unpackedLocalIndices;
    unpackedLocalIndices.resize( numUnpackedIndices );

    globalIndex_array globalIndices;
    unpackedSize += bufferOps::Unpack( buffer, globalIndices );
    localIndex numNewIndices = 0;
    globalIndex_array newGlobalIndices;
    newGlobalIndices.reserve( numUnpackedIndices );
    localIndex const oldSize = this->size();
    for( localIndex a = 0; a < numUnpackedIndices; ++a )
    {
      // check to see if the object already exists by checking for the global
      // index in m_globalToLocalMap. If it doesn't, then add the object
      unordered_map< globalIndex, localIndex >::iterator iterG2L =
        m_globalToLocalMap.find( globalIndices[a] );
      if( iterG2L == m_globalToLocalMap.end() )
      {
        // object does not exist on this domain
        const localIndex newLocalIndex = oldSize + numNewIndices;

        // add the global index of the new object to the globalToLocal map
        m_globalToLocalMap[globalIndices[a]] = newLocalIndex;

        unpackedLocalIndices( a ) = newLocalIndex;

        newGlobalIndices.emplace_back( globalIndices[a] );

        ++numNewIndices;

        GEOSX_ERROR_IF( packList.size() != 0,
                        "ObjectManagerBase::unpack(): packList specified, "
                        "but a new globalIndex is unpacked" );
      }
      else
      {
        // object already exists on this domain
        // get the local index of the node
        localIndex b = iterG2L->second;
        unpackedLocalIndices( a ) = b;
        if( ( sendingRank < rank && m_ghostRank[b] <= -1) || ( sendingRank < m_ghostRank[b] ) )
        {
          m_ghostRank[b] = sendingRank;
        }
      }
    }

    // figure out new size of object container, and resize it
    const localIndex newSize = oldSize + numNewIndices;
    this->resize( newSize );

    // add the new indices to the maps.
    for( int a=0; a<numNewIndices; ++a )
    {
      localIndex const b = oldSize + a;
      m_localToGlobalMap[b] = newGlobalIndices( a );
      m_ghostRank[b] = sendingRank;
    }


    packList = unpackedLocalIndices;
  }


  if( this->hasExtrinsicData< extrinsicMeshData::ParentIndex >() )
  {
    arrayView1d< localIndex > const & parentIndex = this->getExtrinsicData< extrinsicMeshData::ParentIndex >();
    string parentIndicesString;
    unpackedSize += bufferOps::Unpack( buffer, parentIndicesString );
    GEOSX_ERROR_IF( parentIndicesString != extrinsicMeshData::ParentIndex::key(), "ObjectManagerBase::unpack(): label incorrect" );
    unpackedSize += bufferOps::Unpack( buffer,
                                       parentIndex,
                                       packList,
                                       this->m_globalToLocalMap,
                                       this->m_globalToLocalMap );
  }


  if( recursive > 0 )
  {
    string subGroups;
    unpackedSize += bufferOps::Unpack( buffer, subGroups );
    GEOSX_ERROR_IF( subGroups != "SubGroups", "Group::unpack(): group names do not match" );

    decltype( this->getSubGroups().size()) numSubGroups;
    unpackedSize += bufferOps::Unpack( buffer, numSubGroups );
    GEOSX_ERROR_IF( numSubGroups != this->getSubGroups().size(), "Group::unpack(): incorrect number of subGroups" );

    for( auto const & index : this->getSubGroups() )
    {
      GEOSX_UNUSED_VAR( index );
      string subGroupName;
      unpackedSize += bufferOps::Unpack( buffer, subGroupName );
      unpackedSize += this->getGroup< ObjectManagerBase >( subGroupName ).unpackGlobalMaps( buffer, packList, recursive );
    }
  }


  unpackedSize += unpackSets( buffer );


  return unpackedSize;
}


std::set< string > ObjectManagerBase::getPackingExclusionList() const
{
  return m_packingExclusionList;
}


void ObjectManagerBase::excludeWrappersFromPacking( std::set< string > const & wrapperNames )
{
  m_packingExclusionList.insert( wrapperNames.cbegin(), wrapperNames.cend() );
}


localIndex ObjectManagerBase::getNumberOfGhosts() const
{
  localIndex rval = 0;
  for( localIndex i=0; i<size(); ++i )
  {
    if( m_ghostRank[i] > -1 )
    {
      ++rval;
    }
  }
  return rval;
//  return std::count_if( m_ghostRank.begin(), m_ghostRank.end(), [](integer i)->localIndex {return i>-1;} );
}

localIndex ObjectManagerBase::getNumberOfLocalIndices() const
{
  localIndex rval = 0;
  for( localIndex i=0; i<size(); ++i )
  {
    if( m_ghostRank[i] <= -1 )
    {
      ++rval;
    }
  }
  return rval;
  //return std::count_if( m_ghostRank.begin(), m_ghostRank.end(), [](integer i)->localIndex {return i==-1;} );
}

void ObjectManagerBase::setReceiveLists()
{
  for( std::pair< int const, NeighborData > & pair : m_neighborData )
  {
    pair.second.ghostsToReceive().clear();
  }

  for( localIndex a=0; a<size(); ++a )
  {
    if( m_ghostRank[a] > -1 )
    {
      getNeighborData( m_ghostRank[ a ] ).ghostsToReceive().emplace_back( a );
    }
  }
}

integer ObjectManagerBase::splitObject( localIndex const indexToSplit,
                                        int const GEOSX_UNUSED_PARAM( rank ),
                                        localIndex & newIndex )
{
  // if the object index has a zero sized childIndices entry, then this object can be split into two
  // new objects

  if( size()+1 > capacity() )
  {
    reserve( static_cast< localIndex >( size() * m_overAllocationFactor ) );
  }

  // the new indices are tacked on to the end of the arrays
  newIndex = size();
  this->resize( newIndex + 1 );

  // copy the fields
  copyObject( indexToSplit, newIndex );

  if( this->hasExtrinsicData< extrinsicMeshData::ParentIndex >() )
  {
    arrayView1d< localIndex > const & parentIndex = this->getExtrinsicData< extrinsicMeshData::ParentIndex >();
    parentIndex[newIndex] = indexToSplit;
  }

  if( this->hasExtrinsicData< extrinsicMeshData::ChildIndex >() )
  {
    arrayView1d< localIndex > const & childIndex = this->getExtrinsicData< extrinsicMeshData::ChildIndex >();
    childIndex[indexToSplit] = newIndex;
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

void ObjectManagerBase::inheritGhostRankFromParent( std::set< localIndex > const & indices )
{
  arrayView1d< localIndex const > const parentIndex = this->getExtrinsicData< extrinsicMeshData::ParentIndex >();

  for( auto const a : indices )
  {
    m_ghostRank[ a ] = m_ghostRank[ parentIndex[ a ] ];
  }
}


void ObjectManagerBase::copyObject( const localIndex source, const localIndex destination )
{
  for( auto & nameToWrapper: wrappers() )
  {
    WrapperBase * wrapper = nameToWrapper.second;
    if( wrapper->sizedFromParent() )
    {
      wrapper->copy( source, destination );
    }
  }

  for( localIndex i=0; i<m_sets.wrappers().size(); ++i )
  {
    SortedArray< localIndex > & targetSet = m_sets.getReference< SortedArray< localIndex > >( i );

#if !defined(__CUDA_ARCH__)
    targetSet.move( LvArray::MemorySpace::host, true );
#endif

    if( targetSet.count( source ) > 0 )
    {
      targetSet.insert( destination );
    }
  }
}

void ObjectManagerBase::setMaxGlobalIndex()
{
  m_maxGlobalIndex = MpiWrapper::max( m_localMaxGlobalIndex, MPI_COMM_GEOSX );
}

void ObjectManagerBase::cleanUpMap( std::set< localIndex > const & targetIndices,
                                    array1d< SortedArray< localIndex > > & upmap,
                                    arrayView2d< localIndex const > const & downmap )
{
  for( auto const & targetIndex : targetIndices )
  {
    SortedArray< localIndex > eraseList;
    for( auto const & compositeIndex : upmap[targetIndex] )
    {
      bool hasTargetIndex = false;
      for( localIndex a=0; a<downmap.size( 1 ); ++a )
      {
        localIndex const compositeLocalIndex = downmap[compositeIndex][a];
        if( compositeLocalIndex==targetIndex )
        {
          hasTargetIndex=true;
        }
      }
      if( !hasTargetIndex )
      {
        eraseList.insert( compositeIndex );
      }
    }
    for( auto const & val : eraseList )
    {
      upmap[targetIndex].remove( val );
    }
  }
}

void ObjectManagerBase::cleanUpMap( std::set< localIndex > const & targetIndices,
                                    ArrayOfSetsView< localIndex > const & upmap,
                                    arrayView2d< localIndex const > const & downmap )
{
  std::vector< localIndex > eraseList;
  for( localIndex const targetIndex : targetIndices )
  {
    eraseList.clear();
    localIndex pos = 0;
    for( auto const & compositeIndex : upmap[ targetIndex ] )
    {
      bool hasTargetIndex = false;
      for( localIndex a=0; a<downmap.size( 1 ); ++a )
      {
        localIndex const compositeLocalIndex = downmap[compositeIndex][a];
        if( compositeLocalIndex==targetIndex )
        {
          hasTargetIndex=true;
        }
      }

      if( !hasTargetIndex )
      {
        eraseList.emplace_back( pos );
      }

      ++pos;
    }

    localIndex const numUniqueIndices = LvArray::sortedArrayManipulation::makeSortedUnique( eraseList.begin(), eraseList.end() );
    upmap.removeFromSet( targetIndex, eraseList.begin(), eraseList.begin() + numUniqueIndices );
  }
}

void ObjectManagerBase::cleanUpMap( std::set< localIndex > const & targetIndices,
                                    array1d< SortedArray< localIndex > > & upmap,
                                    arrayView1d< arrayView1d< localIndex const > const > const & downmap )
{
  for( auto const & targetIndex : targetIndices )
  {
    SortedArray< localIndex > eraseList;
    for( auto const & compositeIndex : upmap[targetIndex] )
    {
      bool hasTargetIndex = false;
      for( localIndex a=0; a<downmap[compositeIndex].size(); ++a )
      {
        localIndex const compositeLocalIndex = downmap[compositeIndex][a];
        if( compositeLocalIndex==targetIndex )
        {
          hasTargetIndex=true;
        }
      }
      if( !hasTargetIndex )
      {
        eraseList.insert( compositeIndex );
      }
    }
    for( auto const & val : eraseList )
    {
      upmap[targetIndex].remove( val );
    }
  }
}

void ObjectManagerBase::cleanUpMap( std::set< localIndex > const & targetIndices,
                                    ArrayOfSetsView< localIndex > const & upmap,
                                    arrayView1d< arrayView1d< localIndex const > const > const & downmap )
{
  std::vector< localIndex > eraseList;
  for( localIndex const targetIndex : targetIndices )
  {
    eraseList.clear();
    localIndex pos = 0;
    for( localIndex const compositeIndex : upmap[ targetIndex ] )
    {
      bool hasTargetIndex = false;
      for( localIndex a=0; a<downmap[compositeIndex].size(); ++a )
      {
        localIndex const compositeLocalIndex = downmap[compositeIndex][a];
        if( compositeLocalIndex==targetIndex )
        {
          hasTargetIndex=true;
        }
      }

      if( !hasTargetIndex )
      {
        eraseList.emplace_back( pos );
      }

      ++pos;
    }

    localIndex const numUniqueIndices = LvArray::sortedArrayManipulation::makeSortedUnique( eraseList.begin(), eraseList.end() );
    upmap.removeFromSet( targetIndex, eraseList.begin(), eraseList.begin() + numUniqueIndices );
  }
}

void ObjectManagerBase::cleanUpMap( std::set< localIndex > const & targetIndices,
                                    ArrayOfSetsView< localIndex > const & upmap,
                                    ArrayOfArraysView< localIndex const > const & downmap )
{
  std::vector< localIndex > eraseList;
  for( localIndex const targetIndex : targetIndices )
  {
    eraseList.clear();
    for( localIndex const compositeIndex : upmap[ targetIndex ] )
    {
      bool hasTargetIndex = false;
      for( localIndex const compositeLocalIndex : downmap[ compositeIndex ] )
      {
        if( compositeLocalIndex == targetIndex )
        {
          hasTargetIndex = true;
        }
      }

      if( !hasTargetIndex )
      {
        eraseList.emplace_back( compositeIndex );
      }
    }

    localIndex const numUniqueIndices = LvArray::sortedArrayManipulation::makeSortedUnique( eraseList.begin(), eraseList.end() );
    upmap.removeFromSet( targetIndex, eraseList.begin(), eraseList.begin() + numUniqueIndices );
  }
}


void ObjectManagerBase::enforceStateFieldConsistencyPostTopologyChange( std::set< localIndex > const & targetIndices )
{
  arrayView1d< localIndex const > const childFaceIndices = getExtrinsicData< extrinsicMeshData::ChildIndex >();

  for( localIndex const targetIndex : targetIndices )
  {
    localIndex const childIndex = childFaceIndices[targetIndex];
    if( childIndex != -1 )
    {
      this->m_isExternal[targetIndex] = m_isExternal[childIndex];
    }
  }
}


void ObjectManagerBase::moveSets( LvArray::MemorySpace const targetSpace )
{
  m_sets.forWrappers< SortedArray< localIndex > >( [&] ( auto & wrapper )
  {
    SortedArray< localIndex > & set = wrapper.reference();
    set.move( targetSpace );
  } );
}


} /* namespace geosx */
