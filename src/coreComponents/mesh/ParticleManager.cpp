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

#include <map>
#include <vector>

#include "ParticleManager.hpp"

#include "common/TimingMacros.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "generators/ParticleBlockManager.hpp"
#include "mesh/MeshManager.hpp"
#include "schema/schemaUtilities.hpp"

namespace geosx
{

using namespace dataRepository;

class SpatialPartition;

// *********************************************************************************************************************
/**
 * @return
 */

ParticleManager::ParticleManager( string const & name, Group * const parent ):
  ObjectManagerBase( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
  this->registerGroup< Group >( ParticleManager::groupKeyStruct::particleRegionsGroup() );
}

ParticleManager::~ParticleManager()
{
  // TODO Auto-generated destructor stub
}

localIndex ParticleManager::numParticleBlocks() const
{
  localIndex numParticleBlocks = 0;
  this->forParticleSubRegions< ParticleSubRegionBase >( [&]( ParticleSubRegionBase const & )
  {
    numParticleBlocks += 1;
  } );
  return numParticleBlocks;
}

void ParticleManager::resize( integer_array const & numParticles,
                                   string_array const & regionNames,
                                   string_array const & GEOSX_UNUSED_PARAM( particleTypes ) )
{
  localIndex const n_regions = LvArray::integerConversion< localIndex >( regionNames.size());
  for( localIndex reg=0; reg<n_regions; ++reg )
  {
    this->getRegion( reg ).resize( numParticles[reg] );
  }
}

void ParticleManager::setMaxGlobalIndex()
{
  forParticleSubRegions< ParticleSubRegionBase >( [this] ( ParticleSubRegionBase const & subRegion )
  {
    m_localMaxGlobalIndex = std::max( m_localMaxGlobalIndex, subRegion.maxGlobalIndex() );
  } );

  MpiWrapper::allReduce( &m_localMaxGlobalIndex,
                         &m_maxGlobalIndex,
                         1,
                         MPI_MAX,
                         MPI_COMM_GEOSX );
}



Group * ParticleManager::createChild( string const & childKey, string const & childName )
{
  GEOSX_ERROR_IF( !(CatalogInterface::hasKeyName( childKey )),
                  "KeyName ("<<childKey<<") not found in ObjectManager::Catalog" );
  GEOSX_LOG_RANK_0( "Adding Object " << childKey << " named " << childName << " from ObjectManager::Catalog." );

  Group & particleRegions = this->getGroup( ParticleManager::groupKeyStruct::particleRegionsGroup() );
  return &particleRegions.registerGroup( childName,
                                        CatalogInterface::factory( childKey, childName, &particleRegions ) );

}

void ParticleManager::expandObjectCatalogs()
{
  ObjectManagerBase::CatalogInterface::CatalogType const & catalog = ObjectManagerBase::getCatalog();
  for( ObjectManagerBase::CatalogInterface::CatalogType::const_iterator iter = catalog.begin();
       iter!=catalog.end();
       ++iter )
  {
    string const key = iter->first;
    if( key.find( "ParticleRegion" ) != string::npos )
    {
      this->createChild( key, key );
    }
  }
}


void ParticleManager::setSchemaDeviations( xmlWrapper::xmlNode schemaRoot,
                                                xmlWrapper::xmlNode schemaParent,
                                                integer documentationType )
{
  xmlWrapper::xmlNode targetChoiceNode = schemaParent.child( "xsd:choice" );
  if( targetChoiceNode.empty() )
  {
    targetChoiceNode = schemaParent.prepend_child( "xsd:choice" );
    targetChoiceNode.append_attribute( "minOccurs" ) = "0";
    targetChoiceNode.append_attribute( "maxOccurs" ) = "unbounded";
  }

  std::set< string > names;
  this->forParticleRegions( [&]( ParticleRegionBase & particleRegion )
  {
    names.insert( particleRegion.getName() );
  } );

  for( string const & name: names )
  {
    schemaUtilities::SchemaConstruction( getRegion( name ), schemaRoot, targetChoiceNode, documentationType );
  }
}

void ParticleManager::generateMesh( Group & particleBlockManager )
{
  this->forParticleRegions< ParticleRegion >( [&]( auto & particleRegion )
  {
    particleRegion.generateMesh( particleBlockManager.getGroup( keys::particleBlocks ) );
  } );
}

//void ParticleManager::generateAggregates( FaceManager const & faceManager, NodeManager const & nodeManager )
//{
//  this->forParticleRegions< ParticleRegion >( [&]( ParticleRegion & particleRegion )
//  {
//    particleRegion.generateAggregates( faceManager, nodeManager );
//  } );
//}

int ParticleManager::PackSize( string_array const & wrapperNames,
                                    ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackPrivate< false >( junk, wrapperNames, packList );
}

int ParticleManager::Pack( buffer_unit_type * & buffer,
                                string_array const & wrapperNames,
                                ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  return PackPrivate< true >( buffer, wrapperNames, packList );
}

template< bool DOPACK >
int
ParticleManager::PackPrivate( buffer_unit_type * & buffer,
                                   string_array const & wrapperNames,
                                   ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  int packedSize = 0;

//  packedSize += Group::Pack( buffer, wrapperNames, {}, 0, 0);

  packedSize += bufferOps::Pack< DOPACK >( buffer, this->getName() );
  packedSize += bufferOps::Pack< DOPACK >( buffer, numRegions() );

  parallelDeviceEvents events;
  for( typename dataRepository::indexType kReg=0; kReg<numRegions(); ++kReg )
  {
    ParticleRegionBase const & particleRegion = getRegion( kReg );
    packedSize += bufferOps::Pack< DOPACK >( buffer, particleRegion.getName() );

    packedSize += bufferOps::Pack< DOPACK >( buffer, particleRegion.numSubRegions() );

    particleRegion.forParticleSubRegionsIndex< ParticleSubRegionBase >(
      [&]( localIndex const esr, ParticleSubRegionBase const & subRegion )
    {
      packedSize += bufferOps::Pack< DOPACK >( buffer, subRegion.getName() );

      arrayView1d< localIndex const > const particleList = packList[kReg][esr];
      if( DOPACK )
      {
        packedSize += subRegion.pack( buffer, wrapperNames, particleList, 0, false, events );
      }
      else
      {
        packedSize += subRegion.packSize( wrapperNames, particleList, 0, false, events );
      }
    } );
  }

  waitAllDeviceEvents( events );
  return packedSize;
}


int ParticleManager::Unpack( buffer_unit_type const * & buffer,
                                  ElementViewAccessor< arrayView1d< localIndex > > & packList )
{
  return unpackPrivate( buffer, packList );
}

int ParticleManager::Unpack( buffer_unit_type const * & buffer,
                                  ElementReferenceAccessor< array1d< localIndex > > & packList )
{
  return unpackPrivate( buffer, packList );
}

template< typename T >
int ParticleManager::unpackPrivate( buffer_unit_type const * & buffer,
                                         T & packList )
{
  int unpackedSize = 0;

  string name;
  unpackedSize += bufferOps::Unpack( buffer, name );

  GEOSX_ERROR_IF( name != this->getName(), "Unpacked name (" << name << ") does not equal object name (" << this->getName() << ")" );

  localIndex numRegionsRead;
  unpackedSize += bufferOps::Unpack( buffer, numRegionsRead );

  parallelDeviceEvents events;
  for( localIndex kReg=0; kReg<numRegionsRead; ++kReg )
  {
    string regionName;
    unpackedSize += bufferOps::Unpack( buffer, regionName );

    ParticleRegionBase & particleRegion = getRegion( regionName );

    localIndex numSubRegionsRead;
    unpackedSize += bufferOps::Unpack( buffer, numSubRegionsRead );
    particleRegion.forParticleSubRegionsIndex< ParticleSubRegionBase >(
      [&]( localIndex const esr, ParticleSubRegionBase & subRegion )
    {
      string subRegionName;
      unpackedSize += bufferOps::Unpack( buffer, subRegionName );

      /// THIS IS WRONG??
      arrayView1d< localIndex > & particleList = packList[kReg][esr];

      unpackedSize += subRegion.unpack( buffer, particleList, 0, false, events );
    } );
  }

  waitAllDeviceEvents( events );
  return unpackedSize;
}

int ParticleManager::PackGlobalMapsSize( ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackGlobalMapsPrivate< false >( junk, packList );
}

int ParticleManager::PackGlobalMaps( buffer_unit_type * & buffer,
                                          ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  return PackGlobalMapsPrivate< true >( buffer, packList );
}
template< bool DOPACK >
int
ParticleManager::PackGlobalMapsPrivate( buffer_unit_type * & buffer,
                                             ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  int packedSize = 0;

  packedSize += bufferOps::Pack< DOPACK >( buffer, numRegions() );

  for( typename dataRepository::indexType kReg=0; kReg<numRegions(); ++kReg )
  {
    ParticleRegionBase const & particleRegion = getRegion( kReg );
    packedSize += bufferOps::Pack< DOPACK >( buffer, particleRegion.getName() );

    packedSize += bufferOps::Pack< DOPACK >( buffer, particleRegion.numSubRegions() );
    particleRegion.forParticleSubRegionsIndex< ParticleSubRegionBase >(
      [&]( localIndex const esr, ParticleSubRegionBase const & subRegion )
    {
      packedSize += bufferOps::Pack< DOPACK >( buffer, subRegion.getName() );

      arrayView1d< localIndex const > const particleList = packList[kReg][esr];
      if( DOPACK )
      {
        packedSize += subRegion.packGlobalMaps( buffer, particleList, 0 );
      }
      else
      {
        packedSize += subRegion.packGlobalMapsSize( particleList, 0 );
      }
    } );
  }

  return packedSize;
}


int
ParticleManager::UnpackGlobalMaps( buffer_unit_type const * & buffer,
                                        ElementViewAccessor< ReferenceWrapper< localIndex_array > > & packList )
{
  int unpackedSize = 0;

  localIndex numRegionsRead;
  unpackedSize += bufferOps::Unpack( buffer, numRegionsRead );

  packList.resize( numRegionsRead );
  for( localIndex kReg=0; kReg<numRegionsRead; ++kReg )
  {
    string regionName;
    unpackedSize += bufferOps::Unpack( buffer, regionName );

    ParticleRegionBase & particleRegion = getRegion( regionName );

    localIndex numSubRegionsRead;
    unpackedSize += bufferOps::Unpack( buffer, numSubRegionsRead );
    packList[kReg].resize( numSubRegionsRead );
    particleRegion.forParticleSubRegionsIndex< ParticleSubRegionBase >(
      [&]( localIndex const esr, ParticleSubRegionBase & subRegion )
    {
      string subRegionName;
      unpackedSize += bufferOps::Unpack( buffer, subRegionName );

      /// THIS IS WRONG
      localIndex_array & particleList = packList[kReg][esr].get();

      unpackedSize += subRegion.unpackGlobalMaps( buffer, particleList, 0 );
    } );
  }

  return unpackedSize;
}



int ParticleManager::PackUpDownMapsSize( ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsPrivate< false >( junk, packList );
}
int ParticleManager::PackUpDownMapsSize( ElementReferenceAccessor< array1d< localIndex > > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsPrivate< false >( junk, packList );
}

int ParticleManager::PackUpDownMaps( buffer_unit_type * & buffer,
                                          ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  return packUpDownMapsPrivate< true >( buffer, packList );
}
int ParticleManager::PackUpDownMaps( buffer_unit_type * & buffer,
                                          ElementReferenceAccessor< array1d< localIndex > > const & packList ) const
{
  return packUpDownMapsPrivate< true >( buffer, packList );
}

template< bool DOPACK, typename T >
int
ParticleManager::packUpDownMapsPrivate( buffer_unit_type * & buffer,
                                             T const & packList ) const
{
  int packedSize = 0;

  packedSize += bufferOps::Pack< DOPACK >( buffer, numRegions() );

  for( typename dataRepository::indexType kReg=0; kReg<numRegions(); ++kReg )
  {
    ParticleRegionBase const & particleRegion = getRegion( kReg );
    packedSize += bufferOps::Pack< DOPACK >( buffer, particleRegion.getName() );

    packedSize += bufferOps::Pack< DOPACK >( buffer, particleRegion.numSubRegions() );
    particleRegion.forParticleSubRegionsIndex< ParticleSubRegionBase >(
      [&]( localIndex const esr, ParticleSubRegionBase const & subRegion )
    {
      packedSize += bufferOps::Pack< DOPACK >( buffer, subRegion.getName() );

      arrayView1d< localIndex > const particleList = packList[kReg][esr];
      if( DOPACK )
      {
        packedSize += subRegion.packUpDownMaps( buffer, particleList );
      }
      else
      {
        packedSize += subRegion.packUpDownMapsSize( particleList );
      }
    } );
  }

  return packedSize;
}
//template int
//ParticleManager::
//PackUpDownMapsPrivate<true>( buffer_unit_type * & buffer,
//                             ElementViewAccessor<arrayView1d<localIndex>> const & packList ) const;
//template int
//ParticleManager::
//PackUpDownMapsPrivate<false>( buffer_unit_type * & buffer,
//                             ElementViewAccessor<arrayView1d<localIndex>> const & packList ) const;


int
ParticleManager::UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                        ElementReferenceAccessor< localIndex_array > & packList,
                                        bool const overwriteMap )
{
  int unpackedSize = 0;

  localIndex numRegionsRead;
  unpackedSize += bufferOps::Unpack( buffer, numRegionsRead );

  for( localIndex kReg=0; kReg<numRegionsRead; ++kReg )
  {
    string regionName;
    unpackedSize += bufferOps::Unpack( buffer, regionName );

    ParticleRegionBase & particleRegion = getRegion( regionName );

    localIndex numSubRegionsRead;
    unpackedSize += bufferOps::Unpack( buffer, numSubRegionsRead );
    particleRegion.forParticleSubRegionsIndex< ParticleSubRegionBase >(
      [&]( localIndex const kSubReg, ParticleSubRegionBase & subRegion )
    {
      string subRegionName;
      unpackedSize += bufferOps::Unpack( buffer, subRegionName );

      /// THIS IS WRONG
      localIndex_array & particleList = packList[kReg][kSubReg];
      unpackedSize += subRegion.unpackUpDownMaps( buffer, particleList, false, overwriteMap );
    } );
  }

  return unpackedSize;
}




REGISTER_CATALOG_ENTRY( ObjectManagerBase, ParticleManager, string const &, Group * const )
}
