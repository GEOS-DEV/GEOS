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
#include "ParticleBlockManager.hpp"
#include "mesh/MeshManager.hpp"
#include "schema/schemaUtilities.hpp"

namespace geosx
{
using namespace dataRepository;

// *********************************************************************************************************************
/**
 * @return
 */
//START_SPHINX_REFPOS_REG

ParticleManager::ParticleManager( string const & name, Group * const parent ):
  ObjectManagerBase( name, parent ),
  m_referencePosition( 0, 3 )
{
  registerWrapper( viewKeyStruct::referencePositionString(), &m_referencePosition );
  //END_SPHINX_REFPOS_REG
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
  this->forParticleSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase const & )
  {
    numParticleBlocks += 1;
  } );
  return numParticleBlocks;
}

void ParticleManager::resize( integer_array const & numElements,
                                   string_array const & regionNames,
                                   string_array const & GEOSX_UNUSED_PARAM( elementTypes ) )
{
  localIndex const n_regions = LvArray::integerConversion< localIndex >( regionNames.size());
  for( localIndex reg=0; reg<n_regions; ++reg )
  {
    this->getRegion( reg ).resize( numElements[reg] );
  }
}

void ParticleManager::setMaxGlobalIndex()
{
  forParticleSubRegions< ElementSubRegionBase >( [this] ( ElementSubRegionBase const & subRegion )
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
  GEOSX_LOG_RANK_0( "Adding Object " << childKey<<" named "<< childName<<" from ObjectManager::Catalog." );

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
  this->forParticleRegions( [&]( ElementRegionBase & particleRegion )
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
  this->forParticleRegions< ParticleRegion, SurfaceElementRegion >( [&]( auto & particleRegion )
  {
    particleRegion.generateMesh( particleBlockManager.getGroup( keys::particleBlocks ) );
  } );
}

void ParticleManager::generateCellToEdgeMaps( FaceManager const & faceManager )
{
  /*
   * Create cell to edges map
   * I use the existing maps from cells to faces and from faces to edges.
   */
  localIndex faceIndex, edgeIndex;
  int count = 0;
  bool isUnique = true;

  this->forParticleSubRegions< ParticleSubRegion >( [&]( ParticleSubRegion & subRegion )
  {
    FixedOneToManyRelation & cellToEdges = subRegion.edgeList();
    FixedOneToManyRelation const & cellToFaces = subRegion.faceList();
    InterObjectRelation< ArrayOfArrays< localIndex > > const & faceToEdges = faceManager.edgeList();

    //loop over the cells
    for( localIndex kc = 0; kc < subRegion.size(); kc++ )
    {
      count = 0;
      for( localIndex kf = 0; kf < subRegion.numFacesPerElement(); kf++ )
      {
        // loop over edges of each face
        faceIndex = cellToFaces[kc][kf];
        for( localIndex ke = 0; ke < faceToEdges.sizeOfArray( faceIndex ); ke++ )
        {
          isUnique = true;
          edgeIndex = faceToEdges[faceIndex][ke];

          //loop over edges that have already been added to the element.
          for( localIndex kec = 0; kec < count+1; kec++ )
          {
            // make sure that the edge has not been counted yet
            if( cellToEdges( kc, kec ) == edgeIndex )
            {
              isUnique = false;
              break;
            }
          }
          if( isUnique )
          {
            cellToEdges( kc, count ) = edgeIndex;
            count++;
          }

        } // end edge loop
      } // end face loop
    } // end cell loop
  } );
}

//void ParticleManager::generateAggregates( FaceManager const & faceManager, NodeManager const & nodeManager )
//{
//  this->forParticleRegions< ParticleRegion >( [&]( ParticleRegion & particleRegion )
//  {
//    particleRegion.generateAggregates( faceManager, nodeManager );
//  } );
//}

void ParticleManager::generateWells( MeshManager & meshManager,
                                          MeshLevel & meshLevel )
{
  NodeManager & nodeManager = meshLevel.getNodeManager();

  // get the offsets to construct local-to-global maps for well nodes and elements
  nodeManager.setMaxGlobalIndex();
  globalIndex const nodeOffsetGlobal = nodeManager.maxGlobalIndex() + 1;
  localIndex const elemOffsetLocal  = this->getNumberOfElements();
  globalIndex const elemOffsetGlobal = MpiWrapper::sum( elemOffsetLocal );

  globalIndex wellElemCount = 0;
  globalIndex wellNodeCount = 0;

  // construct the wells one by one
  forParticleRegions< WellElementRegion >( [&]( WellElementRegion & wellRegion )
  {

    // get the global well geometry from the well generator
    string const generatorName = wellRegion.getWellGeneratorName();
    InternalWellGenerator const & wellGeometry =
      meshManager.getGroup< InternalWellGenerator >( generatorName );

    // generate the local data (well elements, nodes, perforations) on this well
    // note: each MPI rank knows the global info on the entire well (constructed earlier in InternalWellGenerator)
    // so we only need node and element offsets to construct the local-to-global maps in each wellElemSubRegion
    wellRegion.generateWell( meshLevel, wellGeometry, nodeOffsetGlobal + wellNodeCount, elemOffsetGlobal + wellElemCount );

    // increment counters with global number of nodes and elements
    wellElemCount += wellGeometry.getNumElements();
    wellNodeCount += wellGeometry.getNumNodes();

    string const & subRegionName = wellRegion.getSubRegionName();
    WellElementSubRegion &
    subRegion = wellRegion.getGroup( ElementRegionBase::viewKeyStruct::elementSubRegions() )
                  .getGroup< WellElementSubRegion >( subRegionName );

    globalIndex const numWellElemsGlobal = MpiWrapper::sum( subRegion.size() );

    GEOSX_ERROR_IF( numWellElemsGlobal != wellGeometry.getNumElements(),
                    "Invalid partitioning in well " << subRegionName );

  } );

  // communicate to rebuild global node info since we modified global ordering
  nodeManager.setMaxGlobalIndex();
}

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
    ElementRegionBase const & particleRegion = getRegion( kReg );
    packedSize += bufferOps::Pack< DOPACK >( buffer, particleRegion.getName() );

    packedSize += bufferOps::Pack< DOPACK >( buffer, particleRegion.numSubRegions() );

    particleRegion.forElementSubRegionsIndex< ElementSubRegionBase >(
      [&]( localIndex const esr, ElementSubRegionBase const & subRegion )
    {
      packedSize += bufferOps::Pack< DOPACK >( buffer, subRegion.getName() );

      arrayView1d< localIndex const > const elemList = packList[kReg][esr];
      if( DOPACK )
      {
        packedSize += subRegion.pack( buffer, wrapperNames, elemList, 0, false, events );
      }
      else
      {
        packedSize += subRegion.packSize( wrapperNames, elemList, 0, false, events );
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

    ElementRegionBase & particleRegion = getRegion( regionName );

    localIndex numSubRegionsRead;
    unpackedSize += bufferOps::Unpack( buffer, numSubRegionsRead );
    particleRegion.forElementSubRegionsIndex< ElementSubRegionBase >(
      [&]( localIndex const esr, ElementSubRegionBase & subRegion )
    {
      string subRegionName;
      unpackedSize += bufferOps::Unpack( buffer, subRegionName );

      /// THIS IS WRONG??
      arrayView1d< localIndex > & elemList = packList[kReg][esr];

      unpackedSize += subRegion.unpack( buffer, elemList, 0, false, events );
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
    ElementRegionBase const & particleRegion = getRegion( kReg );
    packedSize += bufferOps::Pack< DOPACK >( buffer, particleRegion.getName() );

    packedSize += bufferOps::Pack< DOPACK >( buffer, particleRegion.numSubRegions() );
    particleRegion.forElementSubRegionsIndex< ElementSubRegionBase >(
      [&]( localIndex const esr, ElementSubRegionBase const & subRegion )
    {
      packedSize += bufferOps::Pack< DOPACK >( buffer, subRegion.getName() );

      arrayView1d< localIndex const > const elemList = packList[kReg][esr];
      if( DOPACK )
      {
        packedSize += subRegion.packGlobalMaps( buffer, elemList, 0 );
      }
      else
      {
        packedSize += subRegion.packGlobalMapsSize( elemList, 0 );
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

    ElementRegionBase & particleRegion = getRegion( regionName );

    localIndex numSubRegionsRead;
    unpackedSize += bufferOps::Unpack( buffer, numSubRegionsRead );
    packList[kReg].resize( numSubRegionsRead );
    particleRegion.forElementSubRegionsIndex< ElementSubRegionBase >(
      [&]( localIndex const esr, ElementSubRegionBase & subRegion )
    {
      string subRegionName;
      unpackedSize += bufferOps::Unpack( buffer, subRegionName );

      /// THIS IS WRONG
      localIndex_array & elemList = packList[kReg][esr].get();

      unpackedSize += subRegion.unpackGlobalMaps( buffer, elemList, 0 );
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
    ElementRegionBase const & particleRegion = getRegion( kReg );
    packedSize += bufferOps::Pack< DOPACK >( buffer, particleRegion.getName() );

    packedSize += bufferOps::Pack< DOPACK >( buffer, particleRegion.numSubRegions() );
    particleRegion.forElementSubRegionsIndex< ElementSubRegionBase >(
      [&]( localIndex const esr, ElementSubRegionBase const & subRegion )
    {
      packedSize += bufferOps::Pack< DOPACK >( buffer, subRegion.getName() );

      arrayView1d< localIndex > const elemList = packList[kReg][esr];
      if( DOPACK )
      {
        packedSize += subRegion.packUpDownMaps( buffer, elemList );
      }
      else
      {
        packedSize += subRegion.packUpDownMapsSize( elemList );
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

    ElementRegionBase & particleRegion = getRegion( regionName );

    localIndex numSubRegionsRead;
    unpackedSize += bufferOps::Unpack( buffer, numSubRegionsRead );
    particleRegion.forElementSubRegionsIndex< ElementSubRegionBase >(
      [&]( localIndex const kSubReg, ElementSubRegionBase & subRegion )
    {
      string subRegionName;
      unpackedSize += bufferOps::Unpack( buffer, subRegionName );

      /// THIS IS WRONG
      localIndex_array & elemList = packList[kReg][kSubReg];
      unpackedSize += subRegion.unpackUpDownMaps( buffer, elemList, false, overwriteMap );
    } );
  }

  return unpackedSize;
}

int ParticleManager::packFracturedElementsSize( ElementViewAccessor< arrayView1d< localIndex > > const & packList,
                                                     string const fractureRegionName ) const
{
  buffer_unit_type * junk = nullptr;
  return packFracturedElementsPrivate< false >( junk, packList, fractureRegionName );
}

int ParticleManager::packFracturedElements( buffer_unit_type * & buffer,
                                                 ElementViewAccessor< arrayView1d< localIndex > > const & packList,
                                                 string const fractureRegionName ) const
{
  return packFracturedElementsPrivate< true >( buffer, packList, fractureRegionName );
}
template< bool DOPACK >
int
ParticleManager::packFracturedElementsPrivate( buffer_unit_type * & buffer,
                                                    ElementViewAccessor< arrayView1d< localIndex > > const & packList,
                                                    string const fractureRegionName ) const
{
  int packedSize = 0;

  SurfaceElementRegion const & embeddedSurfaceRegion =
    this->getRegion< SurfaceElementRegion >( fractureRegionName );
  EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion =
    embeddedSurfaceRegion.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

  arrayView1d< globalIndex const > const embeddedSurfacesLocalToGlobal =
    embeddedSurfaceSubRegion.localToGlobalMap();

  packedSize += bufferOps::Pack< DOPACK >( buffer, numRegions() );

  for( typename dataRepository::indexType kReg=0; kReg<numRegions(); ++kReg )
  {
    ElementRegionBase const & particleRegion = getRegion( kReg );
    packedSize += bufferOps::Pack< DOPACK >( buffer, particleRegion.getName() );

    packedSize += bufferOps::Pack< DOPACK >( buffer, particleRegion.numSubRegions() );
    particleRegion.forElementSubRegionsIndex< ParticleSubRegion >(
      [&]( localIndex const esr, ParticleSubRegion const & subRegion )
    {
      packedSize += bufferOps::Pack< DOPACK >( buffer, subRegion.getName() );

      arrayView1d< localIndex const > const elemList = packList[kReg][esr];
      if( DOPACK )
      {
        packedSize += subRegion.packFracturedElements( buffer, elemList, embeddedSurfacesLocalToGlobal );
      }
      else
      {
        packedSize += subRegion.packFracturedElementsSize( elemList, embeddedSurfacesLocalToGlobal );
      }
    } );
  }

  return packedSize;
}

int
ParticleManager::unpackFracturedElements( buffer_unit_type const * & buffer,
                                               ElementReferenceAccessor< localIndex_array > & packList,
                                               string const fractureRegionName )
{
  int unpackedSize = 0;

  SurfaceElementRegion & embeddedSurfaceRegion =
    this->getRegion< SurfaceElementRegion >( fractureRegionName );
  EmbeddedSurfaceSubRegion & embeddedSurfaceSubRegion =
    embeddedSurfaceRegion.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

  unordered_map< globalIndex, localIndex > embeddedSurfacesGlobalToLocal =
    embeddedSurfaceSubRegion.globalToLocalMap();

  localIndex numRegionsRead;
  unpackedSize += bufferOps::Unpack( buffer, numRegionsRead );

  for( localIndex kReg=0; kReg<numRegionsRead; ++kReg )
  {
    string regionName;
    unpackedSize += bufferOps::Unpack( buffer, regionName );

    ElementRegionBase & particleRegion = getRegion( regionName );

    localIndex numSubRegionsRead;
    unpackedSize += bufferOps::Unpack( buffer, numSubRegionsRead );
    particleRegion.forElementSubRegionsIndex< ParticleSubRegion >(
      [&]( localIndex const kSubReg, ParticleSubRegion & subRegion )
    {
      string subRegionName;
      unpackedSize += bufferOps::Unpack( buffer, subRegionName );

      /// THIS IS WRONG
      localIndex_array & elemList = packList[kReg][kSubReg];
      unpackedSize += subRegion.unpackFracturedElements( buffer, elemList, embeddedSurfacesGlobalToLocal );
    } );
  }

  return unpackedSize;
}


REGISTER_CATALOG_ENTRY( ObjectManagerBase, ParticleManager, string const &, Group * const )
}
