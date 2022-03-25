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
 * @file DomainPartition.cpp
 */

#include "DomainPartition.hpp"

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"



namespace geosx
{
using namespace dataRepository;

DomainPartition::DomainPartition( string const & name,
                                  Group * const parent ):
  Group( name, parent )
{
  this->registerWrapper( "Neighbors", &m_neighbors ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( false );

  this->registerWrapper< SpatialPartition, PartitionBase >( keys::partitionManager ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( false );

  registerGroup( groupKeys.meshBodies );
  registerGroup< constitutive::ConstitutiveManager >( groupKeys.constitutiveManager );
}


DomainPartition::~DomainPartition()
{}

void DomainPartition::initializationOrder( string_array & order )
{
  SortedArray< string > usedNames;
  {
    order.emplace_back( keys::ConstitutiveManager );
    usedNames.insert( keys::ConstitutiveManager );
  }

  {
    order.emplace_back( string( groupKeysStruct::meshBodiesString() ) );
    usedNames.insert( groupKeysStruct::meshBodiesString() );
  }


  for( auto const & subGroup : this->getSubGroups() )
  {
    if( usedNames.count( subGroup.first ) == 0 )
    {
      order.emplace_back( subGroup.first );
    }
  }
}

void DomainPartition::setupCommunications( bool use_nonblocking )
{
  GEOSX_MARK_FUNCTION;

#if defined(GEOSX_USE_MPI)
  if( m_metisNeighborList.empty() )
  {
    PartitionBase & partition1 = getReference< PartitionBase >( keys::partitionManager );
    SpatialPartition & partition = dynamic_cast< SpatialPartition & >(partition1);

    //get communicator, rank, and coordinates
    MPI_Comm cartcomm;
    {
      int reorder = 0;
      MpiWrapper::cartCreate( MPI_COMM_GEOSX, 3, partition.m_Partitions.data(), partition.m_Periodic.data(), reorder, &cartcomm );
      GEOSX_ERROR_IF( cartcomm == MPI_COMM_NULL, "Fail to run MPI_Cart_create and establish communications" );
    }
    int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
    int nsdof = 3;

    MpiWrapper::cartCoords( cartcomm, rank, nsdof, partition.m_coords.data() );

    int ncoords[3];
    addNeighbors( 0, cartcomm, ncoords );

    MpiWrapper::commFree( cartcomm );
  }
  else
  {
    for( integer const neighborRank : m_metisNeighborList )
    {
      m_neighbors.emplace_back( NeighborCommunicator( neighborRank ) );
    }
  }

  // Create an array of the first neighbors.
  array1d< int > firstNeighborRanks;
  for( NeighborCommunicator const & neighbor : m_neighbors )
  {
    firstNeighborRanks.emplace_back( neighbor.neighborRank() );
  }

  int neighborsTag = 54;

  // Send this list of neighbors to all neighbors.
  std::vector< MPI_Request > requests( m_neighbors.size() );
  for( std::size_t i = 0; i < m_neighbors.size(); ++i )
  {
    MpiWrapper::iSend( firstNeighborRanks.toViewConst(), m_neighbors[ i ].neighborRank(), neighborsTag, MPI_COMM_GEOSX, &requests[ i ] );
  }

  // This set will contain the second (neighbor of) neighbors ranks.
  std::set< int > secondNeighborRanks;

  array1d< int > neighborOfNeighborRanks;
  for( std::size_t i = 0; i < m_neighbors.size(); ++i )
  {
    MpiWrapper::recv( neighborOfNeighborRanks, m_neighbors[ i ].neighborRank(), neighborsTag, MPI_COMM_GEOSX, MPI_STATUS_IGNORE );

    // Insert the neighbors of the current neighbor into the set of second neighbors.
    secondNeighborRanks.insert( neighborOfNeighborRanks.begin(), neighborOfNeighborRanks.end() );
  }

  // Remove yourself and all the first neighbors from the second neighbors.
  secondNeighborRanks.erase( MpiWrapper::commRank() );
  for( NeighborCommunicator const & neighbor : m_neighbors )
  {
    secondNeighborRanks.erase( neighbor.neighborRank() );
  }

  for( integer const neighborRank : secondNeighborRanks )
  {
    m_neighbors.emplace_back( NeighborCommunicator( neighborRank ) );
  }

  MpiWrapper::waitAll( requests.size(), requests.data(), MPI_STATUSES_IGNORE );

#endif

  MeshLevel & meshLevel = getMeshBody( 0 ).getMeshLevel( 0 );

  for( NeighborCommunicator const & neighbor : m_neighbors )
  {
    neighbor.addNeighborGroupToMesh( meshLevel );
  }

  NodeManager & nodeManager = meshLevel.getNodeManager();
  FaceManager & faceManager = meshLevel.getFaceManager();
  EdgeManager & edgeManager = meshLevel.getEdgeManager();

  nodeManager.setMaxGlobalIndex();

  CommunicationTools::getInstance().assignGlobalIndices( faceManager, nodeManager, m_neighbors );

  CommunicationTools::getInstance().assignGlobalIndices( edgeManager, nodeManager, m_neighbors );

  CommunicationTools::getInstance().findMatchedPartitionBoundaryObjects( faceManager,
                                                                         m_neighbors );

  CommunicationTools::getInstance().findMatchedPartitionBoundaryObjects( nodeManager,
                                                                         m_neighbors );

  CommunicationTools::getInstance().setupGhosts( meshLevel, m_neighbors, use_nonblocking );


  // TODO Consider moving these calls in ProblemManager::generateMesh()
  // after ghosts are generated, ie after the call to setupCommunications())
  faceManager.sortAllFaceNodes( nodeManager, meshLevel.getElemManager() );
  faceManager.computeGeometry( nodeManager );
}

void DomainPartition::addNeighbors( const unsigned int idim,
                                    MPI_Comm & cartcomm,
                                    int * ncoords )
{
  PartitionBase & partition1 = getReference< PartitionBase >( keys::partitionManager );
  SpatialPartition & partition = dynamic_cast< SpatialPartition & >(partition1);

  if( idim == nsdof )
  {
    bool me = true;
    for( int i = 0; i < nsdof; i++ )
    {
      if( ncoords[i] != partition.m_coords( i ))
      {
        me = false;
        break;
      }
    }
    int const neighborRank = MpiWrapper::cartRank( cartcomm, ncoords );
    if( !me && !std::any_of( m_neighbors.begin(), m_neighbors.end(), [=]( NeighborCommunicator const & nn ) { return nn.neighborRank( ) == neighborRank; } ) )
    {
      m_neighbors.emplace_back( NeighborCommunicator( neighborRank ) );
    }
  }
  else
  {
    const int dim = partition.m_Partitions( LvArray::integerConversion< localIndex >( idim ));
    const bool periodic = partition.m_Periodic( LvArray::integerConversion< localIndex >( idim ));
    for( int i = -1; i < 2; i++ )
    {
      ncoords[idim] = partition.m_coords( LvArray::integerConversion< localIndex >( idim )) + i;
      bool ok = true;
      if( periodic )
      {
        if( ncoords[idim] < 0 )
          ncoords[idim] = dim - 1;
        else if( ncoords[idim] >= dim )
          ncoords[idim] = 0;
      }
      else
      {
        ok = ncoords[idim] >= 0 && ncoords[idim] < dim;
      }
      if( ok )
      {
        addNeighbors( idim + 1, cartcomm, ncoords );
      }
    }
  }
}

} /* namespace geosx */
