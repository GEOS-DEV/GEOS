/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
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
#include "fileIO/silo/SiloFile.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "mpiCommunications/SpatialPartition.hpp"



namespace geosx
{
using namespace dataRepository;

DomainPartition::DomainPartition( std::string const & name,
                                  Group * const parent ):
  Group( name, parent )
{
  this->registerWrapper< array1d<NeighborCommunicator> >(viewKeys.neighbors)->setRestartFlags( RestartFlags::NO_WRITE );
  this->registerWrapper<SpatialPartition,PartitionBase>(keys::partitionManager)->setRestartFlags( RestartFlags::NO_WRITE );

  RegisterGroup( groupKeys.meshBodies );
  RegisterGroup<constitutive::ConstitutiveManager>( groupKeys.constitutiveManager );
  RegisterGroup<CellBlockManager>( keys::cellManager );
}


DomainPartition::~DomainPartition()
{}


void DomainPartition::RegisterDataOnMeshRecursive( Group * const )
{
  Group::RegisterDataOnMeshRecursive( getMeshBodies() );
}


void DomainPartition::InitializationOrder( string_array & order )
{
  SortedArray<string> usedNames;
  {
    order.push_back(keys::ConstitutiveManager);
    usedNames.insert(keys::ConstitutiveManager);
  }

  {
    order.push_back(groupKeysStruct::meshBodiesString);
    usedNames.insert(groupKeysStruct::meshBodiesString);
  }


  for( auto const & subGroup : this->GetSubGroups() )
  {
    if( usedNames.count(subGroup.first) == 0 )
    {
      order.push_back(subGroup.first);
    }
  }
}


void DomainPartition::GenerateSets()
{
  GEOSX_MARK_FUNCTION;

  MeshLevel * const mesh = this->getMeshBody(0)->getMeshLevel(0);
  Group const * const nodeManager = mesh->getNodeManager();

  dataRepository::Group const * const
  nodeSets = nodeManager->GetGroup(ObjectManagerBase::groupKeyStruct::setsString);

  map< string, array1d<bool> > nodeInSet; // map to contain indicator of whether a node is in a set.
  string_array setNames; // just a holder for the names of the sets

  // loop over all wrappers and fill the nodeIndSet arrays for each set
  for( auto & wrapper : nodeSets->wrappers() )
  {
    string name = wrapper.second->getName();
    nodeInSet[name].resize( nodeManager->size() );
    nodeInSet[name] = false;
    Wrapper<SortedArray<localIndex>> const * const setPtr = nodeSets->getWrapper<SortedArray<localIndex>>(name);
    if( setPtr!=nullptr )
    {
      setNames.push_back(name);
      SortedArray<localIndex> const & set = setPtr->reference();
      for( localIndex const a : set )
      {
        nodeInSet[name][a] = true;
      }
    }
  }


  ElementRegionManager * const elementRegionManager = mesh->getElemManager();
  elementRegionManager->forElementSubRegionsComplete( [&]( localIndex const GEOSX_UNUSED_PARAM( er ),
                                                           localIndex const GEOSX_UNUSED_PARAM( esr ),
                                                           ElementRegionBase const * const GEOSX_UNUSED_PARAM( region ),
                                                           auto * const subRegion )
  {
    dataRepository::Group * elementSets = subRegion->sets();

    auto const & elemToNodeMap = subRegion->nodeList();

    for( std::string const & setName : setNames )
    {
      arrayView1d<bool const> const & nodeInCurSet = nodeInSet[setName];

      SortedArray<localIndex> & targetSet = elementSets->registerWrapper< SortedArray<localIndex> >(setName)->reference();
      for( localIndex k = 0 ; k < subRegion->size() ; ++k )
      {
        localIndex const numNodes = subRegion->numNodesPerElement( k );

        localIndex elementInSet = true;
        for ( localIndex i = 0; i < numNodes; ++i )
        {
          if ( !nodeInCurSet( elemToNodeMap[ k ][ i ] ) )
          {
            elementInSet = false;
            break;
          }
        }

        if ( elementInSet )
        {
          targetSet.insert(k);
        }
      }
    }
  });
}


void DomainPartition::SetupCommunications( bool use_nonblocking )
{
  GEOSX_MARK_FUNCTION;
  array1d<NeighborCommunicator> & allNeighbors = this->getReference< array1d<NeighborCommunicator> >( viewKeys.neighbors );

#if defined(GEOSX_USE_MPI)
  if( m_metisNeighborList.empty() )
  {
    PartitionBase   & partition1 = getReference<PartitionBase>(keys::partitionManager);
    SpatialPartition & partition = dynamic_cast<SpatialPartition &>(partition1);

    //get communicator, rank, and coordinates
    MPI_Comm cartcomm;
    {
      int reorder = 0;
      MPI_Cart_create(MPI_COMM_GEOSX, 3, partition.m_Partitions.data(), partition.m_Periodic.data(), reorder, &cartcomm);
      GEOSX_ERROR_IF( cartcomm == MPI_COMM_NULL, "Fail to run MPI_Cart_create and establish communications");
    }
    int const rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
    int nsdof = 3;

    MPI_Cart_coords(cartcomm, rank, nsdof, partition.m_coords.data());

    int ncoords[3];
    AddNeighbors(0, cartcomm, ncoords);

    MPI_Comm_free(&cartcomm);
  }
  else
  {
    for( auto & neighborRank : m_metisNeighborList )
    {
      NeighborCommunicator neighbor;
      neighbor.SetNeighborRank( neighborRank );
      allNeighbors.push_back( std::move(neighbor) );
    }
  }
#endif

  Group * const meshBodies = getMeshBodies();
  MeshBody * const meshBody = meshBodies->GetGroup<MeshBody>(0);
  MeshLevel * const meshLevel = meshBody->GetGroup<MeshLevel>(0);

  for( auto const & neighbor : allNeighbors )
  {
    neighbor.AddNeighborGroupToMesh(meshLevel);
  }

  NodeManager * const nodeManager = meshLevel->getNodeManager();
  FaceManager * const faceManager = meshLevel->getFaceManager();
  EdgeManager * const edgeManager = meshLevel->getEdgeManager();

  nodeManager->SetMaxGlobalIndex();

  CommunicationTools::AssignGlobalIndices( *faceManager, *nodeManager, allNeighbors );

  CommunicationTools::AssignGlobalIndices( *edgeManager, *nodeManager, allNeighbors );

  CommunicationTools::FindMatchedPartitionBoundaryObjects( faceManager,
                                                           allNeighbors );

  CommunicationTools::FindMatchedPartitionBoundaryObjects( nodeManager,
                                                           allNeighbors );

  CommunicationTools::FindGhosts( meshLevel, allNeighbors, use_nonblocking );

  faceManager->SortAllFaceNodes( nodeManager, meshLevel->getElemManager() );
  faceManager->computeGeometry( nodeManager );
}

void DomainPartition::AddNeighbors(const unsigned int idim,
                                   MPI_Comm& cartcomm,
                                   int* ncoords)
{
  PartitionBase   & partition1 = getReference<PartitionBase>(keys::partitionManager);
  SpatialPartition & partition = dynamic_cast<SpatialPartition &>(partition1);
  array1d<NeighborCommunicator> & allNeighbors = this->getReference< array1d<NeighborCommunicator> >( viewKeys.neighbors );

  if (idim == nsdof)
  {
    bool me = true;
    for ( int i = 0 ; i < nsdof ; i++)
    {
      if (ncoords[i] != partition.m_coords(i))
      {
        me = false;
        break;
      }
    }
    if (!me)
    {
      allNeighbors.push_back(NeighborCommunicator());
      int neighborRank = MpiWrapper::Cart_rank(cartcomm, ncoords);
      allNeighbors.back().SetNeighborRank(neighborRank);
    }
  }
  else
  {
    const int dim = partition.m_Partitions( integer_conversion<localIndex>(idim));
    const bool periodic = partition.m_Periodic(integer_conversion<localIndex>(idim));
    for (int i = -1 ; i < 2 ; i++)
    {
      ncoords[idim] = partition.m_coords(integer_conversion<localIndex>(idim)) + i;
      bool ok = true;
      if (periodic)
      {
        if (ncoords[idim] < 0)
          ncoords[idim] = dim - 1;
        else if (ncoords[idim] >= dim)
          ncoords[idim] = 0;
      }
      else
      {
        ok = ncoords[idim] >= 0 && ncoords[idim] < dim;
      }
      if (ok)
      {
        AddNeighbors(idim + 1, cartcomm, ncoords);
      }
    }
  }
}


void DomainPartition::ReadSilo( const SiloFile& siloFile,
                                const int cycleNum,
                                const realT problemTime,
                                const bool isRestart )
{

  ReadFiniteElementMesh( siloFile, cycleNum, problemTime, isRestart );

//  ReadCommonPlanes( siloFile, cycleNum, problemTime, isRestart );
//  ReadCartesianGrid( siloFile, cycleNum, problemTime, isRestart );
//  m_wellboreManager.ReadSilo( siloFile, "WellboreFields", "wellbore_mesh",
//                              DB_NODECENT, cycleNum, problemTime, isRestart );

}


void DomainPartition::ReadFiniteElementMesh( const SiloFile& GEOSX_UNUSED_PARAM( siloFile ),
                                             const int GEOSX_UNUSED_PARAM( cycleNum ),
                                             const realT GEOSX_UNUSED_PARAM( problemTime ),
                                             const bool GEOSX_UNUSED_PARAM( isRestart ) )
{


//  int err = m_feNodeManager->ReadSilo( siloFile, "NodalFields", "volume_mesh",
//                                      DB_NODECENT, cycleNum, problemTime,
// isRestart );
////  err = m_feNodeManager->ReadSilo( siloFile, "NodalFieldsB", "face_mesh",
////                                      DB_NODECENT, cycleNum, problemTime,
// isRestart );
//  if(err)
//    return;
//
//  m_feElementManager->ReadSilo( siloFile, "volume_mesh",
//                               cycleNum, problemTime, isRestart );
//
//  m_feNodeManager->ConstructNodeToElementMap( m_feElementManager );

}





} /* namespace geosx */
