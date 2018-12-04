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
 * @file DomainPartition.cpp
 */

#include "DomainPartition.hpp"

#include "MPI_Communications/SpatialPartition.hpp"
#include "constitutive/ConstitutiveManager.hpp"

#include "fileIO/silo/SiloFile.hpp"

#include "common/TimingMacros.hpp"

#include "common/DataTypes.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "managers/ObjectManagerBase.hpp"
namespace geosx
{
using namespace dataRepository;

DomainPartition::DomainPartition( std::string const & name,
                                  ManagedGroup * const parent ):
  ManagedGroup( name, parent )
{
  this->RegisterViewWrapper< array1d<NeighborCommunicator> >(viewKeys.neighbors)->setRestartFlags( RestartFlags::NO_WRITE );
  this->RegisterViewWrapper<SpatialPartition,PartitionBase>(keys::partitionManager)->setRestartFlags( RestartFlags::NO_WRITE );

  RegisterGroup( groupKeys.meshBodies );
  RegisterGroup<constitutive::ConstitutiveManager>( groupKeys.constitutiveManager );
  RegisterGroup<CellBlockManager>( keys::cellManager );
}


DomainPartition::~DomainPartition()
{}


void DomainPartition::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName("Domain");
  docNode->setSchemaType("UniqueNode");
}


void DomainPartition::InitializationOrder( string_array & order )
{
  set<string> usedNames;
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


void DomainPartition::SetMaps(  )
{
  // ManagedGroup * nodeManager = this->GetGroup(keys::FEM_Nodes);
  // ElementRegionManager * elementRegionManager =
  // this->GetGroup<ElementRegionManager>(keys::FEM_Elements);

  // {
  //  integer_array & elementRegionMap =
  // nodeManager->getReference(keys::elementRegionMap);
  //  integer_array & elementSubRegionMap =
  // nodeManager->getReference(keys::elementSubRegionMap);
  //  integer_array & elementMap = nodeManager->getReference(keys::elementMap);

  //  ManagedGroup * elementRegions =
  // this->GetGroup(dataRepository::keys::elementRegions);

  //  integer elementRegionIndex = 0;
  //  elementRegionManager->forElementRegions( [&](ElementRegion&
  // elementRegion)-> void
  //  {
  //    elementRegion.forCellBlocks( [&]( CellBlockSubRegion & subRegion )->void
  //    {

  //    });
  //    ++elementRegionIndex;
  //  });
  // }
}

void DomainPartition::GenerateSets(  )
{
  MeshLevel * const mesh = this->getMeshBody(0)->getMeshLevel(0);
  ManagedGroup * nodeManager = mesh->getNodeManager();

  dataRepository::ManagedGroup const * nodeSets = nodeManager->GetGroup(dataRepository::keys::sets);

  std::map< string, integer_array > nodeInSet;
  string_array setNames;

  for( auto & viewWrapper : nodeSets->wrappers() )
  {
    string name = viewWrapper.second->getName();
    nodeInSet[name].resize( nodeManager->size() );
    nodeInSet[name] = 0;
    ViewWrapper<set<localIndex>> const * const setPtr = nodeSets->getWrapper<set<localIndex>>(name);
    if( setPtr!=nullptr )
    {
      setNames.push_back(name);
      set<localIndex> const & set = setPtr->reference();
      for( auto const a : set )
      {
        nodeInSet[name][a] = 1;
      }
    }
  }


  ElementRegionManager * elementRegionManager = mesh->getElemManager();

  for( auto & subGroup : elementRegionManager->GetGroup( dataRepository::keys::elementRegions )->GetSubGroups() )
//  elementRegionManager->forElementRegions( [&]( ElementRegion * elementRegion
// )
  {
    ElementRegion * elementRegion = subGroup.second->group_cast<ElementRegion *>();
//    elementRegion->forCellBlocks( [&]( CellBlockSubRegion * subRegion )->void
    for( auto & subRegionIter : elementRegion->GetGroup(dataRepository::keys::cellBlockSubRegions)->GetSubGroups() )
    {
      CellBlockSubRegion * subRegion = subRegionIter.second->group_cast<CellBlockSubRegion *>();
      array2d<localIndex> const & elemsToNodes = subRegion->getWrapper<FixedOneToManyRelation>(subRegion->viewKeys().nodeList)->reference();// getData<array2d<localIndex>>(keys::nodeList);
      dataRepository::ManagedGroup * elementSets = subRegion->GetGroup(dataRepository::keys::sets);
      std::map< string, integer_array > numNodesInSet;

      for( auto & setName : setNames )
      {

        set<localIndex> & targetSet = elementSets->RegisterViewWrapper< set<localIndex> >(setName)->reference();
        for( localIndex k = 0 ; k < subRegion->size() ; ++k )
        {
          integer count = 0;
          for( localIndex a = 0 ; a<elemsToNodes.size(1) ; ++a )
          {
            if( nodeInSet[setName][elemsToNodes[k][a]] == 1 )
            {
              ++count;
            }
          }
          if( count == elemsToNodes.size(1) )
          {
            targetSet.insert(k);
          }
        }
      }
    }
  }
}


void DomainPartition::SetupCommunications()
{
  GEOSX_MARK_FUNCTION;
  PartitionBase   & partition1 = getReference<PartitionBase>(keys::partitionManager);
  SpatialPartition & partition = dynamic_cast<SpatialPartition &>(partition1);
  array1d<NeighborCommunicator> & allNeighbors = this->getReference< array1d<NeighborCommunicator> >( viewKeys.neighbors );

  //get communicator, rank, and coordinates
  MPI_Comm cartcomm;
  {
    int reorder = 0;
    MPI_Cart_create(MPI_COMM_GEOSX, 3, partition.m_Partitions.data(), partition.m_Periodic.data(), reorder, &cartcomm);
  }
  int rank = -1;
  int nsdof = 3;
  MPI_Comm_rank( MPI_COMM_GEOSX, &rank );

  MPI_Comm_rank(cartcomm, &rank);
  MPI_Cart_coords(cartcomm, rank, nsdof, partition.m_coords.data());

  int ncoords[3];
  AddNeighbors(0, cartcomm, ncoords);

  MPI_Comm_free(&cartcomm);


  ManagedGroup * const meshBodies = getMeshBodies();
  MeshBody * const meshBody = meshBodies->GetGroup<MeshBody>(0);
  MeshLevel * const meshLevel = meshBody->GetGroup<MeshLevel>(0);

  for( auto const & neighbor : allNeighbors )
  {
    neighbor.AddNeighborGroupToMesh(meshLevel);
  }



  NodeManager * nodeManager = meshLevel->getNodeManager();
  FaceManager * const faceManager = meshLevel->getFaceManager();

  EdgeManager * const edgeManager = meshLevel->getEdgeManager();

  nodeManager->SetMaxGlobalIndex();

  CommunicationTools::AssignGlobalIndices( *faceManager, *nodeManager, allNeighbors );

  CommunicationTools::AssignGlobalIndices( *edgeManager, *nodeManager, allNeighbors );

  CommunicationTools::FindMatchedPartitionBoundaryObjects( faceManager,
                                                           allNeighbors );

  CommunicationTools::FindMatchedPartitionBoundaryObjects( nodeManager,
                                                           allNeighbors );

  CommunicationTools::FindGhosts( meshLevel, allNeighbors );

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
      int neighborRank;
      MPI_Cart_rank(cartcomm, ncoords, &neighborRank);
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


void DomainPartition::ReadFiniteElementMesh( const SiloFile& siloFile,
                                             const int cycleNum,
                                             const realT problemTime,
                                             const bool isRestart )
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
