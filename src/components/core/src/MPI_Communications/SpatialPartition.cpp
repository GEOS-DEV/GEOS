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
 * @file SpatialPartition.cpp
 * @author settgast1
 * @date Mar 16, 2011
 */

#include "SpatialPartition.hpp"

#include "codingUtilities/Utilities.hpp"

//#include "Common/intrinsic_typedefs.h"
#include <cmath>

#ifdef USE_ATK
#include "slic/slic.hpp"
#endif



namespace
{

// Modulo
// returns a positive value regardless of the sign of numerator
realT Mod(realT num, realT denom)
{
  if( fabs(denom)<fabs(num)*1.0e-14 )
  {
    return num;
  }

  return num - denom * std::floor(num/denom);
}


// MapValueToRange
// returns a periodic value in the range [min, max)
realT MapValueToRange(realT value, realT min, realT max)
{
  return Mod(value-min, max-min)+min;
}



}

namespace geosx
{
using namespace dataRepository;

SpatialPartition::SpatialPartition():
  PartitionBase(),
  m_Partitions(),
  m_Periodic(nsdof),
  m_coords(nsdof),
  m_min(0.0),
  m_max(0.0),
  m_blockSize(1),
  m_gridSize(0.0),
  m_gridMin(0.0),
  m_gridMax(0.0)
{
  m_size = 0;
  m_rank = 0;
  setPartitions(1,1,1);
}

SpatialPartition::~SpatialPartition()
{}

//void SpatialPartition::ReadXML( xmlWrapper::xmlNode const & targetNode )
//{
//  int xpar  = targetNode.attribute("xpar").as_int(1);
//
//}


void SpatialPartition::InitializePostSubGroups( ManagedGroup * const  )
{
  //get size of problem and decomposition
  MPI_Comm_size(MPI_COMM_WORLD, &m_size);

  //check to make sure our dimensions agree
  {
    int check = 1;
    for (unsigned int i = 0 ; i < nsdof ; i++)
    {
      check *= this->m_Partitions(i);
    }
    assert(check == m_size);
  }

  //get communicator, rank, and coordinates
  MPI_Comm cartcomm;
  {
    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, nsdof, m_Partitions.data(), m_Periodic.data(), reorder, &cartcomm);
  }
  MPI_Comm_rank(cartcomm, &m_rank);
  MPI_Cart_coords(cartcomm, m_rank, nsdof, m_coords.data());


  m_color = GetColor();

  //add neighbors
  {
    int ncoords[nsdof];
    m_neighbors.clear();
    AddNeighbors(0, cartcomm, ncoords);
  }

  MPI_Comm_free(&cartcomm);

  //initialize cached requests and status
  m_mpiRequest.resize( 2 * m_neighbors.size() );
  m_mpiStatus.resize( 2 * m_neighbors.size() );
}

void SpatialPartition::InitializeMetis()
{
  //get size of problem and decomposition
  MPI_Comm_size(MPI_COMM_WORLD, &m_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
  //check to make sure our dimensions agree
  {
    assert(m_sizeMetis == m_size);
  }
  //initialize cached requests and status
  m_mpiRequest.resize( 100 );
  m_mpiStatus.resize( 100 );

}

int SpatialPartition::GetColor()
{
  int color = 0;

  if( isOdd(m_coords[0]) )
  {
    color += 1;
  }

  if( isOdd(m_coords[1]) )
  {
    color += 2;
  }

  if( isOdd(m_coords[2]) )
  {
    color += 4;
  }

  m_numColors = 8;

  return color;
}


/*
 * @author walsh24
 *
 * corrects the global ids for the nodes on periodic boundaries, and marks
 * periodic boundaries as domain boundaries
 *
 */
//void SpatialPartition::SetPeriodicDomainBoundaryObjects(DomainPartition&
// domain ){
//
//   gArray1d& nodeGlobalIds = domain.m_feNodeManager.m_localToGlobalMap;
////   gArray1d& faceGlobalIds = domain.m_feFaceManager.m_localToGlobalMap;
////   gArray1d& edgeGlobalIds = domain.m_feEdgeManager.m_localToGlobalMap;
//
//   std::map<globalIndex,localIndex>& nodeGlobalToLocalMap =
// domain.m_feNodeManager.m_globalToLocalMap;
//   array< R1Tensor >& refPositions = *(domain.m_feNodeManager.m_refposition);
//
//   for(unsigned int i =0; i < m_periodicSets.size();++i){
//     int dimension = m_periodicSets[i].m_dimension;
//
//     if( (m_coords[dimension] == 0)  ||
//         (m_coords[dimension] == m_Partitions[dimension]-1) )
//     {
//
//       // Reset global id numbers
//       ///////////////////////////
//
//       array<string>& setnames = m_periodicSets[i].m_setNames;
//       set<localIndex>* theSets[2];
//       theSets[0] = &(domain.m_feNodeManager.m_Sets[setnames[0]]);
//       theSets[1] = &(domain.m_feNodeManager.m_Sets[setnames[1]]);
//
//       PlanarSorter planarSorter(refPositions,dimension);
//
//       if(m_Partitions[dimension] > 1){
//
//         // Multiple partitions
//         //--------------------
//
//         array<int> nbr_coords = m_coords;
//         if(m_coords[dimension] == 0){
//           nbr_coords[dimension] = m_Partitions[dimension]-1;
//         } else {
//           nbr_coords[dimension] = 0;
//         }
//
//         int mySetId = (theSets[0]->size() > 0)? 0 : 1;
//         int nbrSetId = 1-mySetId;
//         if(theSets[nbrSetId]->size() > 0)
//         {
//#ifdef USE_ATK
//           SLIC_ERROR("SpatialPartition::SetPeriodicDomainBoundaryObjects: " +
// setnames[0] + " and " + setnames[1] + " present on same partition\n");
//#endif
//         }
//         set<localIndex>& mySet =  *(theSets[mySetId]);
//
//         // gather local and global ids
//         std::vector<std::pair<localIndex, localIndex>  > myLocalAndGlobalIds;
//
//         for( set<localIndex>::iterator itr=mySet.begin() ; itr!=mySet.end() ; ++itr )
//         {
//           localIndex globalId = nodeGlobalIds[*itr];
//           myLocalAndGlobalIds.push_back(std::pair<localIndex , localIndex>(
// *itr,globalId) );
//         }
//
//         // Sort local/global ids by position in plane
//         std::vector<localIndex>
// mySortedGlobalIds(myLocalAndGlobalIds.size());
//         std::vector<localIndex> nbrSortedGlobalIds;
//
//         std::sort(myLocalAndGlobalIds.begin(),myLocalAndGlobalIds.end(),planarSorter);
//         for(unsigned int ii = 0 ; ii <myLocalAndGlobalIds.size() ; ++ii ){
//           mySortedGlobalIds[ii]  = myLocalAndGlobalIds[ii].second;
//         }
//
//         // Communicate global ids between periodic domains
//         NeighborCommunication& neighbor = m_neighbors[
// neighborCommPtrIndx[nbr_coords] ];
//
//         MPI_Request mpiRequest[2];
//         MPI_Status mpiStatus;
//
//         int receiveSize; int sendSize = mySortedGlobalIds.size();
//
//         neighbor.SendReceive( reinterpret_cast<char*>(&receiveSize),
// sizeof(int),
//                               reinterpret_cast<char*>(&sendSize),
// sizeof(int),
//                               mpiRequest[0], mpiRequest[1] );
//
//         MPI_Waitall( 2, mpiRequest, &mpiStatus );
//
//         nbrSortedGlobalIds.resize(receiveSize);
//         neighbor.SendReceive(
// reinterpret_cast<char*>(&nbrSortedGlobalIds[0]),
// receiveSize*sizeof(localIndex),
//                               reinterpret_cast<char*>(&mySortedGlobalIds[0]),
// sendSize*sizeof(localIndex),
//                               mpiRequest[0], mpiRequest[1] );
//
//         MPI_Waitall( 2, mpiRequest, &mpiStatus );
//
//
//
//         // should have same number of nodes in both sets
//         if(nbrSortedGlobalIds.size() !=  mySortedGlobalIds.size() )
//         {
//#ifdef USE_ATK
//           SLIC_ERROR("SpatialPartition::SetPeriodicDomainBoundaryObjects:
// Size of " + setnames[mySetId] + " does not match size of " +
// setnames[nbrSetId] + " on neighboring partition\n");
//#endif
//         }
//
//
//         // assign new global ids
//         for(unsigned int ii = 0 ; ii <myLocalAndGlobalIds.size() ; ++ii ){
//           localIndex& nd =  myLocalAndGlobalIds[ii].first;
//           nodeGlobalIds[nd] =
// std::min(mySortedGlobalIds[ii],nbrSortedGlobalIds[ii]);
//         }
//
//       } else {
//
//         // Single partition
//         //-----------------
//
//         // Nodes
//         {
//           std::vector< std::vector<std::pair<localIndex, localIndex>  >  >
// setLocalAndGlobalIds(2);
//           for(int a =0; a<2; ++a){
//             // Gather local/global ids
//             for( set<localIndex>::iterator itr=theSets[a]->begin() ;
// itr!=theSets[a]->end() ; ++itr )
//             {
//               localIndex globalId = nodeGlobalIds[*itr];
//               setLocalAndGlobalIds[a].push_back(std::pair<localIndex ,
// localIndex>( *itr,globalId) );
//             }
//             // Sort local/global ids by position in plane
//             std::sort(setLocalAndGlobalIds[a].begin(),setLocalAndGlobalIds[a].end(),planarSorter);
//           }
//
//           // should have same number of nodes in both sets
//           if(setLocalAndGlobalIds[0].size() !=
//  setLocalAndGlobalIds[1].size() )
//           {
//#ifdef USE_ATK
//             SLIC_ERROR("SpatialPartition::SetPeriodicDomainBoundaryObjects:
// Size of " + setnames[0] + " does not match size of " + setnames[1] + " on
// process " +toString(m_rank) +  "\n");
//#endif
//           }
//
//           // assign new global ids and make global to local map point to
// nodes on min boundary
//           for(unsigned int ii = 0 ; ii <setLocalAndGlobalIds[0].size() ; ++ii
// ){
//             localIndex& nd0 =  setLocalAndGlobalIds[0][ii].first;
//             localIndex& nd1 =  setLocalAndGlobalIds[1][ii].first;
//
//             // this could be done once (all nodes in the same set should lie
// on the one boundary)
//             int minBoundarySetIndx = 0;
//             if(  (*domain.m_feNodeManager.m_refposition)[nd1][dimension] <
// (*domain.m_feNodeManager.m_refposition)[nd0][dimension] ){
//               minBoundarySetIndx = 1;
//             }
//             int maxBoundarySetIndx = 1 - minBoundarySetIndx;
//             localIndex localTarget = (minBoundarySetIndx == 0)? nd0 : nd1;
////             localIndex notThelocalTarget = (minBoundarySetIndx == 0)? nd1 :
// nd0;
//
//             // fix up local to global map
//             localIndex minBoundGlobalId =
// setLocalAndGlobalIds[minBoundarySetIndx][ii].second;
//             localIndex maxBoundGlobalId =
// setLocalAndGlobalIds[maxBoundarySetIndx][ii].second;
//
//             nodeGlobalIds[nd0] = minBoundGlobalId;
//             nodeGlobalIds[nd1] = minBoundGlobalId;
//
//             // fix up global to local map
//             nodeGlobalToLocalMap[minBoundGlobalId] = localTarget;
//
//             // not used? in any case make old Global id point to same local
// target
//             nodeGlobalToLocalMap[maxBoundGlobalId] = localTarget;
//           }
//         }
//
//       }
//
//       // set boundary ids
//       /////////////////////
//       std::cout << " Setting new domain boundary:  " << std::endl;
//       domain.m_feNodeManager.SetSubsetValue<int>("isDomainBoundary",setnames[0],1);
//       domain.m_feEdgeManager.SetSubsetValue<int>("isDomainBoundary",setnames[0],1);
//       domain.m_feFaceManager.SetSubsetValue<int>("isDomainBoundary",setnames[0],1);
//
//       domain.m_feNodeManager.SetSubsetValue<int>("isDomainBoundary",setnames[1],1);
//       domain.m_feEdgeManager.SetSubsetValue<int>("isDomainBoundary",setnames[1],1);
//       domain.m_feFaceManager.SetSubsetValue<int>("isDomainBoundary",setnames[1],1);
//
//     }
//   }
//
//}



/*
 * @author walsh24
 *
 * resets the global to local node map, so global ids map to local nodes that
 * lie on minimum boundary
 *
 */
//void SpatialPartition::ResetSinglePartitionGlobalToLocalMap(PhysicalDomainT&
// domain ){
//
//  // loop over matched sets
//  for(unsigned int ps =0; ps < m_periodicSets.size();++ps){
//
//    int dimension = m_periodicSets[ps].m_dimension;
//
//    // single partition
//    if(m_Partitions[dimension]==1){
//
//      array<string>& setnames = m_periodicSets[ps].m_setNames;
//      set<localIndex>* theNodeSets[2];
//      theNodeSets[0] = &(domain.m_feNodeManager.m_Sets[setnames[0]]);
//      theNodeSets[1] = &(domain.m_feNodeManager.m_Sets[setnames[1]]);
//      set<localIndex>* theFaceSets[2];
//      theFaceSets[0] = &(domain.m_feFaceManager.m_Sets[setnames[0]]);
//      theFaceSets[1] = &(domain.m_feFaceManager.m_Sets[setnames[1]]);
//      set<localIndex>* theEdgeSets[2];
//      theEdgeSets[0] = &(domain.m_feEdgeManager.m_Sets[setnames[0]]);
//      theEdgeSets[1] = &(domain.m_feEdgeManager.m_Sets[setnames[1]]);
//
//      gArray1d& nodeLocalToGlobalMap =
// domain.m_feNodeManager.m_localToGlobalMap;
//      gArray1d& faceLocalToGlobalMap =
// domain.m_feFaceManager.m_localToGlobalMap;
//      gArray1d& edgeLocalToGlobalMap =
// domain.m_feEdgeManager.m_localToGlobalMap;
//
//      std::map<globalIndex,localIndex>& nodeGlobalToLocalMap =
// domain.m_feNodeManager.m_globalToLocalMap;
//      std::map<globalIndex,localIndex>& faceGlobalToLocalMap =
// domain.m_feFaceManager.m_globalToLocalMap;
//      std::map<globalIndex,localIndex>& edgeGlobalToLocalMap =
// domain.m_feEdgeManager.m_globalToLocalMap;
//
//      // identify sets on min and max boundaries
//      int minBoundarySetIndx = 0;
//      {
//        localIndex nd0 = *(theNodeSets[0]->begin());
//        localIndex nd1 = *(theNodeSets[1]->begin());
//
//        if(  (*domain.m_feNodeManager.m_refposition)[nd1][dimension] <
// (*domain.m_feNodeManager.m_refposition)[nd0][dimension] ){
//          minBoundarySetIndx = 1;
//        }
//      }
////      int maxBoundarySetIndx = 1 - minBoundarySetIndx;
//
//      // point global to local map to min boundary nodes
//      {
//        set<localIndex>::iterator itr = theNodeSets[minBoundarySetIndx]->begin();
//        set<localIndex>::iterator iend = theNodeSets[minBoundarySetIndx]->end();
//        for( ; itr != iend ; ++itr ){
//          globalIndex g = nodeLocalToGlobalMap[*itr];
//          nodeGlobalToLocalMap[g] = *itr;
//        }
//      }
//
//      // set local to global map for faces and edges
//      //////////////////////////////////////////////
//      // Faces
//      // make local to global maps for boundary faces point to same global
// index
//      // requires nodeLocalToGlobalMap to be set
//      {
//        std::map<set<localIndex>,localIndex> nodeToFaceGlobalMap;
//
//        theFaceSets[0] = &(domain.m_feFaceManager.m_Sets[setnames[0]]);
//        theFaceSets[1] = &(domain.m_feFaceManager.m_Sets[setnames[1]]);
//
//
//        for( set<localIndex>::iterator itr = theFaceSets[0]->begin(); itr !=
// theFaceSets[0]->end(); ++itr){
//          localIndex faceLocal = *itr;
//          lArray1d& faceNodes =
// domain.m_feFaceManager.m_toNodesRelation[faceLocal];
//
//          set<localIndex> gnodes;
//          for(unsigned ii = 0; ii < faceNodes.size(); ++ii){
//            gnodes.insert(nodeLocalToGlobalMap[faceNodes[ii]]);
//          }
//          localIndex gf = faceLocalToGlobalMap[faceLocal];
//          nodeToFaceGlobalMap[gnodes] = gf;
//          faceGlobalToLocalMap[gf] = faceLocal;
//        }
//
//        for( set<localIndex>::iterator itr = theFaceSets[1]->begin(); itr !=
// theFaceSets[1]->end(); ++itr){
//          localIndex faceLocal = *itr;
//          lArray1d& faceNodes =
// domain.m_feFaceManager.m_toNodesRelation[faceLocal];
//
//          set<localIndex> gnodes;
//          for(unsigned ii = 0; ii < faceNodes.size(); ++ii){
//            gnodes.insert(nodeLocalToGlobalMap[faceNodes[ii]]);
//          }
//          faceLocalToGlobalMap[faceLocal] = nodeToFaceGlobalMap[gnodes];
//        }
//      }
//
//      // Edges
//      // make local to global maps for boundary edges point to same global
// index
//      // requires nodeLocalToGlobalMap to be set
//      {
//        std::map<set<localIndex>,localIndex> nodeToEdgeGlobalMap;
//
//        theEdgeSets[0] = &(domain.m_feEdgeManager.m_Sets[setnames[0]]);
//        theEdgeSets[1] = &(domain.m_feEdgeManager.m_Sets[setnames[1]]);
//
//        for( set<localIndex>::iterator itr = theEdgeSets[0]->begin(); itr !=
// theEdgeSets[0]->end(); ++itr){
//          localIndex edgeLocal = *itr;
//          lArray2d& edgeNodes = domain.m_feEdgeManager.m_toNodesRelation;
//
//          set<localIndex> gnodes;
//          for(unsigned ii = 0; ii < 2; ++ii){
//            gnodes.insert(nodeLocalToGlobalMap[ edgeNodes(edgeLocal,ii) ]);
//          }
//
//          localIndex ge = edgeLocalToGlobalMap[edgeLocal];
//          nodeToEdgeGlobalMap[gnodes] = ge;
//          edgeGlobalToLocalMap[ge] = edgeLocal;
//        }
//
//
//        for( set<localIndex>::iterator itr = theEdgeSets[1]->begin(); itr !=
// theEdgeSets[1]->end(); ++itr){
//          localIndex edgeLocal = *itr;
//          lArray2d& edgeNodes = domain.m_feEdgeManager.m_toNodesRelation;
//
//          set<localIndex> gnodes;
//          for(unsigned ii = 0; ii < 2; ++ii){
//            gnodes.insert(nodeLocalToGlobalMap[ edgeNodes(edgeLocal,ii) ]);
//          }
//          edgeLocalToGlobalMap[edgeLocal] = nodeToEdgeGlobalMap[gnodes];
//        }
//      }
//
//      // set global to local maps for faces and edges
//      //////////////////////////////////////////////
//
//      // point global to local map to min boundary faces
//      {
//        set<localIndex>::iterator itr = theFaceSets[minBoundarySetIndx]->begin();
//        set<localIndex>::iterator iend = theFaceSets[minBoundarySetIndx]->end();
//        for( ; itr != iend ; ++itr ){
//          globalIndex g = faceLocalToGlobalMap[*itr];
//          faceGlobalToLocalMap[g] = *itr;
//        }
//      }
//
//      // point global to local map to min boundary edges
//      {
//        set<localIndex>::iterator itr = theEdgeSets[minBoundarySetIndx]->begin();
//        set<localIndex>::iterator iend = theEdgeSets[minBoundarySetIndx]->end();
//        for( ; itr != iend ; ++itr ){
//          globalIndex g = edgeLocalToGlobalMap[*itr];
//          edgeGlobalToLocalMap[g] = *itr;
//        }
//      }
//
//    }
//  }
//
//}

/**
 * @author walsh24
 *
 * Creates ghosts objects when required on the same partition
 *
 **/
//void SpatialPartition::CreateSinglePartitionGhostObjects(PhysicalDomainT&
// domain,
//                                                         const bool
// contactActive,
//                                                         const int
// elementGhostingDepth )
//{
//
//  array<PhysicalDomainT::ObjectDataStructureKeys> objectNames;
//  NeighborCommunication::SyncNames(objectNames);
//
//
//  // loop over matched sets
//  for(unsigned int ps =0; ps < m_periodicSets.size();++ps){
//
//    int dimension = m_periodicSets[ps].m_dimension;
//
//    // single partition
//    if(m_Partitions[dimension]==1){
//
//      if (m_rank  == 0)  std::cout << "Setting ghost objects for single
// periodic partitions." << std::endl;
//
//      m_hasLocalGhosts = true;
//
//      realT myMax = m_max[dimension];
//      realT myMin = m_min[dimension];
//      realT myMid = 0.5*(myMax+myMin);
//      realT gridLength = m_gridSize[dimension];
//
//      array<string>& setnames = m_periodicSets[ps].m_setNames;
//      set<localIndex>* theNodeSets[2];
//      theNodeSets[0] = &(domain.m_feNodeManager.m_Sets[setnames[0]]);
//      theNodeSets[1] = &(domain.m_feNodeManager.m_Sets[setnames[1]]);
//      set<localIndex>* theFaceSets[2];
//      theFaceSets[0] = &(domain.m_feFaceManager.m_Sets[setnames[0]]);
//      theFaceSets[1] = &(domain.m_feFaceManager.m_Sets[setnames[1]]);
//      set<localIndex>* theEdgeSets[2];
//      theEdgeSets[0] = &(domain.m_feEdgeManager.m_Sets[setnames[0]]);
//      theEdgeSets[1] = &(domain.m_feEdgeManager.m_Sets[setnames[1]]);
//
//      gArray1d& nodeLocalToGlobalMap =
// domain.m_feNodeManager.m_localToGlobalMap;
//      gArray1d& faceLocalToGlobalMap =
// domain.m_feFaceManager.m_localToGlobalMap;
//      gArray1d& edgeLocalToGlobalMap =
// domain.m_feEdgeManager.m_localToGlobalMap;
//
//      std::map<globalIndex,localIndex>& nodeGlobalToLocalMap =
// domain.m_feNodeManager.m_globalToLocalMap;
//      std::map<globalIndex,localIndex>& faceGlobalToLocalMap =
// domain.m_feFaceManager.m_globalToLocalMap;
//      std::map<globalIndex,localIndex>& edgeGlobalToLocalMap =
// domain.m_feEdgeManager.m_globalToLocalMap;
//
//
//      // pack ghost objects
//      /////////////////////
//
//      // must perform two self communications - one to lower boundary and one
// to upper
//      // need to pack both before unpacking - otherwise newly created boundary
// will be sent also
//      // also need to make sure that globaltolocal index points to the correct
// boundary when unpacking
//      // otherwise will unpack but link to lower boundary.
//      NeighborCommunication selfCommunication[2];
//
//      for(int a =0; a < 2; ++a){
//        selfCommunication[a].SetDomain(domain);
//        selfCommunication[a].Initialize(m_rank,m_rank,-1);
//        lArray1d& setNodes =
// selfCommunication[a].tempNeighborData.matchedIndices[PhysicalDomainT::FiniteElementNodeManager];
//        // gather boundary nodes and point global to local map to the same set
//        for( set<localIndex>::iterator itr=theNodeSets[a]->begin() ;
// itr!=theNodeSets[a]->end() ; ++itr )
//        {
//          globalIndex g = nodeLocalToGlobalMap[*itr];
//          nodeGlobalToLocalMap[g] = *itr;
//
//          setNodes.push_back(*itr );
//        }
//        for( set<localIndex>::iterator itr=theFaceSets[a]->begin() ;
// itr!=theFaceSets[a]->end() ; ++itr )
//        {
//          globalIndex g = faceLocalToGlobalMap[*itr];
//          faceGlobalToLocalMap[g] = *itr;
//        }
//        for( set<localIndex>::iterator itr=theEdgeSets[a]->begin() ;
// itr!=theEdgeSets[a]->end() ; ++itr )
//        {
//          globalIndex g = edgeLocalToGlobalMap[*itr];
//          edgeGlobalToLocalMap[g] = *itr;
//        }
//
//        // pack ghosts
//        //TODO Stuart, does this need to also communicate the FE?
//        selfCommunication[a].FindPackGhostsDiscreteElement( contactActive,
// elementGhostingDepth);
//#ifdef SRC_EXTERNAL
//        selfCommunication[a].FindPackGhostsFaultElement( elementGhostingDepth
// );
//#endif
//      }
//
//
//      // unpack ghosts
//      ////////////////
//
//      // need to check that global index of set points to receiving set
//      for(int a =0; a < 2; ++a){
//
//        // change global to local map to point to receiving set
//        int nota = 1-a;
//        for( set<localIndex>::iterator itr=theNodeSets[nota]->begin() ;
// itr!=theNodeSets[nota]->end() ; ++itr )
//        {
//          globalIndex g = nodeLocalToGlobalMap[*itr];
//          nodeGlobalToLocalMap[g] = *itr;
//        }
//        for( set<localIndex>::iterator itr=theFaceSets[nota]->begin() ;
// itr!=theFaceSets[nota]->end() ; ++itr )
//        {
//          globalIndex g = faceLocalToGlobalMap[*itr];
//          faceGlobalToLocalMap[g] = *itr;
//        }
//        for( set<localIndex>::iterator itr=theEdgeSets[nota]->begin() ;
// itr!=theEdgeSets[nota]->end() ; ++itr )
//        {
//          globalIndex g = edgeLocalToGlobalMap[*itr];
//          edgeGlobalToLocalMap[g] = *itr;
//        }
//
//        localIndex firstNodeIndex = m_domain->m_feNodeManager.DataLengths();
//        localIndex firstFaceIndex = m_domain->m_feFaceManager.DataLengths();
//        localIndex firstEdgeIndex = m_domain->m_feEdgeManager.DataLengths();
//
//        const lArray1d& nodeList =
// selfCommunication[a].GetSendLocalIndices(PhysicalDomainT::FiniteElementNodeManager);
//        const lArray1d& edgeList =
// selfCommunication[a].GetSendLocalIndices(PhysicalDomainT::FiniteElementEdgeManager);
//        const lArray1d& faceList =
// selfCommunication[a].GetSendLocalIndices(PhysicalDomainT::FiniteElementFaceManager);
//
//        // copy element Lists - is updated later
//        std::map<std::string,lArray1d> elementLists =
// selfCommunication[a].GetElementRegionSendLocalIndices();
//
//        // map old local indicies to new local indices
//        ///////////////////////////////////////////////
//        std::map<localIndex,localIndex> newNodeMap;
//        std::map<localIndex,localIndex> newEdgeMap;
//        std::map<localIndex,localIndex> newFaceMap;
//        std::map<localIndex,localIndex> oldNodeMap;
//        std::map<localIndex,localIndex> oldEdgeMap;
//        std::map<localIndex,localIndex> oldFaceMap;
//
//        localIndex numNewNodes = 0;
//        for( lArray1d::size_type ii = 0; ii < nodeList.size();++ii){
//          localIndex oldIndx = nodeList[ii];
//          localIndex newIndx;
//          if(theNodeSets[a]->count(oldIndx) == 0){
//            newIndx = numNewNodes + firstNodeIndex;
//            numNewNodes++;
//          } else {
//            newIndx = nodeGlobalToLocalMap[nodeLocalToGlobalMap[oldIndx]];
//          }
//          newNodeMap[oldIndx] = newIndx;
//          oldNodeMap[newIndx] = oldIndx;
//        }
//
//        localIndex numNewEdges = 0;
//        for( lArray1d::size_type ii = 0; ii < edgeList.size();++ii){
//          localIndex oldIndx = edgeList[ii];
//          localIndex newIndx;
//          if(theEdgeSets[a]->count(oldIndx) == 0){
//            newIndx = numNewEdges + firstEdgeIndex;
//            numNewEdges++;
//          } else {
//            newIndx = edgeGlobalToLocalMap[edgeLocalToGlobalMap[oldIndx]];
//            edgeLocalToGlobalMap[newIndx] = edgeLocalToGlobalMap[oldIndx];
//          }
//          newEdgeMap[oldIndx] = newIndx;
//          oldEdgeMap[newIndx] = oldIndx;
//        }
//
//        localIndex numNewFaces = 0;
//        for( lArray1d::size_type ii = 0; ii < faceList.size();++ii){
//          localIndex oldIndx = faceList[ii];
//          localIndex newIndx;
//          if(theFaceSets[a]->count(oldIndx) == 0){
//            newIndx = numNewFaces + firstFaceIndex;
//            numNewFaces++;
//          } else {
//            newIndx = faceGlobalToLocalMap[faceLocalToGlobalMap[oldIndx]];
//            faceLocalToGlobalMap[newIndx] = faceLocalToGlobalMap[oldIndx];
//          }
//          newFaceMap[oldIndx] = newIndx;
//          oldFaceMap[newIndx] = oldIndx;
//        }
//
//
//        // create new nodes
//        ///////////////////
//        m_domain->m_feNodeManager.resize(firstNodeIndex+numNewNodes);
//        array<R1Tensor>& referencePosition =
// m_domain->m_feNodeManager.GetFieldData<FieldInfo::referencePosition>();
//
//        for(unsigned  kn = firstNodeIndex; kn <
// m_domain->m_feNodeManager.DataLengths();++kn){
//            localIndex oldNode = oldNodeMap[kn];
//            nodeLocalToGlobalMap[kn] = nodeLocalToGlobalMap[oldNode];
//
//            m_domain->m_feNodeManager.CopyObjectWithExcludedSets( oldNode, kn,
// setnames );
//
//            if( referencePosition[kn][dimension] >  myMid){
//              referencePosition[kn][dimension] -= gridLength;
//            } else {
//              referencePosition[kn][dimension] += gridLength;
//            }
//        }
//
//        // create new edges
//        ////////////////////
//        m_domain->m_feEdgeManager.resize(firstEdgeIndex+numNewEdges);
//        for(unsigned  ke = firstEdgeIndex; ke <
// m_domain->m_feEdgeManager.DataLengths();++ke){
//            localIndex oldEdge = oldEdgeMap[ke];
//            edgeLocalToGlobalMap[ke] = edgeLocalToGlobalMap[oldEdge];
//
//            m_domain->m_feEdgeManager.CopyObjectWithExcludedSets( oldEdge, ke,
// setnames );
//
//            lArray2d& edgesToNodes =
// m_domain->m_feEdgeManager.m_toNodesRelation;
//            for(unsigned jj = 0; jj < 2;++jj){
//              edgesToNodes(ke,jj) = newNodeMap[ edgesToNodes(oldEdge,jj) ];
//            }
//        }
//
//        // create new faces
//        ///////////////////
//        m_domain->m_feFaceManager.resize(firstFaceIndex+numNewFaces);
//
//        for(unsigned kf = firstFaceIndex; kf <
// m_domain->m_feFaceManager.m_numFaces;++kf){
//            localIndex oldFace = oldFaceMap[kf];
//
//            faceLocalToGlobalMap[kf] = faceLocalToGlobalMap[oldFace];
//
//            m_domain->m_feFaceManager.CopyObjectWithExcludedSets( oldFace,
// kf,setnames );
//
//            lArray1d& faceToEdges =
// m_domain->m_feFaceManager.m_toEdgesRelation[kf];
//            faceToEdges =
// m_domain->m_feFaceManager.m_toEdgesRelation(oldFace);
//            for(unsigned jj = 0; jj < faceToEdges.size();++jj){
//              faceToEdges[jj] = newEdgeMap[ faceToEdges[jj] ];
//            }
//            lArray1d& faceToNodes =
// m_domain->m_feFaceManager.m_toNodesRelation[kf];
//            faceToNodes =
// m_domain->m_feFaceManager.m_toNodesRelation(oldFace);
//            for(unsigned jj = 0; jj < faceToNodes.size();++jj){
//              faceToNodes[jj] = newNodeMap[ faceToNodes[jj] ];
//            }
//        }
//
//        // create new elements
//        //////////////////////
//        localIndex totalNumElements = 0;
//        for( std::map<std::string,lArray1d>::iterator itr=elementLists.begin()
// ;
//              itr!=elementLists.end() ; ++itr ){
//
//          ElementRegionT& elementRegion =
//  m_domain->m_feElementManager.m_ElementRegions[itr->first];
//          lArray1d& oldElems = itr->second;
//
//          localIndex firstIndex = elementRegion.DataLengths();
//          localIndex newSize = firstIndex + oldElems.size();
//          elementRegion.resize(newSize);
//          //std::cout << "elems sizes "<< firstIndex << " " << newSize << " "
// << elementRegion.m_numElems << std::endl;
//
//          lArray2d& elemToNodeMap = elementRegion.m_toNodesRelation;
//          lArray2d& elemToFaceMap = elementRegion.m_toFacesRelation;
//
//          gArray1d& elemLocalToGlobalMap = elementRegion.m_localToGlobalMap;
//
//
//          for(unsigned jj =0; jj < oldElems.size();++jj){
//            localIndex oldElem = oldElems[jj];
//            localIndex newElem = firstIndex + jj;
//
//            elementRegion.CopyObjectWithExcludedSets( oldElem, newElem ,
// setnames);
//
//            for(unsigned kk =0; kk < elementRegion.m_numNodesPerElem; ++kk){
//              elemToNodeMap(newElem,kk)
//                   = newNodeMap[ elemToNodeMap(oldElem,kk) ];
//            }
//
//            for( int kk =0; kk < elementRegion.m_numFacesPerElement; ++kk){
//              elemToFaceMap(newElem,kk)
//                   = newFaceMap[ elemToFaceMap(oldElem,kk) ];
//            }
//
//            elemLocalToGlobalMap[newElem] = elemLocalToGlobalMap[oldElem];
//
//            oldElems[jj] = newElem; // make element list refer to new elements
//
//          }
//
//          totalNumElements += elementRegion.m_numElems;
//        }
//
//        m_domain->m_feElementManager.m_numElems = totalNumElements;
//
//
//
//        // (6d) update maps in the anti-dependency direction
//        //      must update the "upward" pointing maps, as the unpacking only
// updated "downward" pointing maps
//        domain.m_feNodeManager.AddToNodeToElementMap(
// domain.m_feElementManager, elementLists );
//        domain.m_feFaceManager.AddToFaceToElementMap(
// domain.m_feElementManager, elementLists );
//
//
//        //Faces
//        {
//          lArray1d newFaceIndices(numNewFaces);
//          for(unsigned kf(firstFaceIndex), ii(0); kf <
// m_domain->m_feFaceManager.m_numFaces;++kf,++ii){
//            newFaceIndices[ii] = kf;
//          }
//
//          domain.m_feNodeManager.AddToVariableOneToManyFromInverse(
// "nodeToFaceMap", domain.m_feFaceManager.m_toNodesRelation,  newFaceIndices);
//          domain.m_feEdgeManager.AddToVariableOneToManyFromInverse(
// "edgesToFaces", domain.m_feFaceManager.m_toEdgesRelation, newFaceIndices);
//        }
//
//        //Edges
//        {
//          lArray1d newEdgeIndices(numNewEdges);
//          for(unsigned ke(firstEdgeIndex), ii(0); ke <
// m_domain->m_feEdgeManager.DataLengths();++ke,++ii){
//            newEdgeIndices[ii] = ke;
//          }
//
//          domain.m_feNodeManager.AddToNodeToEdgeMap( domain.m_feEdgeManager,
// newEdgeIndices );
//        }
//
//      }
//
//
//
//    }
//  }
//
//}

/**
 * @author walsh24
 * Labels ghostRank of objects on the same partition
 *
 **/
//void SpatialPartition::SetSinglePartitionGhostArrays(PhysicalDomainT& domain )
//{
//  // check for single partition
//  bool isSinglePartition = false;
//  for(unsigned int p =0; p < m_Partitions.size();++p){
//    if(m_Partitions[p]==1){
//      isSinglePartition = true;
//      break;
//    }
//  }
//
//  if(isSinglePartition){
//
//    const localIndex n = NeighborCommunication::NumberOfSyncNames();
//    array<PhysicalDomainT::ObjectDataStructureKeys> objectNames;
//    NeighborCommunication::SyncNames(objectNames);
//
//    for(localIndex i = 0; i < n; ++i){
//
//      if(objectNames[i] != PhysicalDomainT::FiniteElementElementManager)
//      {
//        array<Field<FieldInfo::ghostRank>::Type>& ghostRankCurr
//            = domain.GetObjectDataStructure(objectNames[i]).GetFieldData<FieldInfo::ghostRank>();
//        ObjectDataStructureBaseT& object =
//  domain.GetObjectDataStructure(objectNames[i]);
//
//        gArray1d& localToGlobalMap = object.m_localToGlobalMap;
//        std::map<globalIndex,localIndex>& globalToLocalMap =
// object.m_globalToLocalMap;
//
//        for( localIndex l=0 ; l!= localToGlobalMap.size() ; ++l )
//        {
//          globalIndex g = localToGlobalMap[l];
//          localIndex lb = globalToLocalMap[g];
//          // if map from local -> global -> local does not return the same
// object, it is a ghost.
//          if(lb != l){
//            ghostRankCurr[l] = m_rank;
//            ghostRankCurr[lb] = -1;
//             m_localGhostSources[objectNames[i]].push_back(lb);
//             m_localGhosts[objectNames[i]].push_back(l);
//
//          }
//        }
//      }
//      else
//      {
//        //--per usual: elements are treated as a special case
//        for( std::map< ElementManagerT::RegKeyType, ElementRegionT >::iterator
// iregion=domain.m_feElementManager.m_ElementRegions.begin() ;
//            iregion!=domain.m_feElementManager.m_ElementRegions.end() ;
// ++iregion )
//        {
//          const std::string& elemRegionName = iregion->first;
//          ElementRegionT& elemRegion = iregion->second;
//          array<Field<FieldInfo::ghostRank>::Type>& ghostRankCurr =
// elemRegion.GetFieldData<FieldInfo::ghostRank>();
//
//
//          gArray1d& localToGlobalMap = elemRegion.m_localToGlobalMap;
//          std::map<globalIndex,localIndex>&  globalToLocalMap =
// elemRegion.m_globalToLocalMap;
//
//          for( localIndex l=0 ; l!= ghostRankCurr.size() ; ++l )
//          {
//            globalIndex g = localToGlobalMap[l];
//            localIndex lb = globalToLocalMap[g];
//            // if map from local -> global -> local does not return the same
// object, it is a ghost.
//            if(lb != l){
//              ghostRankCurr[l] = m_rank;
//              ghostRankCurr[lb] = -1;
//              m_elementRegionsLocalGhostSources[elemRegionName].push_back(lb);
//              m_elementRegionsLocalGhosts[elemRegionName].push_back(l);
//
//            }
//          }
//        }
//      }
//
//    }
//  }
//}
//
///**
// * @author walsh24
// *
// *
//**/
//void
// SpatialPartition::CorrectReferencePositionsForPeriodicBoundaries(PhysicalDomainT&
// domain){
//  array< R1Tensor >& refPositions = *(domain.m_feNodeManager.m_refposition);
//
//  for(unsigned int i =0; i < m_periodicSets.size();++i){
//
//    int dimension = m_periodicSets[i].m_dimension;
//
//    if( (m_coords[dimension] == 0)  ||
//        (m_coords[dimension] == (m_Partitions[dimension]-1) ) )
//    {
//
////      svector& setnames = m_periodicSets[i].m_setNames;
////      set<localIndex>* theSets[2];
////      theSets[0] = &(domain.m_feNodeManager.m_Sets[setnames[0]]);
////      theSets[1] = &(domain.m_feNodeManager.m_Sets[setnames[1]]);
//
//      realT myMax = m_max[dimension];
//      realT myMin = m_min[dimension];
//      realT gridLength = m_gridSize[dimension];
//
//      if(m_Partitions[dimension]>1){
//
//        // Multiple partitions
//        //--------------------
//
//        if( m_coords[dimension] == 0)  // lower boundary
//        {
//
//          for( localIndex ii =0; ii < refPositions.size() ; ++ii )
//          {
//            if( refPositions[ii][dimension] > myMax + 1e-16){
//
//              realT oldPos = refPositions[ii][dimension];
//              realT newPos = oldPos - gridLength;
//
//              // This is a bit of a hack
//              // The periodic boundary conditions are enforced if shifting the
// ghost node will bring it closer to the
//              // minimum edge of the partition than it is to the maximum edge.
//              // Needed if there are two partitions along one axis to
// distinguish ghost nodes on maximum and minimum borders.
//              if( oldPos-myMax > myMin-newPos) refPositions[ii][dimension] =
// newPos;
//            }
//          }
//         } else { // upper boundary
//
//          for( localIndex ii =0; ii < refPositions.size() ; ++ii )
//          {
//            if( refPositions[ii][dimension] < myMin - 1e-16){
//
//              realT oldPos = refPositions[ii][dimension];
//              realT newPos = oldPos + gridLength;
//
//              // This is a bit of a hack
//              // The periodic boundary conditions are enforced if shifting the
// ghost node will bring it closer to the
//              // maximum edge of the partition than it is to the minimum edge.
//              // Needed if there are two partitions along one axis to
// distinguish ghost nodes on maximum and minimum borders.
//              if(myMin-oldPos > newPos-myMax) refPositions[ii][dimension] =
// newPos;
//
//            }
//          }
//        }
//
//      } else {
//        // Empty - Single partition node reference positions are modified on
// creation
//      }
//
//    } // is boundary partition
//  } // loop over boundary sets
//}

void SpatialPartition::AddNeighbors(const unsigned int idim,
                                    MPI_Comm& cartcomm,
                                    int* ncoords)
{

  if (idim == nsdof)
  {
    bool me = true;
    for ( unsigned int i = 0 ; i < nsdof ; i++)
    {
      if (ncoords[i] != this->m_coords(i))
      {
        me = false;
        break;
      }
    }
    if (!me)
    {
      m_neighbors.push_back(NeighborCommunicator());
      int rank;
      MPI_Cart_rank(cartcomm, ncoords, &rank);
      m_neighbors.back().SetNeighborRank( rank );
//      m_neighbors.back().Initialize( rank, this->m_rank, this->m_size );

//      array<int> nbrcoords(nsdof);
//      for(unsigned int i =0; i < nsdof; ++i) nbrcoords[i] = ncoords[i];
//      neighborCommPtrIndx[nbrcoords] = m_neighbors.size()-1;
    }
  }
  else
  {
    const int dim = this->m_Partitions(idim);
    const bool periodic = this->m_Periodic(idim);
    for (int i = -1 ; i < 2 ; i++)
    {
      ncoords[idim] = this->m_coords(idim) + i;
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

void SpatialPartition::AddNeighborsMetis(gSet& neighborList)
{
  gSet::iterator itNeighbor = neighborList.begin();
  for ( ; itNeighbor != neighborList.end() ; itNeighbor++)
  {
    m_neighbors.push_back(NeighborCommunicator());
    m_neighbors.back().SetNeighborRank( integer_conversion<int>(*itNeighbor) );

//    m_neighbors.back().Initialize( integer_conversion<int>(*itNeighbor), this->m_rank, this->m_size );
  }
}

/**
 * <SpatialPartition xpar="1" ypar="3" zpar="1"
 *                   periodic_y_axis="1">
 *   <PeriodicSets>
 *      <PeriodicSet setnames="Ymax Ymin" dimension="1" />
 *   </PeriodicSets>
 *  </SpatialPartition>
 **/
//void SpatialPartition::ReadXMLInput(TICPP::HierarchicalDataNode& hdn)
//{
//
//  m_ghostDepth = hdn.GetAttributeOrDefault<unsigned int>("ghostDepth",1);
//  //R.W. Read Metis partitioning information
//  //Currently implemented as an alternative to the Cartesian system
//  //if <SpatialPartition Metispar=numofpar> is given, will use Metis instead
// of Cartesian
//  //Should not effect the use of previous implementations
//  if (hdn.HasAttribute("Metispar"))
//  {
//    m_sizeMetis = hdn.GetAttributeOrDefault<unsigned int>("Metispar",1);
//    //assign default values for m_Partitions and m_Periodic
//
//  }
//  else
//  {
//
//  if(m_Partitions(0) * m_Partitions(1) * m_Partitions(2) == 0){ // i.e. has
// not previously been set
//    m_Partitions(0) = hdn.GetAttributeOrDefault<unsigned int>("xpar",1);
//    m_Partitions(1) = hdn.GetAttributeOrDefault<unsigned int>("ypar",1);
//    m_Partitions(2) = hdn.GetAttributeOrDefault<unsigned int>("zpar",1);
//    m_size = m_Partitions(0) * m_Partitions(1) * m_Partitions(2);
//  } else {
//	// add partition data to hdn (to allow xml to be written to file after
// overridden from commandline)
//	hdn.AddAttributePair("xpar", toString( m_Partitions(0)) );
//	hdn.AddAttributePair("ypar", toString( m_Partitions(1)) );
//	hdn.AddAttributePair("zpar", toString( m_Partitions(2)) );
//  }
//
//    if( m_Partitions(0) > 1 )
//    {
//      m_PartitionLocations[0] = hdn.GetAttributeVector<realT>("xloc",",");
//      if( !(m_PartitionLocations[0].empty()) && (m_Partitions(0)-1) !=
// static_cast<int>( m_PartitionLocations[0].size() ) )
//      {
//#ifdef USE_ATK
//        SLIC_ERROR( "SpatialPartition::ReadXMLInput(): number of x-partition
// locations does not equal number of partitions - 1\n");
//#endif
//      }
//    }
//    if( m_Partitions(1) > 1 )
//    {
//      m_PartitionLocations[1] = hdn.GetAttributeVector<realT>("yloc",",");
//      if( !(m_PartitionLocations[1].empty()) && (m_Partitions(1)-1) !=
// static_cast<int>( m_PartitionLocations[1].size() ) )
//      {
//#ifdef USE_ATK
//        SLIC_ERROR( "SpatialPartition::ReadXMLInput(): number of y-partition
// locations does not equal number of partitions - 1\n");
//#endif
//      }
//    }
//    if( m_Partitions(2) > 1 )
//    {
//      m_PartitionLocations[2] = hdn.GetAttributeVector<realT>("zloc",",");
//      if( !(m_PartitionLocations[2].empty()) && (m_Partitions(2)-1) !=
// static_cast<int>( m_PartitionLocations[2].size() ) )
//      {
//#ifdef USE_ATK
//        SLIC_ERROR( "SpatialPartition::ReadXMLInput(): number of z-partition
// locations does not equal number of partitions - 1\n");
//#endif
//      }
//    }
//
//
//    // periodic boundaries
//    m_Periodic(0) = hdn.GetAttributeOrDefault<unsigned
// int>("periodic_x_axis",0);
//    m_Periodic(1) = hdn.GetAttributeOrDefault<unsigned
// int>("periodic_y_axis",0);
//    m_Periodic(2) = hdn.GetAttributeOrDefault<unsigned
// int>("periodic_z_axis",0);
//
//    TICPP::HierarchicalDataNode* periodicSetsNode =
// hdn.GetChild("PeriodicSets");
//    if (periodicSetsNode != NULL)
//    {
//      if (m_rank  == 0) std::cout << "Reading PeriodicSets:" << std::endl;
//      for (TICPP::HierarchicalDataNode* psNode = periodicSetsNode->Next(true);
// psNode; psNode
//      = periodicSetsNode->Next())
//      {
//        PeriodicSet periodicSet;
//        periodicSet.ReadXML(*psNode);
//        m_periodicSets.push_back(periodicSet);
//
//        if( m_Periodic( periodicSet.m_dimension ) != 1 ){
//#ifdef USE_ATK
//          SLIC_ERROR( "SpatialPartition::ReadXMLInput(): Periodic set
// requested for non-periodic dimension " + toString(periodicSet.m_dimension) +
// " \n");
//#endif
//        }
//      }
//    }
//  }
//
//  const realT bufferSize =
// hdn.GetAttributeOrDefault<realT>("contactGhostBuffer",0.0);
//  SetContactGhostRange( bufferSize );
//}
//
//void SpatialPartition::getSizes(R1Tensor& min, R1Tensor& max) const
//{
//  min = m_gridMin;
//  max = m_gridMax;
//}
//void SpatialPartition::getPartitionSizes(R1Tensor& min, R1Tensor& max) const
//{
//  min = m_min;
//  max = m_max;
//}
//void SpatialPartition::getPartitionGeometricalBoundary(R1Tensor& min,
// R1Tensor& max) const
//{
//  min = m_min;
//  max = m_max;
//}
//
//void SpatialPartition::UpdatePartitionBoundingBox(NodeManager& nodeManager)
//{
//  const array<integer>& isGhost =
// nodeManager.GetFieldData<FieldInfo::ghostRank>();
//  m_xBoundingBoxMin = (*nodeManager.m_refposition)[0];
//  m_xBoundingBoxMin += (*nodeManager.m_displacement)[0];
//  m_xBoundingBoxMax = m_xBoundingBoxMin;
//
//  for (localIndex i = 1; i != nodeManager.DataLengths(); ++i)
//  {
//    if (isGhost[i] < 0)
//    {
//      R1Tensor x = (*nodeManager.m_refposition)[i];
//      x += (*nodeManager.m_displacement)[i];
//      for (localIndex dim = 0; dim < nsdof; ++dim)
//      {
//        m_xBoundingBoxMin[dim] = std::min(m_xBoundingBoxMin[dim],x[dim]);
//        m_xBoundingBoxMax[dim] = std::max(m_xBoundingBoxMax[dim],x[dim]);
//      }
//    }
//  }
//}

void SpatialPartition::GetPartitionBoundingBox(R1Tensor& xmin, R1Tensor& xmax )
{
  xmin = m_xBoundingBoxMin;
  xmax = m_xBoundingBoxMax;
}

/**
 * @param min global minimum spatial dimensions
 * @param max global maximum spatial dimensions
 **/
void SpatialPartition::setSizes( const R1Tensor& min, const R1Tensor& max )
{

  {
    //get size of problem and decomposition
    MPI_Comm_size(MPI_COMM_WORLD, &m_size);

    //check to make sure our dimensions agree
    {
      int check = 1;
      for (unsigned int i = 0 ; i < nsdof ; i++)
      {
        check *= this->m_Partitions(i);
      }
      assert(check == m_size);
    }

    //get communicator, rank, and coordinates
    MPI_Comm cartcomm;
    {
      int reorder = 0;
      MPI_Cart_create(MPI_COMM_WORLD, nsdof, m_Partitions.data(), m_Periodic.data(), reorder, &cartcomm);
    }
    MPI_Comm_rank(cartcomm, &m_rank);
    MPI_Cart_coords(cartcomm, m_rank, nsdof, m_coords.data());


    m_color = GetColor();

    //add neighbors
    {
      int ncoords[nsdof];
      m_neighbors.clear();
      AddNeighbors(0, cartcomm, ncoords);
    }

    MPI_Comm_free(&cartcomm);

  }



  // global values
  m_gridMin = min;
  m_gridMax = max;
  m_gridSize = max;
  m_gridSize -= min;

  // block values
  m_blockSize = m_gridSize;

  m_min = min;
  for( unsigned int i=0 ; i<nsdof ; ++i )
  {
    const int nloc = m_Partitions(i) - 1;
    const localIndex nlocl = static_cast<localIndex>(nloc);
    if( m_PartitionLocations[i].empty() )
    {
      // the default "even" spacing
      m_blockSize(i) /= m_Partitions(i);
      m_min(i) += m_coords(i) * m_blockSize(i);
      m_max(i) = min(i) + (m_coords(i) + 1) * m_blockSize(i);

      m_PartitionLocations[i].resize(nlocl);
      localIndex j = 0;
      for(array<real64>::iterator it = m_PartitionLocations[i].begin() ; it != m_PartitionLocations[i].end() ; ++it, ++j)
      {
        *it = (j+1) * m_blockSize(i);
      }
    }
    else if( nlocl == m_PartitionLocations[i].size() )
    {
      const int parIndex = m_coords[i];
      if( parIndex == 0 )
      {
        m_min[i] = min[i];
        m_max[i] = m_PartitionLocations[i][parIndex];
      }
      else if( parIndex == nloc )
      {
        m_min[i] = m_PartitionLocations[i][parIndex-1];
        m_max[i] = max[i];
      }
      else
      {
        m_min[i] = m_PartitionLocations[i][parIndex-1];
        m_max[i] = m_PartitionLocations[i][parIndex];
      }
    }
    else
    {
#ifdef USE_ATK
      SLIC_ERROR( "SpatialPartition::setSizes(): number of partition locations does not equal number of partitions - 1\n");
#endif
    }
  }
}

void SpatialPartition::setGlobalDomainSizes( const R1Tensor& min, const R1Tensor& max )
{
  // global values
  // without updating partition sizes.  We need this in mesh generator when we
  // have extension zones.
  m_gridMin = min;
  m_gridMax = max;
  m_gridSize = max;
  m_gridSize -= min;
}

void SpatialPartition::SetPartitionGeometricalBoundary( R1Tensor& min, R1Tensor& max )
{
  // We need this in mesh generator when we have extension zones.
  m_min = min;
  m_max = max;
}


bool SpatialPartition::IsCoordInPartition( const realT& coord, const int dir )
{
  bool rval = true;
  const int i = dir;
  if(m_Periodic(i))
  {
    if(m_Partitions(i) != 1 )
    {
      realT localCenter = MapValueToRange(coord,m_gridMin(i),m_gridMax(i));
      rval = rval && localCenter >= m_min(i) && localCenter < m_max(i);
    }

  }
  else
  {
    rval = rval && (m_Partitions(i)==1 || (coord >= m_min(i) && coord < m_max(i)));
  }

  return rval;
}

bool SpatialPartition::IsCoordInPartition( const R1Tensor& elemCenter )
{
  bool rval = true;
  for(unsigned int i = 0 ; i < nsdof ; i++)
  {
    if(m_Periodic(i))
    {

      if(m_Partitions(i) != 1 )
      {
        realT localCenter = MapValueToRange(elemCenter(i),m_gridMin(i),m_gridMax(i));
        rval = rval && localCenter >= m_min(i) && localCenter < m_max(i);
      }

    }
    else
    {
      rval = rval && (m_Partitions(i)==1 || (elemCenter(i) >= m_min(i) && elemCenter(i) < m_max(i)));
    }
  }
  return rval;
}

bool SpatialPartition::IsCoordInPartition( const R1Tensor& elemCenter, const int numDistPartition  )
{
  bool rval = true;
  R1Tensor m_xBoundingBoxMinTemp, m_xBoundingBoxMaxTemp;
  for(unsigned int i = 0 ; i < nsdof ; i++)
  {
    m_xBoundingBoxMinTemp(i) = m_min(i) - numDistPartition*m_blockSize(i);
    m_xBoundingBoxMaxTemp(i) = m_max(i) + numDistPartition*m_blockSize(i);
  }

  for(unsigned int i = 0 ; i < nsdof ; i++)
  {
    if(m_Periodic(i))
    {

      if(m_Partitions(i) != 1 )
      {
        realT localCenter = MapValueToRange(elemCenter(i),m_gridMin(i),m_gridMax(i));
        rval = rval && localCenter >= m_xBoundingBoxMinTemp(i) && localCenter <= m_xBoundingBoxMaxTemp(i);
      }

    }
    else
    {

      rval = rval && (m_Partitions(i)==1 || (elemCenter(i) >= m_xBoundingBoxMinTemp(i) && elemCenter(i) <= m_xBoundingBoxMaxTemp(i)));
    }
  }
  return rval;
}

bool SpatialPartition::IsCoordInPartitionClosed( const R1Tensor& elemCenter )
// A variant with intervals closed at both ends
{
  bool rval = true;
  for(unsigned int i = 0 ; i < nsdof ; i++)
  {
    if(m_Periodic(i))
    {

      if(m_Partitions(i) != 1 )
      {
        realT localCenter = MapValueToRange(elemCenter(i),m_gridMin(i),m_gridMax(i));
        rval = rval && localCenter >= m_min(i) && localCenter < m_max(i);
      }

    }
    else
    {
      rval = rval && (m_Partitions(i)==1 || (elemCenter(i) >= m_min(i) && elemCenter(i) <= m_max(i)));
    }
  }
  return rval;
}

bool SpatialPartition::IsCoordInPartitionBoundingBox( const R1Tensor& elemCenter )

{
  bool rval = true;
  for(unsigned int i = 0 ; i < nsdof ; i++)
  {
    if(m_Periodic(i))
    {

      if(m_Partitions(i) != 1 )
      {
        realT localCenter = MapValueToRange(elemCenter(i),m_gridMin(i),m_gridMax(i));
        rval = rval && localCenter >= m_xBoundingBoxMin(i) && localCenter <= m_xBoundingBoxMax(i);
      }

    }
    else
    {
      rval = rval && (m_Partitions(i)==1 || (elemCenter(i) >= m_xBoundingBoxMin(i) && elemCenter(i) <= m_xBoundingBoxMax(i)));
    }
  }
  return rval;
}



void SpatialPartition::SetContactGhostRange( const realT bufferSize )
{
  m_contactGhostMin = m_min;
  m_contactGhostMin -= bufferSize;

  m_contactGhostMax = m_max;
  m_contactGhostMax += bufferSize;
}

bool SpatialPartition::IsCoordInContactGhostRange( const R1Tensor& elemCenter )
{
  bool rval = true;
  for(unsigned int i = 0 ; i < nsdof ; i++)
  {
    if(m_Periodic(i))
    {
      if(m_Partitions(i) != 1 )
      {

        realT minBuffer = m_min(i)-m_contactGhostMin(i);
        realT maxBuffer = m_contactGhostMax(i)-m_max(i);
        realT localCenterA = MapValueToRange(elemCenter(i),m_gridMin(i)-minBuffer,m_gridMax(i)-minBuffer);
        realT localCenterB = MapValueToRange(elemCenter(i),m_gridMin(i)+maxBuffer,m_gridMax(i)+maxBuffer);

        rval = rval && (m_Partitions(i)==1 || (localCenterA >= m_contactGhostMin(i) && localCenterB < m_contactGhostMax(i)));
      }
    }
    else
    {
      rval = rval && (m_Partitions(i)==1 || (elemCenter(i) >= m_contactGhostMin(i) && elemCenter(i) < m_contactGhostMax(i)));
    }
  }
  return rval;
}

//
//void SpatialPartition::WriteSiloDerived( SiloFile& siloFile )
//{
//  siloFile.DBWriteWrapper("m_Partitions",m_Partitions);
//  siloFile.DBWriteWrapper("m_Periodic",m_Periodic);
//  siloFile.DBWriteWrapper("m_coords",m_coords);
//  siloFile.DBWriteWrapper("m_PartitionLocations0",m_PartitionLocations[0]);
//  siloFile.DBWriteWrapper("m_PartitionLocations1",m_PartitionLocations[1]);
//  siloFile.DBWriteWrapper("m_PartitionLocations2",m_PartitionLocations[2]);
//  siloFile.DBWriteWrapper("m_blockSize",m_blockSize);
//  siloFile.DBWriteWrapper("m_min",m_min);
//  siloFile.DBWriteWrapper("m_max",m_max);
//  siloFile.DBWriteWrapper("m_gridSize",m_gridSize);
//  siloFile.DBWriteWrapper("m_gridMin",m_gridMin);
//  siloFile.DBWriteWrapper("m_gridMax",m_gridMax);
//
//}
//
//void SpatialPartition::ReadSiloDerived( const SiloFile& siloFile )
//{
//  siloFile.DBReadWrapper("m_Partitions",m_Partitions);
//  siloFile.DBReadWrapper("m_Periodic",m_Periodic);
//  siloFile.DBReadWrapper("m_coords",m_coords);
//  siloFile.DBReadWrapper("m_PartitionLocations0",m_PartitionLocations[0]);
//  siloFile.DBReadWrapper("m_PartitionLocations1",m_PartitionLocations[1]);
//  siloFile.DBReadWrapper("m_PartitionLocations2",m_PartitionLocations[2]);
//  siloFile.DBReadWrapper("m_blockSize",m_blockSize);
//  siloFile.DBReadWrapper("m_min",m_min);
//  siloFile.DBReadWrapper("m_max",m_max);
//  siloFile.DBReadWrapper("m_gridSize",m_gridSize);
//  siloFile.DBReadWrapper("m_gridMin",m_gridMin);
//  siloFile.DBReadWrapper("m_gridMax",m_gridMax);
//
//
//}
}
