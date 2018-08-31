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
 * @file PartitionBase.h
 * @author settgast1
 * @date Mar 16, 2011
 */

#ifndef PARTITIONBASE_H_
#define PARTITIONBASE_H_

#include "common/DataTypes.hpp"
#include <mpi.h>
#include "NeighborCommunicator.hpp"


class oBinStream;
class iBinStream;


namespace geosx
{

namespace dataRepository
{
class ManagedGroup;
}

class DomainPartition;
class ObjectManagerBase;

class PartitionBase
{

public:

  virtual ~PartitionBase();


  void SetDomain( DomainPartition * domain );


  virtual void InitializePostSubGroups( dataRepository::ManagedGroup * const ) = 0;


  virtual bool IsCoordInPartition( const R1Tensor& elemCenter ) = 0;
  virtual bool IsCoordInPartition( const R1Tensor& elemCenter,
                                   const int numDistPartition ) = 0;
  virtual bool IsCoordInPartition( const realT& coord, const int dir ) = 0;

  virtual void setSizes( const R1Tensor& min, const R1Tensor& max ) = 0;

  virtual void setPartitions( unsigned int xPartitions,
                              unsigned int yPartitions,
                              unsigned int zPartitions ) = 0;

  virtual bool IsCoordInContactGhostRange( const R1Tensor& elemCenter ) = 0;

//  virtual void ReadXML( xmlWrapper::xmlNode const & targetNode ) = 0;

  //virtual void AssignGlobalIndices( DomainPartition * domain );

//  virtual void FindMatchedBoundaryIndices( string const & key,
//                                           const ObjectManagerBase& object );


//  virtual void SetUpNeighborLists( DomainPartition * domain,
//                                   const bool contactActive );

  void SetRankOfNeighborNeighbors();

//  virtual void ResetNeighborLists( PhysicalDomainT& domain,
//                                   const int elementGhostingDepth );

//  virtual void ModifyGhostsAndNeighborLists( const ModifiedObjectLists& modifiedObjects );

  template< typename T >
  void SendReceive( const array1d<array1d<T> >& sendArray, array1d<array1d<T> >& recvArray );

//  void SynchronizeFields( const std::map<std::string, string_array >& fieldNames,
//                          const CommRegistry::commID commID = CommRegistry::genericComm01 );

  void SetOwnedByRank( const std::map< std::string, globalIndex_array>& localBoundaryGlobalIndices,
                       std::map<std::string, std::map< globalIndex, int > >& boundaryOwnership );

  void SetGhostArrays( DomainPartition * domain );

  localIndex_array GetFaceSendIndices();

  virtual void SetContactGhostRange( const double bufferSize ) = 0;
//
//  void SetBufferSizes( const std::map<string, string_array >& fieldNames,
//                       const CommRegistry::commID commID  );
//
//  int NumberOfNeighbors( ) {return integer_conversion<int>(m_neighbors.size());}

  int m_size;
  int m_sizeMetis;
  int m_rank;

  virtual int GetColor() = 0;

  int Color() const {return m_color;}
  int NumColor() const {return m_numColors;}

//  void WriteSilo( SiloFile& siloFile );
//
//  void ReadSilo( const SiloFile& siloFile );
  void DeleteExcessNeighbors();
  void GraphBasedColoring();

protected:
  PartitionBase();
  PartitionBase( const unsigned int numPartitions, const unsigned int thisPartiton );
//
  array1d<NeighborCommunicator> m_neighbors;

  array1d<MPI_Request> m_mpiRequest;
  array1d<MPI_Status> m_mpiStatus;

  R1Tensor m_contactGhostMin;
  R1Tensor m_contactGhostMax;

  int m_color;
  int m_numColors;

  DomainPartition * const m_domain;

public:
  realT m_t1;
  realT m_t2;
  realT m_t3;
  realT m_t4;

  bool m_hasLocalGhosts;
  std::map<std::string, localIndex_array> m_localGhosts;
  std::map< std::string, localIndex_array> m_elementRegionsLocalGhosts;

  std::map<std::string, localIndex_array> m_localGhostSources;
  std::map< std::string, localIndex_array> m_elementRegionsLocalGhostSources;

  int m_ghostDepth;

private:
//  virtual void AssignGlobalIndices( ObjectDataStructureBaseT& object, const ObjectDataStructureBaseT&
// compositionObject );

  void CommunicateRequiredObjectIndices();

//  virtual void WriteSiloDerived( SiloFile& siloFile ) = 0;
//
//  virtual void ReadSiloDerived( const SiloFile& siloFile ) = 0;


};

}

#endif /* PARTITIONBASE_H_ */
