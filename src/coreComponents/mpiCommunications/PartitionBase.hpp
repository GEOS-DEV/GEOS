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
 * @file PartitionBase.hpp
 */

#ifndef GEOSX_MPICOMMUNICATIONS_PARTITIONBASE_HPP_
#define GEOSX_MPICOMMUNICATIONS_PARTITIONBASE_HPP_

#include "common/DataTypes.hpp"

#include "mpiCommunications/NeighborCommunicator.hpp"


class oBinStream;
class iBinStream;


namespace geosx
{

namespace dataRepository
{
class Group;
}

class DomainPartition;
class ObjectManagerBase;

class PartitionBase
{

public:

  virtual ~PartitionBase();


  void SetDomain( DomainPartition * domain );


  virtual bool IsCoordInPartition( const R1Tensor & elemCenter ) = 0;
  virtual bool IsCoordInPartition( const R1Tensor & elemCenter,
                                   const int numDistPartition ) = 0;
  virtual bool IsCoordInPartition( const realT & coord, const int dir ) = 0;

  virtual void setSizes( const R1Tensor & min, const R1Tensor & max ) = 0;

  virtual void setPartitions( unsigned int xPartitions,
                              unsigned int yPartitions,
                              unsigned int zPartitions ) = 0;

  virtual bool IsCoordInContactGhostRange( const R1Tensor & elemCenter ) = 0;

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
  void SendReceive( const array1d< array1d< T > > & sendArray, array1d< array1d< T > > & recvArray );

//  void SynchronizeFields( const std::map<std::string, string_array >& fieldNames,
//                          const CommRegistry::commID commID = CommRegistry::genericComm01 );

  void SetOwnedByRank( const std::map< std::string, globalIndex_array > & localBoundaryGlobalIndices,
                       std::map< std::string, std::map< globalIndex, int > > & boundaryOwnership );

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

  virtual void InitializePostSubGroups( dataRepository::Group * const ) = 0;

//
  std::vector< NeighborCommunicator > m_neighbors;

  array1d< MPI_Request > m_mpiRequest;
  array1d< MPI_Status > m_mpiStatus;

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
  std::map< std::string, localIndex_array > m_localGhosts;
  std::map< std::string, localIndex_array > m_elementRegionsLocalGhosts;

  std::map< std::string, localIndex_array > m_localGhostSources;
  std::map< std::string, localIndex_array > m_elementRegionsLocalGhostSources;

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

#endif /* GEOSX_MPICOMMUNICATIONS_PARTITIONBASE_HPP_ */
