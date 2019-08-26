/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
 * CommunicationTools.hpp
 *
 *  Created on: Jan 6, 2018
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MPI_COMMUNICATIONS_COMMUNICATIONTOOLS_HPP_
#define SRC_COMPONENTS_CORE_SRC_MPI_COMMUNICATIONS_COMMUNICATIONTOOLS_HPP_

#include "common/DataTypes.hpp"
#include "mpi.h"
#include <set>
namespace geosx
{


class ObjectManagerBase;
class NeighborCommunicator;
class MeshLevel;
class ElementRegionManager;

class MPI_iCommData;


class CommunicationTools
{
public:

  enum class Reduction { Sum, Min, Max };

  CommunicationTools();
  ~CommunicationTools();

  template< typename T >
  static MPI_Datatype getMpiType();

  static MPI_Op getMpiOp( Reduction const op );

  static void AssignGlobalIndices( ObjectManagerBase & object,
                                   ObjectManagerBase const & compositionObject,
                                   array1d<NeighborCommunicator> & neighbors );

  static void AssignNewGlobalIndices( ObjectManagerBase & object,
                                      std::set<localIndex> const & indexList );

  static void
  AssignNewGlobalIndices( ElementRegionManager & elementManager,
                          std::map< std::pair<localIndex,localIndex>, std::set<localIndex> > const & newElems );

  static void FindGhosts( MeshLevel * const meshLevel,
                          array1d<NeighborCommunicator> & neighbors );

  static int MPI_Size( MPI_Comm const & comm );
  static int MPI_Rank( MPI_Comm const & comm );

  static std::set<int> & getFreeCommIDs();
  static int reserveCommID();
  static void releaseCommID( int & ID );

  static void FindMatchedPartitionBoundaryObjects( ObjectManagerBase * const group,
                                                   array1d<NeighborCommunicator> & allNeighbors );

  static void SynchronizeFields( const std::map<string, string_array >& fieldNames,
                                 MeshLevel * const mesh,
                                 array1d<NeighborCommunicator> & allNeighbors );

  static void SynchronizePackSendRecvSizes( const std::map<string, string_array >& fieldNames,
                                            MeshLevel * const mesh,
                                            array1d<NeighborCommunicator> & neighbors,
                                            MPI_iCommData & icomm );

  static void SynchronizePackSendRecv( const std::map<string, string_array >& fieldNames,
                                       MeshLevel * const mesh,
                                       array1d<NeighborCommunicator> & allNeighbors,
                                       MPI_iCommData & icomm );

  static void SynchronizeUnpack( MeshLevel * const mesh,
                                 array1d<NeighborCommunicator> & neighbors,
                                 MPI_iCommData & icomm );

  template<typename T>
  static void allGather( T const myValue, array1d<T> & allValues );

  /**
   * @brief Compute exclusive prefix sum and full sum
   * @tparam T type of local (rank) value
   * @tparam U type of global (sum) value
   * @param value the local value
   * @return a pair where first is the prefix sum, second is the full sum
   */
  template< typename U, typename T >
  static std::pair<U, U> PrefixSum( T const & value );

  template<typename T>
  static void Broadcast( T & value, int srcRank = 0 );

  template< typename T >
  static T Reduce( T const & value, Reduction const op);

  template<typename T>
  static T Sum( T const & value );

  template<typename T>
  static T Min( T const & value );

  template<typename T>
  static T Max( T const & value );
};


class MPI_iCommData
{
public:

  MPI_iCommData():
    size(0),
    commID(-1),
    sizeCommID(-1),
    fieldNames(),
    mpiSendBufferRequest(),
    mpiRecvBufferRequest(),
    mpiSendBufferStatus(),
    mpiRecvBufferStatus()
  {
    commID = CommunicationTools::reserveCommID();
    sizeCommID = CommunicationTools::reserveCommID();
  }

  ~MPI_iCommData()
  {
    if( commID >= 0 )
    {
      CommunicationTools::releaseCommID(commID);
    }

    if( sizeCommID >= 0 )
    {
      CommunicationTools::releaseCommID(sizeCommID);
    }

  }

  void resize( localIndex numMessages )
  {
    mpiSendBufferRequest.resize( numMessages );
    mpiRecvBufferRequest.resize( numMessages );
    mpiSendBufferStatus.resize( numMessages );
    mpiRecvBufferStatus.resize( numMessages );
    mpiSizeSendBufferRequest.resize( numMessages );
    mpiSizeRecvBufferRequest.resize( numMessages );
    mpiSizeSendBufferStatus.resize( numMessages );
    mpiSizeRecvBufferStatus.resize( numMessages );
    size = static_cast<int>(numMessages);
  }

  int size;
  int commID;
  int sizeCommID;
  std::map<string, string_array > fieldNames;

  array1d<MPI_Request> mpiSendBufferRequest;
  array1d<MPI_Request> mpiRecvBufferRequest;
  array1d<MPI_Status>  mpiSendBufferStatus;
  array1d<MPI_Status>  mpiRecvBufferStatus;

  array1d<MPI_Request> mpiSizeSendBufferRequest;
  array1d<MPI_Request> mpiSizeRecvBufferRequest;
  array1d<MPI_Status>  mpiSizeSendBufferStatus;
  array1d<MPI_Status>  mpiSizeRecvBufferStatus;
};


template<> inline MPI_Datatype CommunicationTools::getMpiType<double>()         { return MPI_DOUBLE; }
template<> inline MPI_Datatype CommunicationTools::getMpiType<int>()            { return MPI_INT; }
template<> inline MPI_Datatype CommunicationTools::getMpiType<long int>()       { return MPI_LONG; }
template<> inline MPI_Datatype CommunicationTools::getMpiType<long long int>()  { return MPI_LONG_LONG; }

inline MPI_Op CommunicationTools::getMpiOp( Reduction const op )
{
  switch( op )
  {
    case Reduction::Sum:
      return MPI_SUM;
    case Reduction::Min:
      return MPI_MIN;
    case Reduction::Max:
      return MPI_MAX;
    default:
      GEOS_ERROR( "Unsupported reduction op" );
  }
}

template<typename T>
void CommunicationTools::allGather( T const myValue, array1d<T> & allValues )
{
#ifdef GEOSX_USE_MPI
  int const mpiSize = MPI_Size( MPI_COMM_GEOSX );
  allValues.resize( mpiSize );

  MPI_Datatype const MPI_TYPE = getMpiType<T>();

  MPI_Allgather( &myValue, 1, MPI_TYPE, allValues.data(), 1, MPI_TYPE, MPI_COMM_GEOSX );

#else
  allValues.resize(1);
  allValues[0] = myValue;
#endif
}

template< typename U, typename T >
std::pair<U, U> CommunicationTools::PrefixSum( T const & value )
{
  array1d<T> gather;
  allGather( value, gather );

  std::pair<U, U> rval( 0, 0 );

  int const mpiSize = MPI_Size( MPI_COMM_GEOSX );
  int const mpiRank = MPI_Rank( MPI_COMM_GEOSX );

  for( localIndex p = 0 ; p < mpiSize ; ++p )
  {
    if( p < mpiRank )
    {
      rval.first += static_cast<U>( gather[p] );
    }
    rval.second += static_cast<U>( gather[p] );
  }

  return rval;
}

template<typename T>
void CommunicationTools::Broadcast( T & value, int srcRank )
{
#ifdef GEOSX_USE_MPI
  MPI_Datatype const mpiType = getMpiType<T>();
  MPI_Bcast( &value, 1, mpiType, srcRank, MPI_COMM_GEOSX );
#endif
}

template< typename T >
T CommunicationTools::Reduce( T const & value, Reduction const op )
{
  T result = value;
#ifdef GEOSX_USE_MPI
  MPI_Allreduce( &value, &result, 1, getMpiType<T>(), getMpiOp( op ), MPI_COMM_GEOSX );
#endif
  return result;
}

template< typename T >
T CommunicationTools::Sum( T const & value )
{
  return CommunicationTools::Reduce( value, Reduction::Sum );
}

template< typename T >
T CommunicationTools::Min( T const & value )
{
  return CommunicationTools::Reduce( value, Reduction::Min );
}

template< typename T >
T CommunicationTools::Max( T const & value )
{
  return CommunicationTools::Reduce( value, Reduction::Max );
}

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MPI_COMMUNICATIONS_COMMUNICATIONTOOLS_HPP_ */
