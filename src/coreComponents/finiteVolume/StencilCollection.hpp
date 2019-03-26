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
 * @file StencilCollection.hpp
 *
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXSTENCIL_HPP_
#define SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXSTENCIL_HPP_

#include <common/DataTypes.hpp>
#include <rajaInterface/GEOS_RAJA_Interface.hpp>

namespace geosx
{

/**
 * @class StencilCollection
 *
 * Class representing a spatial discretization stencil.
 * The idea is to hide implementation (i.e. memory layout) behind a lambda-based interface
 * which evaluates user-provided functions on the stencil
 *
 * TODO:
 * - group by by stencil size and/or element region for efficient execution
 */
template <typename INDEX, typename WEIGHT>
class StencilCollection
{
public:

  /**
   * @brief Number of points the flux is between (normally 2)
   */
  static localIndex constexpr NUM_POINT_IN_FLUX = 2;

  /**
   * @brief Maximum number of points in a stencil (required to use static arrays in kernels)
   */
  static localIndex constexpr MAX_STENCIL_SIZE = 18;

  // provide aliases for template type parameters
  using index_type  = INDEX;
  using weight_type = WEIGHT;

  /**
   * @struct Accessor
   *
   * A wrapper of a single connection stencil
   * Hides implementation details of stencil representation
   * (e.g. array of structs vs struct of arrays vs CSR-like storage)
   */
  class Accessor;

  explicit StencilCollection();

  explicit StencilCollection( localIndex numConn, localIndex avgStencilSize );

  /// return the size of the stencil collection (i.e. number of connections)
  localIndex numConnections() const;

  /// resize the collection
  void reserve(localIndex numConn, localIndex avgStencilSize);

  /// add data for one connection
  void add( localIndex const numPts,
            INDEX  const * indices,
            WEIGHT const * weights,
            localIndex const connectorIndex );

  /// zero out connections
  void zero( localIndex const connectorIndex,
             INDEX const cells[2] );

  /// called after adding connections is done to compress the data and release unused memory
  void compress();

  /// evaluate a user function on each connection
  template <typename POLICY=stencilPolicy, typename LAMBDA>
  void forAll( LAMBDA && lambda ) const;

  Accessor operator[]( localIndex iconn ) const;

  struct Entry
  {
    INDEX  index;
    WEIGHT weight;
  };

  csArrayView2d<Entry const> getConnections() const { return m_connections; }

private:

  csArray2d<Entry> m_connections;
  map<localIndex, localIndex> m_connectorIndices;

};

// *** implementation ***

template <typename INDEX, typename WEIGHT>
class StencilCollection<INDEX, WEIGHT>::Accessor
{
public:

  Accessor( csArrayView2d<Entry const> const & connections, localIndex const index )
    : m_size( connections.size(index) ),
      m_entries( connections[index] )
  {}

  /// return the stencil size
  localIndex size() const { return m_size; }

  /// return the point index of connected point i
  INDEX index( localIndex const i ) const { return m_entries[i].index; }

  /// apply a user-defined function on the two connected cells only
  template <typename LAMBDA>
  void forConnected( LAMBDA && lambda ) const
  {
    for (localIndex i = 0; i < StencilCollection<INDEX, WEIGHT>::NUM_POINT_IN_FLUX; ++i)
    {
      lambda( m_entries[i].index, i );
    }
  }

  /// apply a user-defined function on the stencil
  template <typename LAMBDA>
  void forAll( LAMBDA && lambda ) const
  {
    for (localIndex i = 0; i < size(); ++i)
    {
      lambda( m_entries[i].index, m_entries[i].weight, i );
    }
  }

private:

  localIndex m_size;
  StencilCollection<INDEX, WEIGHT>::Entry const * m_entries;
};

template<typename INDEX, typename WEIGHT>
StencilCollection<INDEX, WEIGHT>::StencilCollection()
  : StencilCollection( 0, 0 )
{

}

template<typename INDEX, typename WEIGHT>
StencilCollection<INDEX, WEIGHT>::StencilCollection( localIndex numConn, localIndex avgStencilSize )
  : m_connections()
{
  reserve(numConn, avgStencilSize);
}

template<typename INDEX, typename WEIGHT>
localIndex StencilCollection<INDEX, WEIGHT>::numConnections() const
{
  return m_connections.size();
}

template<typename INDEX, typename WEIGHT>
void StencilCollection<INDEX, WEIGHT>::reserve( localIndex numConn, localIndex avgStencilSize )
{
  m_connections.reserveNumArrays( numConn );
  m_connections.reserveValues( numConn * avgStencilSize );
}

template<typename INDEX, typename WEIGHT>
template<typename POLICY, typename LAMBDA>
void StencilCollection<INDEX, WEIGHT>::forAll( LAMBDA && lambda ) const
{
  csArrayView2d<Entry const> const & connections = getConnections();
  forall_in_range<POLICY>( 0, numConnections(), GEOSX_LAMBDA ( localIndex const index )
  {
    lambda( Accessor( connections, index ), index );
  });
}

template<typename INDEX, typename WEIGHT>
void StencilCollection<INDEX, WEIGHT>::add( localIndex const numPts,
                                                    INDEX  const * indices,
                                                    WEIGHT const * weights,
                                                    localIndex const connectorIndex )
{
  GEOS_ERROR_IF( numPts >= MAX_STENCIL_SIZE, "Maximum stencil size exceeded" );

  stackArray1d<Entry, MAX_STENCIL_SIZE> entries(numPts);
  for (localIndex i = 0; i < numPts; ++i)
  {
    entries[i] = { indices[i], weights[i] };
  }

  m_connections.appendArray( entries.data(), numPts );
  m_connectorIndices[connectorIndex] = m_connections.size() - 1;
}

template<typename INDEX, typename WEIGHT>
void StencilCollection<INDEX, WEIGHT>::zero( localIndex const connectorIndex,
                                                     INDEX const cells[2] )
{
  localIndex const connectionListIndex = m_connectorIndices.at( connectorIndex );

  Entry * const entries = m_connections[connectionListIndex];

  if( ( entries[0].index == cells[0] && entries[1].index == cells[1] ) ||
      ( entries[0].index == cells[1] && entries[1].index == cells[0] ) )
  {
    for (localIndex i = 0; i < m_connections.size(connectionListIndex); ++i)
    {
      entries[i].weight = 0; // TODO remove entries altogether?
    }
  }
}


template<typename INDEX, typename WEIGHT>
void StencilCollection<INDEX, WEIGHT>::compress()
{
  // nothing for the moment
}

template<typename INDEX, typename WEIGHT>
typename StencilCollection<INDEX, WEIGHT>::Accessor
StencilCollection<INDEX, WEIGHT>::operator[](localIndex iconn) const
{
  return Accessor( m_connections, iconn );
}

}


#endif //SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXSTENCIL_HPP_
