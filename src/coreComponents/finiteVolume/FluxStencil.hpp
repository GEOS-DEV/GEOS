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

/**
 * @file FluxStencil.hpp
 *
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXSTENCIL_HPP_
#define SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXSTENCIL_HPP_

#include <common/DataTypes.hpp>
#include "codingUtilities/Utilities.hpp"

namespace geosx
{

/**
 * @class FluxStencil
 *
 * Class representing a spatial discretization stencil.
 */
template <typename INDEX, typename WEIGHT>
class FluxStencil
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

  explicit FluxStencil();

  explicit FluxStencil( localIndex const numConn,
                        localIndex const avgStencilSize );

  /// return the size of the stencil collection (i.e. number of connections)
  localIndex numConnections() const;

  /// resize the collection
  void reserve(localIndex const numConn,
               localIndex const avgStencilSize);

  /// add data for one connection
  void add( localIndex const numPts,
            INDEX  const * const indices,
            WEIGHT const * const weights,
            localIndex const connectorIndex );

  /// zero out connections
  bool zero( localIndex const connectorIndex );

  /// called after adding connections is done to compress the data and release unused memory
  void compress();

  struct Entry
  {
    INDEX  index;
    WEIGHT weight;
  };

  ArrayOfArraysView<Entry const, true> getConnections() const { return m_connections; }

  const localIndex_array &getConnectorMeshIndices() const { return m_connectorMeshIndices; }
  
private:

  ArrayOfArrays<Entry> m_connections;
  map<localIndex, localIndex> m_connectorIndices;
  localIndex_array m_connectorMeshIndices; 
};

template<typename INDEX, typename WEIGHT>
FluxStencil<INDEX, WEIGHT>::FluxStencil()
  : FluxStencil( 0, 0 )
{

}

template<typename INDEX, typename WEIGHT>
FluxStencil<INDEX, WEIGHT>::FluxStencil( localIndex const numConn,
                                                     localIndex const avgStencilSize )
  : m_connections()
{
  reserve(numConn, avgStencilSize);
}

template<typename INDEX, typename WEIGHT>
localIndex FluxStencil<INDEX, WEIGHT>::numConnections() const
{
  return m_connections.size();
}

template<typename INDEX, typename WEIGHT>
void FluxStencil<INDEX, WEIGHT>::reserve( localIndex const numConn,
                                          localIndex const avgStencilSize )
{
  m_connections.reserve( numConn );
  m_connections.reserveValues( numConn * avgStencilSize );
}

template<typename INDEX, typename WEIGHT>
void FluxStencil<INDEX, WEIGHT>::add( localIndex const numPts,
                                      INDEX  const * const indices,
                                      WEIGHT const * const weights,
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
  m_connectorMeshIndices.push_back(connectorIndex);  
}

template<typename INDEX, typename WEIGHT>
bool FluxStencil<INDEX, WEIGHT>::zero( localIndex const connectorIndex )
{
  return
  executeOnMapValue( m_connectorIndices, connectorIndex, [&]( localIndex const connectionListIndex )
  {
    Entry * const entries = m_connections[connectionListIndex];
    for (localIndex i = 0; i < m_connections.sizeOfArray( connectionListIndex ); ++i)
    {
      entries[i].weight = 0; // TODO remove entries altogether?
    }
//    m_connections.resizeArray( connectionListIndex, 0 );
  });
}


template<typename INDEX, typename WEIGHT>
void FluxStencil<INDEX, WEIGHT>::compress()
{
  // nothing for the moment
}

}


#endif //SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXSTENCIL_HPP_
