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
 * @file FluxStencil.hpp
 *
 */

#ifndef GEOSX_FINITEVOLUME_FLUXSTENCIL_HPP_
#define GEOSX_FINITEVOLUME_FLUXSTENCIL_HPP_

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

private:

  ArrayOfArrays<Entry> m_connections;
  map<localIndex, localIndex> m_connectorIndices;

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
  GEOSX_ERROR_IF( numPts >= MAX_STENCIL_SIZE, "Maximum stencil size exceeded" );

  stackArray1d<Entry, MAX_STENCIL_SIZE> entries(numPts);
  for (localIndex i = 0; i < numPts; ++i)
  {
    entries[i] = { indices[i], weights[i] };
  }

  m_connections.appendArray( entries.data(), numPts );
  m_connectorIndices[connectorIndex] = m_connections.size() - 1;
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


#endif //GEOSX_FINITEVOLUME_FLUXSTENCIL_HPP_
