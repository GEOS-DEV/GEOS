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
            INDEX  const * const elementRegionIndices,
            INDEX  const * const elementSubRegionIndices,
            INDEX  const * const elementIndices,
            WEIGHT const * const weights,
            localIndex const connectorIndex );

  /// zero out connections
  bool zero( localIndex const connectorIndex );

  /// called after adding connections is done to compress the data and release unused memory
  void compress();

  ArrayOfArraysView<INDEX const, true> getElementRegionIndices() const { return m_elementRegionIndices; }
  ArrayOfArraysView<INDEX const, true> getElementSubRegionIndices() const { return m_elementRegionIndices; }
  ArrayOfArraysView<INDEX const, true> getElementIndices() const { return m_elementIndices; }
  ArrayOfArraysView<WEIGHT const, true> getWeights() const { return m_weights; }

private:

  ArrayOfArrays<INDEX>  m_elementRegionIndices;
  ArrayOfArrays<INDEX>  m_elementSubRegionIndices;
  ArrayOfArrays<INDEX>  m_elementIndices;
  ArrayOfArrays<WEIGHT> m_weights;

  
  map<localIndex, localIndex> m_connectorIndices;

};

template<typename INDEX, typename WEIGHT>
FluxStencil<INDEX, WEIGHT>::FluxStencil()
  : FluxStencil( 0, 0 )
{

}

template<typename INDEX, typename WEIGHT>
FluxStencil<INDEX, WEIGHT>::FluxStencil( localIndex const numConn,
                                         localIndex const avgStencilSize ):
  m_elementRegionIndices(),
  m_elementSubRegionIndices(),
  m_elementIndices(),
  m_weights(),
  m_connectorIndices()
{
  reserve(numConn, avgStencilSize);
}

template<typename INDEX, typename WEIGHT>
localIndex FluxStencil<INDEX, WEIGHT>::numConnections() const
{
  return m_weights.size();
}

template<typename INDEX, typename WEIGHT>
void FluxStencil<INDEX, WEIGHT>::reserve( localIndex const numConn,
                                          localIndex const avgStencilSize )
{
  m_elementRegionIndices.reserve( numConn );
  m_elementSubRegionIndices.reserve( numConn );
  m_elementIndices.reserve( numConn );
  m_weights.reserve( numConn );

  m_elementRegionIndices.reserveValues( numConn * avgStencilSize );
  m_elementSubRegionIndices.reserveValues( numConn * avgStencilSize );
  m_elementIndices.reserveValues( numConn * avgStencilSize );
  m_weights.reserveValues( numConn * avgStencilSize );
}

template<typename INDEX, typename WEIGHT>
void FluxStencil<INDEX, WEIGHT>::add( localIndex const numPts,
                                      INDEX  const * const elementRegionIndices,
                                      INDEX  const * const elementSubRegionIndices,
                                      INDEX  const * const elementIndices,
                                      WEIGHT const * const weights,
                                      localIndex const connectorIndex )
{
  GEOS_ERROR_IF( numPts >= MAX_STENCIL_SIZE, "Maximum stencil size exceeded" );

  m_elementRegionIndices.appendArray( elementRegionIndices, numPts );
  m_elementSubRegionIndices.appendArray( elementSubRegionIndices, numPts );
  m_elementIndices.appendArray( elementIndices, numPts );
  m_weights.appendArray( weights, numPts );

  m_connectorIndices[connectorIndex] = m_weights.size() - 1;
}

template<typename INDEX, typename WEIGHT>
bool FluxStencil<INDEX, WEIGHT>::zero( localIndex const connectorIndex )
{
  return
  executeOnMapValue( m_connectorIndices, connectorIndex, [&]( localIndex const connectionListIndex )
  {
    for (localIndex i = 0; i < m_weights.sizeOfArray( connectionListIndex ); ++i)
    {
      m_weights[connectionListIndex][i] = 0; // TODO remove entries altogether?
    }
//    m_connections.resizeArray( connectionListIndex, 0 );
  });
}


template<typename INDEX, typename WEIGHT>
void FluxStencil<INDEX, WEIGHT>::compress()
{
  // nothing for the moment
}


/**
 * @struct A structure containing a single cell (element) identifier triplet
 */
struct CellDescriptor
{
  localIndex region;
  localIndex subRegion;
  localIndex index;

  bool operator==( CellDescriptor const & other )
  {
    return( region==other.region && subRegion==other.subRegion && index==other.index );
  }
};

/**
 * @struct A structure describing an arbitrary point participating in a stencil
 *
 * Nodal and face center points are identified by local mesh index.
 * Cell center points are identified by a triplet <region,subregion,index>.
 *
 * The sad reality is, a boundary flux MPFA stencil may be comprised of a mix of
 * cell and face centroids, so we have to discriminate between them at runtime
 */
struct PointDescriptor
{
  enum class Tag { CELL, FACE, NODE };

  Tag tag;

  union
  {
    localIndex nodeIndex;
    localIndex faceIndex;
    CellDescriptor cellIndex;
  };
};


}


#endif //SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXSTENCIL_HPP_
