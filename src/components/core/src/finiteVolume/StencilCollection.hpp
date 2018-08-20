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

/*
 * @file StencilCollection.hpp
 *
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXSTENCIL_HPP_
#define SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXSTENCIL_HPP_

#include <common/DataTypes.hpp>
#include <rajaInterface/GEOS_RAJA_Policies.hpp>

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
 * - optimize for contiguous memory layout
 * - group by by stencil size and/or element region for efficient execution
 */
template <typename IndexType, typename WeightType>
class StencilCollection
{
public:

  // provide aliases for template type parameters
  using index_type = IndexType;
  using weight_type = WeightType;

  /**
   * @struct Accessor
   *
   * A wrapper of a single connection stencil
   * Hides implementation details of stencil representation
   * (e.g. array of structs vs struct of arrays vs CSR-like storage)
   */
  class Accessor;

  explicit StencilCollection();

  explicit StencilCollection(localIndex numConn, localIndex avgStencilSize);

  /// return the size of the stencil collection (i.e. number of connections)
  localIndex numConnections() const;

  /// resize the collection
  void reserve(localIndex numConn, localIndex avgStencilSize);

  /// add data for one connection
  void add(IndexType const cells[2],
           array1d<IndexType> const & stencilCells,
           array1d<WeightType> const & stencilWeights);

  /// called after adding connections is done to compress the data and release unused memory
  void compress();

  /// evaluate a user function on each connection
  template <typename POLICY=stencilPolicy, typename LAMBDA=void>
  void forAll(LAMBDA && lambda) const;

  Accessor operator[](localIndex iconn) const;

private:

  /**
   * @struct A structure describing a single (generally multi-point) FV connection stencil
   */
  struct Connection
  {
    IndexType         connectedPointIndices[2]; ///< identifiers of connected points
    array1d<IndexType>  stencilPointIndices;      ///< identifiers of points in stencil
    array1d<WeightType> stencilWeights;           ///< stencil weights (e.g. transmissibilities)

    void resize(localIndex const size) { stencilPointIndices.resize(size);
      stencilWeights.resize(size);    }
  };

  /// array that holds the list of CellConnection objects.
  array1d<Connection> m_connectionList;

};

// *** implementation ***

template <typename IndexType, typename WeightType>
class StencilCollection<IndexType, WeightType>::Accessor
{
public:

  Accessor(StencilCollection<IndexType, WeightType> const & stencil, localIndex index)
    : m_conn(stencil.m_connectionList[index]) {}

  /// return the stencil size
  localIndex size() const { return m_conn.stencilPointIndices.size(); }

  /// apply a user-defined function on the two connected cells only
  template <typename LAMBDA>
  void forConnected(LAMBDA && lambda) const
  {
    for (localIndex i = 0; i < 2; ++i)
      lambda(m_conn.connectedPointIndices[i], i);
  }

  /// apply a user-defined function on the stencil
  template <typename LAMBDA>
  void forAll(LAMBDA && lambda) const
  {
    for (localIndex i = 0; i < size(); ++i)
      lambda(m_conn.stencilPointIndices[i], m_conn.stencilWeights[i], i);
  }

private:

  StencilCollection<IndexType, WeightType>::Connection const & m_conn;

};

template<typename IndexType, typename WeightType>
StencilCollection<IndexType, WeightType>::StencilCollection()
  : StencilCollection(0, 0)
{

}

template<typename IndexType, typename WeightType>
StencilCollection<IndexType, WeightType>::StencilCollection(localIndex numConn, localIndex avgStencilSize)
  : m_connectionList()
{
  reserve(numConn, avgStencilSize);
}

template<typename IndexType, typename WeightType>
localIndex StencilCollection<IndexType, WeightType>::numConnections() const
{
  return m_connectionList.size();
}

template<typename IndexType, typename WeightType>
void StencilCollection<IndexType, WeightType>::reserve(localIndex numConn, localIndex avgStencilSize)
{
  m_connectionList.reserve(numConn);
}

template<typename IndexType, typename WeightType>
template<typename POLICY, typename LAMBDA>
void StencilCollection<IndexType, WeightType>::forAll(LAMBDA &&lambda) const
{
  RAJA::RangeSegment seg(0, numConnections());
  RAJA::forall<POLICY>(seg, [=] (localIndex index) mutable -> void
  {
    lambda(Accessor(*this, index));
  });
}

template<typename IndexType, typename WeightType>
void StencilCollection<IndexType, WeightType>::add(IndexType const cells[2],
                                                   const array1d<IndexType> & stencilCells,
                                                   const array1d<WeightType> & stencilWeights)
{
  Connection conn = { { cells[0], cells[1] }, stencilCells, stencilWeights };
  m_connectionList.push_back(conn);
}

template<typename IndexType, typename WeightType>
void StencilCollection<IndexType, WeightType>::compress()
{
  // nothing for the moment
}

template<typename IndexType, typename WeightType>
typename StencilCollection<IndexType, WeightType>::Accessor
StencilCollection<IndexType, WeightType>::operator[](localIndex iconn) const
{
  return Accessor(*this, iconn);
}

}


#endif //SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXSTENCIL_HPP_
