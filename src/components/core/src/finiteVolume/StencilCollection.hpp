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
 * @class FluxStencil
 *
 * Class representing a spatial discretization stencil.
 * The idea is to hide implementation (i.e. memory layout) behind a lambda-based interface
 * which evaluates user-provided functions on the stencil
 *
 * TODO:
 * - optimize for contiguous memory layout
 * - group by by stencil size and/or element region for efficient execution
 */
class StencilCollection
{
public:

  explicit StencilCollection() : m_faceConnectors() {}

  /**
 * @struct A structure containing a single cell (element) identifier triplet
 */
  struct CellDescriptor
  {
    localIndex region;
    localIndex subRegion;
    localIndex index;
  };

  /**
   * @struct A structure describing a single (generally multi-point) FV connection stencil
   */
  struct CellConnection
  {
    CellDescriptor        connectedCellIndices[2]; ///< identifiers of connected cells
    array<CellDescriptor> stencilCellIndices;      ///< identifiers of cells in stencil
    array<real64>         stencilWeights;          ///< stencil weights (e.g. transmissibilities)

    void resize(localIndex const size) { stencilCellIndices.resize(size);
                                         stencilWeights.resize(size);    }
  };

  /**
   * @struct A wrapper of a single connection stencil
   * Hides implementation details of stencil representation
   * (e.g. array of structs vs struct of arrays vs CSR-like storage)
   */
  class Accessor
  {
  public:

    /// return the stencil size
    localIndex size() const { return m_conn.stencilCellIndices.size(); }

    /// apply a user-defined function on the two connected cells only
    template <typename LAMBDA>
    void forConnected(LAMBDA && lambda) const
    {
      for (localIndex i = 0; i < 2; ++i)
        lambda(m_conn.connectedCellIndices[i].region,
               m_conn.connectedCellIndices[i].subRegion,
               m_conn.connectedCellIndices[i].index,
               i);
    }

    /// apply a user-defined function on the stencil
    template <typename LAMBDA>
    void forAll(LAMBDA && lambda) const
    {
      for (localIndex i = 0; i < size(); ++i)
        lambda(m_conn.stencilCellIndices[i].region,
               m_conn.stencilCellIndices[i].subRegion,
               m_conn.stencilCellIndices[i].index,
               m_conn.stencilWeights[i], i);
    }

  private:

    friend class StencilCollection;

    Accessor(StencilCollection const & stencil, localIndex index)
      : m_conn(stencil.m_faceConnectors[index]) {}

    CellConnection const & m_conn;

  };

  /// return the size of the stencil collection (i.e. number of connections)
  localIndex numConnections() const { return m_faceConnectors.size(); }

  /// resize the collection
  void resize(localIndex newsize) { m_faceConnectors.resize(newsize); }

  // TODO: this interface smells, it's temporary
  /// set data for one connection
  void set(localIndex index, CellDescriptor connCells[2], array<CellDescriptor> cells, array<real64> weights);

  /// evaluate a user function on each connection (the function receives a StencilAccessor object)
  template <typename POLICY=stencilPolicy, typename LAMBDA=void>
  void forAll(LAMBDA && lambda) const
  {
    RAJA::RangeSegment seg(0, numConnections());
    RAJA::forall<POLICY>(seg, [&] (localIndex index) -> void
    {
      lambda(Accessor(*this, index));
    });
  }

private:

  friend class StencilView;

  /// array that holds the list of CellConnection objects.
  array<CellConnection> m_faceConnectors;

};

}


#endif //SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXSTENCIL_HPP_
