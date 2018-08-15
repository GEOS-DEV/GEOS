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
 * @file FluxApproximationBase.hpp
 *
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXAPPROXIMATIONBASE_HPP_
#define SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXAPPROXIMATIONBASE_HPP_

#include "dataRepository/ManagedGroup.hpp"
#include "finiteVolume/StencilCollection.hpp"
#include "managers/DomainPartition.hpp"

namespace geosx
{

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

/**
 * @class FluxApproximationBase
 *
 * Base class for various flux approximation classes.
 * Stores the main and boundary stencils, construction is implemented in derived classes.
 * Main stencil is the one for cell-to-cell fluxes.
 * Boundary stencils are for Dirichlet boundary conditions
 */
class FluxApproximationBase : public dataRepository::ManagedGroup
{
public:

  // necessary declarations for factory instantiation of derived classes
  using CatalogInterface = cxx_utilities::CatalogInterface<FluxApproximationBase, string const &, ManagedGroup * const >;
  static typename CatalogInterface::CatalogType& GetCatalog();

  // typedefs for stored stencil types
  using CellStencil     = StencilCollection<CellDescriptor, real64>;
  using BoundaryStencil = StencilCollection<PointDescriptor, real64>;

  void FillDocumentationNode() override;

  FluxApproximationBase() = delete;

  FluxApproximationBase(string const & name, dataRepository::ManagedGroup * const parent);

  /// provides const access to the cell stencil collection
  CellStencil const & getStencil() const;

  /// provides access to the cell stencil collection
  CellStencil & getStencil();

  /// return a boundary stencil by face set name
  BoundaryStencil const & getBoundaryStencil(string const & setName) const;

  /// return a boundary stencil by face set name
  BoundaryStencil & getBoundaryStencil(string const & setName);

  /// check if a stencil exists
  bool hasBoundaryStencil(string const & setName) const;

  /// call a user-provided function for each boundary stencil
  template<typename LAMBDA>
  void forBoundaryStencils(LAMBDA && lambda) const;

  /// triggers computation of the stencil, implemented in derived classes
  void compute(DomainPartition * domain);

  struct viewKeyStruct
  {
    static constexpr auto fieldNameString         = "fieldName";
    static constexpr auto boundaryFieldNameString = "boundaryFieldName";
    static constexpr auto coeffNameString         = "coefficientName";
    static constexpr auto cellStencilString       = "cellStencil";
  };

  struct groupKeyStruct
  {
    static constexpr auto boundarySetDataString = "BoundarySetData";
  };

protected:

  /// actual computation of the cell-to-cell stencil, to be overridden by implementations
  virtual void computeMainStencil(DomainPartition * domain, CellStencil & stencil) = 0;

  /// actual computation of the boundary stencil, to be overridden by implementations
  virtual void computeBoundaryStencil(DomainPartition * domain, set<localIndex> const & faceSet, BoundaryStencil & stencil) = 0;

  /// pointer to boundary set manager
  dataRepository::ManagedGroup * m_boundarySetData;

  /// name of the primary solution field
  string m_fieldName;

  /// name of the boundary field (used to filter boundary conditions)
  string m_boundaryFieldName;

  /// name of the coefficient field
  string m_coeffName;

};

template<typename LAMBDA>
void FluxApproximationBase::forBoundaryStencils(LAMBDA && lambda) const
{
  this->forViewWrappersByType<BoundaryStencil>([&] (auto const & vw) -> void
  {
    if (vw.getName() != viewKeyStruct::cellStencilString)
      lambda(vw.reference());
  });
}

} // namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXAPPROXIMATIONBASE_HPP_
