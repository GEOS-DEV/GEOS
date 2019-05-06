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
 * @file FluxApproximationBase.hpp
 *
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXAPPROXIMATIONBASE_HPP_
#define SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXAPPROXIMATIONBASE_HPP_

#include "dataRepository/ManagedGroup.hpp"
#include "finiteVolume/FluxStencil.hpp"
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
  using CellStencil     = FluxStencil<CellDescriptor, real64>;
  using BoundaryStencil = FluxStencil<PointDescriptor, real64>;

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
  void forCellStencils(LAMBDA && lambda) const;

  /// call a user-provided function for each boundary stencil
  template<typename LAMBDA>
  void forBoundaryStencils(LAMBDA && lambda) const;

  /// triggers computation of the stencil, implemented in derived classes
  void compute( DomainPartition const & domain );

  struct viewKeyStruct
  {
    static constexpr auto fieldNameString             = "fieldName";
    static constexpr auto boundaryFieldNameString     = "boundaryFieldName";
    static constexpr auto coeffNameString             = "coefficientName";
    static constexpr auto targetRegionsString         = "targetRegions";
    static constexpr auto areaRelativeToleranceString = "areaRelTol";

    static constexpr auto cellStencilString           = "cellStencil";
    static constexpr auto fractureStencilString       = "fractureStencil";
  };

  struct groupKeyStruct
  {
  };

protected:

  virtual void InitializePostInitialConditions_PreSubGroups( ManagedGroup * const rootGroup ) override;

  /// actual computation of the cell-to-cell stencil, to be overridden by implementations
  virtual void computeCellStencil( DomainPartition const & domain,
                                   CellStencil & stencil ) = 0;

  virtual void computeFractureStencil( DomainPartition const & domain,
                                       CellStencil & fractureStencil,
                                       CellStencil & cellStencil ) = 0;

  /// actual computation of the boundary stencil, to be overridden by implementations
  virtual void computeBoundaryStencil( DomainPartition const & domain,
                                       set<localIndex> const & faceSet,
                                       BoundaryStencil & stencil ) = 0;

  /// name of the primary solution field
  string m_fieldName;

  /// name of the boundary field (used to filter boundary conditions)
  string m_boundaryFieldName;

  /// name of the coefficient field
  string m_coeffName;

  /// names of target regions to build the stencil for
  string_array m_targetRegions;

  /// relative tolerance
  real64 m_areaRelTol;

};

template<typename LAMBDA>
void FluxApproximationBase::forCellStencils(LAMBDA && lambda) const
{
  this->forViewWrappersByType<CellStencil>([&] (auto const & vw) -> void
  {
    lambda(vw.reference());
  });
}

template<typename LAMBDA>
void FluxApproximationBase::forBoundaryStencils(LAMBDA && lambda) const
{
  this->forViewWrappersByType<BoundaryStencil>([&] (auto const & vw) -> void
  {
    lambda(vw.reference());
  });
}

} // namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXAPPROXIMATIONBASE_HPP_
