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
#include "finiteVolume/StencilCollection.hpp"
#include "managers/DomainPartition.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
  static constexpr auto FVstencil = "FVstencil";
}
}

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
  enum class Tag { CELL, FACE, NODE, PERF };

  Tag tag;

  union
  {
    localIndex nodeIndex;
    localIndex faceIndex;
    localIndex perfIndex;
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
  using CellStencil = StencilCollection<CellDescriptor, real64>;
  using FaceStencil = StencilCollection<PointDescriptor, real64>;
  using WellStencil = StencilCollection<PointDescriptor, real64>;

  FluxApproximationBase() = delete;

  FluxApproximationBase(string const & name, dataRepository::ManagedGroup * const parent);

  /// provides const access to the cell stencil collection
  CellStencil const & getCellStencil() const;

  /// provides access to the cell stencil collection
  CellStencil & getCellStencil();

  /// return a boundary stencil by face set name
  FaceStencil const & getFaceStencil(string const & setName) const;

  /// return a boundary stencil by face set name
  FaceStencil & getFaceStencil(string const & setName);

  /// check if a stencil exists
  bool hasFaceStencil(string const & setName) const;

  /// return a boundary stencil by face set name
  WellStencil const & getWellStencil(string const & wellName) const;

  /// return a boundary stencil by face set name
  WellStencil & getWellStencil(string const & wellName);

  /// check if a stencil exists
  bool hasWellStencil(string const & wellName) const;

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

    dataRepository::ViewKey fieldName         = { fieldNameString };
    dataRepository::ViewKey boundaryFieldName = { boundaryFieldNameString };
    dataRepository::ViewKey coeffName         = { coeffNameString };
    dataRepository::ViewKey cellStencil       = { cellStencilString };

  } viewKeysFABase;

  struct groupKeyStruct
  {

    static constexpr auto faceStencilsString = "faceStencils";
    static constexpr auto wellStencilsString = "wellStencils";

    dataRepository::ViewKey faceStencils = { faceStencilsString };
    dataRepository::ViewKey wellStencils = { wellStencilsString };

  } groupKeysFABase;

  /// actual computation of the cell-to-cell stencil, to be overridden by implementations
  virtual void computeCellStencil(DomainPartition const * domain,
                                  CellStencil & stencil) const = 0;

  /// actual computation of the boundary stencil, to be overridden by implementations
  virtual void computeFaceStencil( DomainPartition const * domain,
                                   set<localIndex> const & faceSet,
                                   FaceStencil & stencil ) const = 0;

  /// actual computation of well-to-cell stencil, to be overridden by implementations
  virtual void computeWellStencil( DomainPartition const * domain,
                                   WellBase const * well,
                                   WellStencil & stencil ) const = 0;

protected:

  void InitializePostInitialConditions_PreSubGroups(ManagedGroup * const rootGroup) override;

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
  this->GetGroup( groupKeysFABase.faceStencils )->forViewWrappersByType<FaceStencil>([&] (auto const & vw) -> void
  {
    if (vw.getName() != viewKeyStruct::cellStencilString)
      lambda(vw.reference());
  });
}

} // namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXAPPROXIMATIONBASE_HPP_
