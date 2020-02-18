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
 * @file FluxApproximationBase.hpp
 *
 */

#ifndef GEOSX_FINITEVOLUME_FLUXAPPROXIMATIONBASE_HPP_
#define GEOSX_FINITEVOLUME_FLUXAPPROXIMATIONBASE_HPP_

#include "dataRepository/Group.hpp"
#include "finiteVolume/FluxStencil.hpp"
#include "CellElementStencilTPFA.hpp"
#include "FaceElementStencil.hpp"
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
class FluxApproximationBase : public dataRepository::Group
{
public:

  // necessary declarations for factory instantiation of derived classes
  using CatalogInterface = dataRepository::CatalogInterface<FluxApproximationBase, string const &, Group * const >;
  static typename CatalogInterface::CatalogType& GetCatalog();

  // typedefs for stored stencil types
  using BoundaryStencil = FluxStencil<PointDescriptor, real64>;

  FluxApproximationBase() = delete;

  FluxApproximationBase(string const & name, dataRepository::Group * const parent);

  /// return a boundary stencil by face set name
  BoundaryStencil const & getBoundaryStencil(string const & setName) const;

  /// return a boundary stencil by face set name
  BoundaryStencil & getBoundaryStencil(string const & setName);

  /// check if a stencil exists
  bool hasBoundaryStencil(string const & setName) const;

  /// call a user-provided function for each boundary stencil
  template<typename LAMBDA>
  void forCellStencils(LAMBDA && lambda) const;


  template<typename TYPE, typename ... TYPES, typename LAMBDA>
  void forStencils(LAMBDA && lambda) const;

  /// call a user-provided function for each boundary stencil
  template<typename LAMBDA>
  void forBoundaryStencils(LAMBDA && lambda) const;

  /// triggers computation of the stencil, implemented in derived classes
  void compute( DomainPartition const & domain );

  virtual void addToFractureStencil( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                     string const & GEOSX_UNUSED_PARAM( faceElementRegionName ),
                                     bool const GEOSX_UNUSED_PARAM(initFlag) ) {}


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

  string_array const & targetRegions() const { return m_targetRegions; }
  string_array &       targetRegions()       { return m_targetRegions; }

protected:

  virtual void InitializePostInitialConditions_PreSubGroups( Group * const rootGroup ) override;

  /// actual computation of the cell-to-cell stencil, to be overridden by implementations
  virtual void computeCellStencil( DomainPartition const & domain ) = 0;

  /// actual computation of the boundary stencil, to be overridden by implementations
  virtual void computeBoundaryStencil( DomainPartition const & domain,
                                       SortedArray<localIndex> const & faceSet,
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
//TODO remove dependence on CellElementStencilTPFA and FaceElementStencil
  this->forWrappers<CellElementStencilTPFA,FaceElementStencil>([&] (auto const * const wrapper) -> void
  {
    lambda(wrapper->reference());
  });
}

template<typename TYPE, typename ... TYPES, typename LAMBDA>
void FluxApproximationBase::forStencils(LAMBDA && lambda) const
{
  this->forWrappers<TYPE,TYPES...>([&] (auto const * const wrapper) -> void
  {
    lambda(wrapper->reference());
  });
}


template<typename LAMBDA>
void FluxApproximationBase::forBoundaryStencils(LAMBDA && lambda) const
{
  this->forWrappers<BoundaryStencil>([&] (auto const * const wrapper) -> void
  {
    lambda(wrapper->reference());
  });
}

} // namespace geosx

#endif //GEOSX_FINITEVOLUME_FLUXAPPROXIMATIONBASE_HPP_
