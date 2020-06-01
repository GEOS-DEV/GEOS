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
 * @struct CellDescriptor
 * @brief A structure containing a single cell (element) identifier triplet.
 */
struct CellDescriptor
{
  /// region index
  localIndex region;
  /// subregion index
  localIndex subRegion;
  /// cell index
  localIndex index;

  /**
   * @brief Comparison operator between two CellDescriptors.
   * @param[in] other the CellDescriptor to compare with
   * @return true if they represent the same mesh element
   */
  bool operator==( CellDescriptor const & other )
  {
    return( region==other.region && subRegion==other.subRegion && index==other.index );
  }
};

/**
 * @struct PointDescriptor
 * @brief A structure describing an arbitrary point participating in a stencil.
 *
 * Nodal and face center points are identified by local mesh index.
 * Cell center points are identified by a triplet <region,subregion,index>.
 *
 * The sad reality is, a boundary flux MPFA stencil may be comprised of a mix of
 * cell and face centroids, so we have to discriminate between them at runtime
 */
struct PointDescriptor
{
  /// Enum to classify the variable location
  enum class Tag { CELL, FACE, NODE };

  /// The tag
  Tag tag;

  /// union to characterize a PointDescriptor
  union
  {
    /// node index
    localIndex nodeIndex;
    /// face index
    localIndex faceIndex;
    /// CellDescriptor index
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

  /// Alias for CatalogInterface, necessary declarations for factory instantiation of derived classes
  using CatalogInterface = dataRepository::CatalogInterface< FluxApproximationBase, string const &, Group * const >;
  /**
   * @brief Return the data type in the data repository.
   * @return the data type in the data repository
   */
  static typename CatalogInterface::CatalogType & GetCatalog();

  /// Alias for stored stencil types
  using BoundaryStencil = FluxStencil< PointDescriptor, real64 >;

  FluxApproximationBase() = delete;

  /**
   * @brief Constructor.
   * @param name the name of the FluxApproximationBase in the data repository
   * @param parent the parent group of this group.
   */
  FluxApproximationBase( string const & name, dataRepository::Group * const parent );

  /**
   * @brief Return a boundary stencil by face set name.
   * @param[in] setName the face set name
   * @return the boundary stencil by face set name
   */
  BoundaryStencil const & getBoundaryStencil( string const & setName ) const;

  /**
   * @copydoc getBoundaryStencil( string const & ) const
   */
  BoundaryStencil & getBoundaryStencil( string const & setName );

  /**
   * @brief Check if a stencil exists.
   * @param[in] setName the face set name
   * @return true if a stencil exists
   */
  bool hasBoundaryStencil( string const & setName ) const;

  /**
   * @brief Call a user-provided function for each stencil.
   * @tparam LAMBDA The type of lambda function passed into the parameter list.
   * @param[in] lambda The LAMBDA function
   */
  template< typename LAMBDA >
  void forAllStencils( LAMBDA && lambda ) const;

  /**
   * @brief Call a user-provided function for the each stencil according to the provided TYPE.
   * @tparam TYPE The type to be passed to forWrappers
   * @tparam TYPES Other types to be passed to forWrappers
   * @tparam LAMBDA The type of lambda function passed into the parameter list.
   * @param[in] lambda The LAMBDA function
   */
  template< typename TYPE, typename ... TYPES, typename LAMBDA >
  void forStencils( LAMBDA && lambda ) const;

  /**
   * @brief Call a user-provided function for each boundary stencil.
   * @tparam LAMBDA The type of lambda function passed into the parameter list.
   * @param[in] lambda The LAMBDA function
   */
  template< typename LAMBDA >
  void forBoundaryStencils( LAMBDA && lambda ) const;

  /**
   * @brief Triggers computation of the stencil, implemented in derived classes.
   * @param[in,out] domain The domain on which to perform the stencil computation
   */
  void compute( DomainPartition & domain );

  /**
   * @brief Add a new fracture stencil.
   * @param[in,out] domain The domain on which to add the fracture stencil
   * @param[in] faceElementRegionName the face element region name
   * @param[in] initFlag if true initialize physical fields, like pressure
   */
  virtual void addToFractureStencil( DomainPartition & domain,
                                     string const & faceElementRegionName,
                                     bool const initFlag )
  {
    GEOSX_UNUSED_VAR( domain );
    GEOSX_UNUSED_VAR( faceElementRegionName );
    GEOSX_UNUSED_VAR( initFlag );
  }

  virtual void addEDFracToFractureStencil( DomainPartition & GEOSX_UNUSED_PARAM( domain ),
                                           string const & GEOSX_UNUSED_PARAM( embeddedSurfaceRegionName ) ) {}

  /**
   * @brief View keys.
   */
  struct viewKeyStruct
  {
    /// The key for fieldName
    static constexpr auto fieldNameString             = "fieldName";
    /// The key for boundaryFieldName
    static constexpr auto boundaryFieldNameString     = "boundaryFieldName";
    /// The key for coefficientName
    static constexpr auto coeffNameString             = "coefficientName";
    /// The key for targetRegions
    static constexpr auto targetRegionsString         = "targetRegions";
    /// The key for areaRelTol
    static constexpr auto areaRelativeToleranceString = "areaRelTol";
    /// The key for cellStencil
    static constexpr auto cellStencilString           = "cellStencil";
    /// The key for fractureStencil
    static constexpr auto fractureStencilString       = "fractureStencil";
  };

  struct groupKeyStruct
  {};

  /**
   * @brief Returns the target region name.
   * @return the target region name
   */
  string_array const & targetRegions() const { return m_targetRegions; }
  /**
   * @copydoc targetRegions() const
   */
  string_array & targetRegions()       { return m_targetRegions; }

protected:

  /**
   * @brief Called by InitializePostInitialConditions() prior to initializing sub-Groups.
   * @param[in] rootGroup A group that is passed in to the initialization functions
   *            in order to facilitate the initialization.
   */
  virtual void InitializePostInitialConditions_PreSubGroups( Group * const rootGroup ) override;

  /**
   * @brief Actual computation of the cell-to-cell stencil, to be overridden by implementations.
   * @param[in] domain the domain on which to perform the computation
   */
  virtual void computeCellStencil( DomainPartition const & domain ) = 0;

  /**
   * @brief Actual computation of the boundary stencil, to be overridden by implementations.
   * @param[in] domain the domain on which to perform the computation
   * @param[in] faceSet set of faces
   * @param[out] stencil the boundary stencil
   */
  virtual void computeBoundaryStencil( DomainPartition const & domain,
                                       SortedArrayView< localIndex const > const & faceSet,
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

template< typename LAMBDA >
void FluxApproximationBase::forAllStencils( LAMBDA && lambda ) const
{
//TODO remove dependence on CellElementStencilTPFA and FaceElementStencil
  this->forWrappers< CellElementStencilTPFA, FaceElementStencil >( [&] ( auto const & wrapper )
  {
    lambda( wrapper.reference());
  } );
}

template< typename TYPE, typename ... TYPES, typename LAMBDA >
void FluxApproximationBase::forStencils( LAMBDA && lambda ) const
{
  this->forWrappers< TYPE, TYPES... >( [&] ( auto const & wrapper )
  {
    lambda( wrapper.reference());
  } );
}


template< typename LAMBDA >
void FluxApproximationBase::forBoundaryStencils( LAMBDA && lambda ) const
{
  this->forWrappers< BoundaryStencil >( [&] ( auto const & wrapper )
  {
    lambda( wrapper.reference());
  } );
}

} // namespace geosx

#endif //GEOSX_FINITEVOLUME_FLUXAPPROXIMATIONBASE_HPP_
