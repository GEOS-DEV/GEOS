/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FiniteElementBase.hpp
 */


#ifndef GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_FINITEELEMENTBASE_HPP_
#define GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_FINITEELEMENTBASE_HPP_

#include "common/DataTypes.hpp"
#include "common/GeosxMacros.hpp"
#include "finiteElement/PDEUtilities.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/EdgeManager.hpp"
#include "mesh/FaceManager.hpp"

namespace geos
{
namespace finiteElement
{

/// @brief Device-compatible (non virtual) Base class for all finite element formulations.
template< int NUM_SUPPORT_POINTS,
          int NUM_FACES,
          int NUM_QUADRATURE_POINTS >
class FiniteElementBase_impl
{
public:

  /// The number of nodes per element.
  constexpr static localIndex numNodes = NUM_SUPPORT_POINTS;

  /// The number of support points per element.
  constexpr static localIndex numSupportPoints = NUM_SUPPORT_POINTS;

  /// The maximum number of support points per element.
  constexpr static localIndex maxSupportPoints = numSupportPoints;

  /// The number of faces per element.
  constexpr static localIndex numFaces = NUM_FACES;

  /// The number of quadrature points per element.
  constexpr static localIndex numQuadraturePoints = NUM_QUADRATURE_POINTS;

  /// Number of sampling points.
  constexpr static int numSamplingPointsPerDirection = 10;

  /// The number of sampling points per element.
  constexpr static int numSamplingPoints = numSamplingPointsPerDirection * numSamplingPointsPerDirection * numSamplingPointsPerDirection;

#ifdef __CUDACC__
  #pragma diag_push
  #pragma nv_diag_suppress 20012
#endif
  /// Default constructor.
  GEOS_HOST_DEVICE FiniteElementBase_impl() = default;
  /// Default destructor.
  GEOS_HOST_DEVICE ~FiniteElementBase_impl() = default;
  /// Default copy constructor.
  GEOS_HOST_DEVICE FiniteElementBase_impl( FiniteElementBase_impl const & ) = default;
  /**
   * @brief Default copy assignment operator.
   * @return A reference to this object.
   */
  GEOS_HOST_DEVICE FiniteElementBase_impl & operator=( FiniteElementBase_impl const & ) = default;
  /// Default move constructor.
  GEOS_HOST_DEVICE FiniteElementBase_impl( FiniteElementBase_impl && ) = default;
  /**
   * @brief Default move assignment operator.
   * @return A reference to this object.
   */
  GEOS_HOST_DEVICE FiniteElementBase_impl & operator=( FiniteElementBase_impl && ) = default;
#ifdef __CUDACC__
  #pragma diag_pop
#endif

  /**
   * @struct StackVariables
   * @brief Kernel variables allocated on the stack.
   *
   * Contains variables that will be allocated on the stack. Used only by Virtual Element classes to
   * hold the computed projections of basis functions
   */
  struct StackVariables
  {};

  /**
   * @struct MeshData
   * @brief Variables used to initialize the class.
   */
  template< typename SUBREGION_TYPE >
  struct MeshData
  {};

  /**
   * @brief Get the number of quadrature points.
   * @return The number of quadrature points.
   */
  GEOS_HOST_DEVICE
  static localIndex getNumQuadraturePoints()
  {
    return numQuadraturePoints;
  }

  /**
   * @brief Get the number of quadrature points.
   * @param stack Stack variables as filled by @ref setupStack.
   * @return The number of quadrature points.
   */
  template< typename STACK_VARIABLES_TYPE >
  GEOS_HOST_DEVICE
  static localIndex getNumQuadraturePoints( STACK_VARIABLES_TYPE const & stack )
  {
    GEOS_UNUSED_VAR( stack );
    return numQuadraturePoints;
  }

  /**
   * @brief Get the number of support points.
   * @return The number of support points.
   */
  GEOS_HOST_DEVICE
  static localIndex getNumSupportPoints()
  {
    return numNodes;
  }

  /**
   * @brief Get the number of support points.
   * @param stack Object that holds stack variables.
   * @return The number of support points.
   */
  template< typename STACK_VARIABLES_TYPE >
  GEOS_HOST_DEVICE
  static localIndex getNumSupportPoints( STACK_VARIABLES_TYPE const & stack )
  {
    GEOS_UNUSED_VAR( stack );
    return numNodes;
  }

  /**
   * @brief Get the maximum number of support points.
   * @return The maximum number of support points.
   */
  GEOS_HOST_DEVICE
  static localIndex getMaxSupportPoints()
  {
    return maxSupportPoints;
  }

  /**
   * @brief Get the Sampling Point Coord In the Parent Space
   *
   * @param linearIndex linear index of the sampling point
   * @param samplingPointCoord coordinates of the sampling point
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void getSamplingPointCoordInParentSpace( int const & linearIndex,
                                                  real64 (& samplingPointCoord)[3] )
  {
    GEOS_UNUSED_VAR( linearIndex, samplingPointCoord );
    GEOS_ERROR( " Element type not supported." );
  }



  /**
   * @brief Method to fill a MeshData object.
   * @param nodeManager The node manager.
   * @param edgeManager The edge manager.
   * @param faceManager The face manager.
   * @param cellSubRegion The cell sub-region for which the element has to be initialized.
   * @param meshData MeshData struct to be filled.
   */
  template< typename SUBREGION_TYPE,
            typename MESH_DATA_TYPE >
  static void fillMeshData( NodeManager const & nodeManager,
                            EdgeManager const & edgeManager,
                            FaceManager const & faceManager,
                            SUBREGION_TYPE const & cellSubRegion,
                            MESH_DATA_TYPE & meshData )
  {
    GEOS_UNUSED_VAR( nodeManager,
                     edgeManager,
                     faceManager,
                     cellSubRegion,
                     meshData );
  }

  /**
   * @brief Abstract initialization method.
   * @details It calls the fillMeshData method of the specific element implementation.
   * @tparam LEAF Type of the derived finite element implementation.
   * @param nodeManager The node manager.
   * @param edgeManager The edge manager.
   * @param faceManager The face manager.
   * @param cellSubRegion The cell sub-region for which the element has to be initialized.
   * @param meshData The struct to be filled according to the @p LEAF class needs.
   */
  template< typename LEAF,
            typename SUBREGION_TYPE,
            typename MESH_DATA_TYPE >
  static void initialize( NodeManager const & nodeManager,
                          EdgeManager const & edgeManager,
                          FaceManager const & faceManager,
                          SUBREGION_TYPE const & cellSubRegion,
                          MESH_DATA_TYPE & meshData
                          )
  {
    LEAF::template fillMeshData< SUBREGION_TYPE >( nodeManager,
                                                   edgeManager,
                                                   faceManager,
                                                   cellSubRegion,
                                                   meshData );
  }


  /**
   * @brief Empty setup method.
   * @param cellIndex The index of the cell with respect to the cell sub region.
   * @param meshData MeshData struct filled by @ref fillMeshData.
   * @param stack Object that holds stack variables.
   */
  template< typename MESH_DATA_TYPE,
            typename STACK_VARIABLES_TYPE >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void setupStack( localIndex const & cellIndex,
                          MESH_DATA_TYPE const & meshData,
                          STACK_VARIABLES_TYPE & stack )
  {
    GEOS_UNUSED_VAR( cellIndex,
                     meshData,
                     stack );

  }


  /**
   * @brief Abstract setup method, possibly computing cell-dependent properties.
   * @tparam LEAF Type of the derived finite element implementation.
   * @param cellIndex The index of the cell with respect to the cell sub region to which the element
   * has been initialized previously (see @ref initialize).
   * @param meshData A MeshData object previously filled.
   * @param stack Object that holds stack variables.
   */
  template< typename LEAF,
            typename MESH_DATA_TYPE >
  GEOS_HOST_DEVICE
  void setup( localIndex const & cellIndex,
              MESH_DATA_TYPE const & meshData,
              typename LEAF::StackVariables & stack ) const
  {
    LEAF::setupStack( cellIndex, meshData, stack );
  }


  /**
   * @brief Getter for the function space.
   * @tparam The number of components per support point (i.e., 1 if
   *   scalar variable, 3 if vector variable)
   * @return The function space.
   */
  template< int N >
  GEOS_HOST_DEVICE
  constexpr static PDEUtilities::FunctionSpace getFunctionSpace()
  {
    if constexpr ( N == 1 )
    {
      return PDEUtilities::FunctionSpace::H1;
    }
    else if constexpr ( N == 3 )
    {
      return PDEUtilities::FunctionSpace::H1vector;
    }
    else
    {
      static_assert( N == 1 || N == 3, "Unsupported number of components per support point" );
    }
  }



  /**
   * @brief Empty method, here for compatibility with methods that require a stabilization of the
   * grad-grad bilinear form.
   * @tparam NUMDOFSPERTRIALSUPPORTPOINT Number of degrees of freedom for each support point.
   * @tparam MAXSUPPORTPOINTS Maximum number of support points allowed for this element.
   * @tparam UPPER If true only the upper triangular part of @p matrix is modified.
   * @param stack Stack variables as filled by @ref setupStack.
   * @param matrix The matrix that needs to be stabilized.
   * @param scaleFactor Scaling of the stabilization matrix.
   */
  template< localIndex NUMDOFSPERTRIALSUPPORTPOINT,
            localIndex MAXSUPPORTPOINTS,
            bool UPPER,
            typename STACK_VARIABLES_TYPE >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void addGradGradStabilization( STACK_VARIABLES_TYPE const & stack,
                                        real64 ( & matrix )[MAXSUPPORTPOINTS * NUMDOFSPERTRIALSUPPORTPOINT][MAXSUPPORTPOINTS * NUMDOFSPERTRIALSUPPORTPOINT],
                                        real64 const & scaleFactor )
  {
    GEOS_UNUSED_VAR( stack,
                     matrix,
                     scaleFactor );
  }



  /**
   * @brief Add stabilization of grad-grad bilinear form to input matrix.
   * @tparam LEAF Type of the derived finite element implementation.
   * @tparam NUMDOFSPERTRIALSUPPORTPOINT Number of degrees of freedom for each support point.
   * @tparam UPPER If true only the upper triangular part of @p matrix is modified.
   * @param stack Stack variables created by a call to @ref setup.
   * @param matrix The input matrix to which values have to be added.
   * @param scaleFactor Optional scaling of the stabilization matrix. Defaults to 1.0.
   */
  template< typename LEAF, localIndex NUMDOFSPERTRIALSUPPORTPOINT, bool UPPER = false >
  GEOS_HOST_DEVICE
  void addGradGradStabilizationMatrix( typename LEAF::StackVariables const & stack,
                                       real64 ( & matrix )[LEAF::maxSupportPoints * NUMDOFSPERTRIALSUPPORTPOINT][LEAF::maxSupportPoints * NUMDOFSPERTRIALSUPPORTPOINT],
                                       real64 const scaleFactor = 1.0 ) const
  {
    LEAF::template addGradGradStabilization< NUMDOFSPERTRIALSUPPORTPOINT,
                                             LEAF::maxSupportPoints,
                                             UPPER >( stack,
                                                      matrix,
                                                      scaleFactor );
  }

  /**
   * @brief Empty method, here for compatibility with methods that require a stabilization of the
   * grad-grad bilinear form.
   * @details This method is intended to be used with @p targetVector being the residual and @p dofs
   * being the degrees of freedom of the previous solution.
   * @tparam NUMDOFSPERTRIALSUPPORTPOINT Number of degrees of freedom for each support point.
   * @param stack Stack variables as filled by @ref setupStack.
   * @param dofs The degrees of freedom of the function where the stabilization operator has to be
   * evaluated.
   * @param targetVector The input vector to which values have to be added, seen in chunks of length
   * @p NUMDOFSPERTRIALSUPPORTPOINT.
   * @param scaleFactor Scaling of the stabilization matrix.
   */
  template< localIndex NUMDOFSPERTRIALSUPPORTPOINT,
            localIndex MAXSUPPORTPOINTS,
            typename STACK_VARIABLES_TYPE >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void addEvaluatedGradGradStabilization( STACK_VARIABLES_TYPE const & stack,
                                                 real64 const ( &dofs )[MAXSUPPORTPOINTS][NUMDOFSPERTRIALSUPPORTPOINT],
                                                 real64 ( & targetVector )[MAXSUPPORTPOINTS][NUMDOFSPERTRIALSUPPORTPOINT],
                                                 real64 const scaleFactor )
  {
    GEOS_UNUSED_VAR( stack );
    GEOS_UNUSED_VAR( dofs );
    GEOS_UNUSED_VAR( targetVector );
    GEOS_UNUSED_VAR( scaleFactor );
  }

  /**
   * @brief Add a grad-grad stabilization operator evaluated at a provided vector of dofs to input
   * vector.
   * @details This method is used to modify a residual consistently when the jacobian includes a
   * stabilization term.
   * @tparam LEAF Type of the derived finite element implementation.
   * @tparam NUMDOFSPERTRIALSUPPORTPOINT Number of degrees of freedom for each support point.
   * @param stack Stack variables created by a call to @ref setup.
   * @param dofs The vector of dofs to evaluate the stabilization.
   * @param targetVector The input vector to which values have to be added, seen in chunks of length
   * @p NUMDOFSPERTRIALSUPPORTPOINT.
   * @param scaleFactor Optional scaling of the stabilization matrix. Defaults to 1.0.
   */
  template< typename LEAF, localIndex NUMDOFSPERTRIALSUPPORTPOINT >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void
  addEvaluatedGradGradStabilizationVector( typename LEAF::StackVariables const & stack,
                                           real64 const ( &dofs )[LEAF::maxSupportPoints]
                                           [NUMDOFSPERTRIALSUPPORTPOINT],
                                           real64 ( & targetVector )[LEAF::maxSupportPoints]
                                           [NUMDOFSPERTRIALSUPPORTPOINT],
                                           real64 const scaleFactor = 1.0 ) const
  {
    LEAF::template
    addEvaluatedGradGradStabilization< NUMDOFSPERTRIALSUPPORTPOINT, LEAF::maxSupportPoints >( stack,
                                                                                              dofs,
                                                                                              targetVector,
                                                                                              scaleFactor );
  }


};


/**
 * @brief Base class for FEM element implementations.
 */
class FiniteElementBase
{
public:

  /**
   * @brief Default constructor.
   * @param numSupportPoints The number of support points.
   * @param maxSupportPoints The maximum number of support points.
   * @param numQuadraturePoints The number of quadrature points.
   */
  FiniteElementBase( localIndex const numSupportPoints,
                     localIndex const maxSupportPoints,
                     localIndex const numQuadraturePoints ):
    m_numSupportPoints( numSupportPoints ),
    m_maxSupportPoints( maxSupportPoints ),
    m_numQuadraturePoints( numQuadraturePoints )
  { }

  /**
   * @brief Destructor
   */
  GEOS_HOST_DEVICE
  virtual ~FiniteElementBase() = default;

  /**
   * @brief Getter for the number of quadrature points per element.
   * @return The number of quadrature points per element.
   */
  localIndex getNumQuadraturePoints() const
  {
    return m_numQuadraturePoints;
  }

  /**
   * @brief Getter for the number of support points per element.
   * @return The number of support points per element.
   */
  localIndex getNumSupportPoints() const { return m_numSupportPoints; };

  /**
   * @brief Get the maximum number of support points for this element.
   * @details This should be used to know the size of pre-allocated objects whose size depend on the
   * number of support points.
   * @return The number of maximum support points for this element.
   */
  localIndex getMaxSupportPoints() const { return m_maxSupportPoints; };

private:
  localIndex const m_numSupportPoints;
  localIndex const m_maxSupportPoints;
  localIndex const m_numQuadraturePoints;
};

} // namespace geos::finiteElement
} // namespace geos

#endif //GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_FINITEELEMENTBASE_HPP_
