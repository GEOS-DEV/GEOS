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

#if defined(GEOS_USE_DEVICE)
#define CALC_FEM_SHAPE_IN_KERNEL
#endif



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


/**
 * @brief Base class for FEM element implementations.
 */
class FiniteElementBase
{
public:

  /// Default Constructor
  FiniteElementBase() = default;


  /// Number of sampling points.
  constexpr static int numSamplingPointsPerDirection = 10;

  /**
   * @brief Copy Constructor
   * @param source The object to copy.
   */
  GEOS_HOST_DEVICE
  FiniteElementBase( FiniteElementBase const & source ):
#ifdef CALC_FEM_SHAPE_IN_KERNEL
    m_viewGradN(),
    m_viewDetJ()
#else
    m_viewGradN( source.m_viewGradN ),
    m_viewDetJ( source.m_viewDetJ )
#endif
  {
    GEOS_UNUSED_VAR( source ); // suppress warning when CALC_FEM_SHAPE_IN_KERNEL is defined
  }

  /// Default Move constructor
  FiniteElementBase( FiniteElementBase && ) = default;

  /**
   * @brief Deleted copy assignment operator
   * @return deleted
   */
  FiniteElementBase & operator=( FiniteElementBase const & ) = delete;

  /**
   * @brief Deleted move assignment operator
   * @return deleted
   */
  FiniteElementBase & operator=( FiniteElementBase && ) = delete;

  /**
   * @brief Destructor
   */
  GEOS_HOST_DEVICE
  virtual ~FiniteElementBase()
  {}

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
   * @brief Method to fill a MeshData object.
   * @param nodeManager The node manager.
   * @param edgeManager The edge manager.
   * @param faceManager The face manager.
   * @param cellSubRegion The cell sub-region for which the element has to be initialized.
   * @param meshData MeshData struct to be filled.
   */
  template< typename SUBREGION_TYPE >
  static void fillMeshData( NodeManager const & nodeManager,
                            EdgeManager const & edgeManager,
                            FaceManager const & faceManager,
                            SUBREGION_TYPE const & cellSubRegion,
                            MeshData< SUBREGION_TYPE > & meshData )
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
  template< typename LEAF, typename SUBREGION_TYPE >
  static void initialize( NodeManager const & nodeManager,
                          EdgeManager const & edgeManager,
                          FaceManager const & faceManager,
                          SUBREGION_TYPE const & cellSubRegion,
                          typename LEAF::template MeshData< SUBREGION_TYPE > & meshData
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
  template< typename SUBREGION_TYPE >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void setupStack( localIndex const & cellIndex,
                          MeshData< SUBREGION_TYPE > const & meshData,
                          StackVariables & stack )
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
  template< typename LEAF, typename SUBREGION_TYPE >
  GEOS_HOST_DEVICE
  void setup( localIndex const & cellIndex,
              typename LEAF::template MeshData< SUBREGION_TYPE > const & meshData,
              typename LEAF::StackVariables & stack ) const
  {
    LEAF::setupStack( cellIndex, meshData, stack );
  }

  /**
   * @brief Virtual getter for the number of quadrature points per element.
   * @return The number of quadrature points per element.
   */
  GEOS_HOST_DEVICE
  virtual localIndex getNumQuadraturePoints() const = 0;

  /**
   * @brief Virtual getter for the number of support points per element.
   * @return The number of support points per element.
   */
  GEOS_HOST_DEVICE
  virtual localIndex getNumSupportPoints() const = 0;

  /**
   * @brief An helper struct to determine the function space.
   * @tparam N The number of components per support point (i.e., 1 if
   *   scalar variable, 3 if vector variable)
   */
  template< int N >
  struct FunctionSpaceHelper
  {};

  /**
   * @brief Getter for the function space.
   * @tparam The number of components per support point (i.e., 1 if
   *   scalar variable, 3 if vector variable)
   * @return The function space.
   */
  template< int N >
  GEOS_HOST_DEVICE
  constexpr static PDEUtilities::FunctionSpace getFunctionSpace();

  /**
   * @brief Get the maximum number of support points for this element.
   * @details This should be used to know the size of pre-allocated objects whose size depend on the
   * number of support points.
   * @return The number of maximum support points for this element.
   */
  GEOS_HOST_DEVICE
  virtual localIndex getMaxSupportPoints() const = 0;

  /**
   * @brief Get the shape function gradients.
   * @tparam LEAF Type of the derived finite element implementation.
   * @param k The element index.
   * @param q The quadrature point index.
   * @param X Array of coordinates as the reference for the gradients.
   * @param gradN Return array of the shape function gradients.
   * @return The determinant of the Jacobian transformation matrix.
   *
   * This function calls the function to calculate shape function gradients.
   */
  template< typename LEAF >
  GEOS_HOST_DEVICE
  real64 getGradN( localIndex const k,
                   localIndex const q,
                   real64 const (&X)[LEAF::maxSupportPoints][3],
                   real64 ( &gradN )[LEAF::maxSupportPoints][3] ) const;

  /**
   * @brief Get the shape function gradients.
   * @tparam LEAF Type of the derived finite element implementation.
   * @param k The element index.
   * @param q The quadrature point index.
   * @param X Array of coordinates as the reference for the gradients.
   * @param stack Stack variables relative to the element @p k created by a call to @ref setup.
   * @param gradN Return array of the shape function gradients.
   * @return The determinant of the Jacobian transformation matrix.
   *
   * This function calls the function to calculate shape function gradients.
   */
  template< typename LEAF >
  GEOS_HOST_DEVICE
  real64 getGradN( localIndex const k,
                   localIndex const q,
                   real64 const (&X)[LEAF::maxSupportPoints][3],
                   typename LEAF::StackVariables const & stack,
                   real64 ( &gradN )[LEAF::maxSupportPoints][3] ) const;

  /**
   * @brief Get the shape function gradients.
   * @tparam LEAF Type of the derived finite element implementation.
   * @param k The element index.
   * @param q The quadrature point index.
   * @param X dummy variable.
   * @param gradN Return array of the shape function gradients.
   * @return The determinant of the Jacobian transformation matrix.
   *
   * This function returns pre-calculated shape function gradients.
   */
  template< typename LEAF >
  GEOS_HOST_DEVICE
  real64 getGradN( localIndex const k,
                   localIndex const q,
                   int const X,
                   real64 ( &gradN )[LEAF::maxSupportPoints][3] ) const;
  /**
   * @brief Get the shape function gradients.
   * @tparam LEAF Type of the derived finite element implementation.
   * @param k The element index.
   * @param q The quadrature point index.
   * @param X dummy variable.
   * @param stack Stack variables relative to the element @p k created by a call to @ref setup.
   * @param gradN Return array of the shape function gradients.
   * @return The determinant of the Jacobian transformation matrix.
   *
   * This function returns pre-calculated shape function gradients.
   */
  template< typename LEAF >
  GEOS_HOST_DEVICE
  real64 getGradN( localIndex const k,
                   localIndex const q,
                   int const X,
                   typename LEAF::StackVariables const & stack,
                   real64 ( &gradN )[LEAF::maxSupportPoints][3] ) const;


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
  template< localIndex NUMDOFSPERTRIALSUPPORTPOINT, localIndex MAXSUPPORTPOINTS, bool UPPER >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void addGradGradStabilization( StackVariables const & stack,
                                        real64 ( & matrix )
                                        [MAXSUPPORTPOINTS * NUMDOFSPERTRIALSUPPORTPOINT]
                                        [MAXSUPPORTPOINTS * NUMDOFSPERTRIALSUPPORTPOINT],
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
                                       real64 ( & matrix )
                                       [LEAF::maxSupportPoints * NUMDOFSPERTRIALSUPPORTPOINT]
                                       [LEAF::maxSupportPoints * NUMDOFSPERTRIALSUPPORTPOINT],
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
  template< localIndex NUMDOFSPERTRIALSUPPORTPOINT, localIndex MAXSUPPORTPOINTS >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void addEvaluatedGradGradStabilization( StackVariables const & stack,
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






  /**
   * @brief Sets m_viewGradN equal to an input view.
   * @param source The view to assign to m_viewGradN.
   */
  void setGradNView( arrayView4d< real64 const > const & source )
  {
    GEOS_ERROR_IF_NE_MSG( source.size( 1 ),
                          getNumQuadraturePoints(),
                          "2nd-dimension of gradN array does not match number of quadrature points" );
    GEOS_ERROR_IF_NE_MSG( source.size( 2 ),
                          getMaxSupportPoints(),
                          "3rd-dimension of gradN array does not match number of support points" );
    GEOS_ERROR_IF_NE_MSG( source.size( 3 ),
                          3,
                          "4th-dimension of gradN array does not match 3" );

    m_viewGradN = source;
  }

  /**
   * @brief Sets m_viewDetJ equal to an input view.
   * @param source The view to assign to m_viewDetJ.
   */
  void setDetJView( arrayView2d< real64 const > const & source )
  {
    GEOS_ERROR_IF_NE_MSG( source.size( 1 ),
                          getNumQuadraturePoints(),
                          "2nd-dimension of gradN array does not match number of quadrature points" );
    m_viewDetJ = source;
  }

  /**
   * @brief Getter for m_viewGradN
   * @return A new arrayView copy of m_viewGradN.
   */
  arrayView4d< real64 const > getGradNView() const
  {
    return m_viewGradN;
  }

  /**
   * @brief Getter for m_viewDetJ
   * @return A new arrayView copy of m_viewDetJ.
   */
  arrayView2d< real64 const > getDetJView() const
  {
    return m_viewDetJ;
  }


protected:
  /// View to potentially hold pre-calculated shape function gradients.
  arrayView4d< real64 const > m_viewGradN;

  /// View to potentially hold pre-calculated weighted jacobian transformation
  /// determinants.
  arrayView2d< real64 const > m_viewDetJ;
};

/// @cond Doxygen_Suppress

//*************************************************************************************************
//***** Definitions *******************************************************************************
//*************************************************************************************************

template<>
struct FiniteElementBase::FunctionSpaceHelper< 1 >
{
  GEOS_HOST_DEVICE
  constexpr static PDEUtilities::FunctionSpace getFunctionSpace()
  {
    return PDEUtilities::FunctionSpace::H1;
  }
};

template<>
struct FiniteElementBase::FunctionSpaceHelper< 3 >
{
  GEOS_HOST_DEVICE
  constexpr static PDEUtilities::FunctionSpace getFunctionSpace()
  {
    return PDEUtilities::FunctionSpace::H1vector;
  }
};

template< int N >
GEOS_HOST_DEVICE
constexpr PDEUtilities::FunctionSpace FiniteElementBase::getFunctionSpace()
{
  return FunctionSpaceHelper< N >::getFunctionSpace();
}

template< typename LEAF >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 FiniteElementBase::getGradN( localIndex const k,
                                    localIndex const q,
                                    real64 const (&X)[LEAF::maxSupportPoints][3],
                                    real64 (& gradN)[LEAF::maxSupportPoints][3] ) const
{
  GEOS_UNUSED_VAR( k );
  return LEAF::calcGradN( q, X, gradN );
}

template< typename LEAF >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 FiniteElementBase::getGradN( localIndex const k,
                                    localIndex const q,
                                    real64 const (&X)[LEAF::maxSupportPoints][3],
                                    typename LEAF::StackVariables const & stack,
                                    real64 ( & gradN )[LEAF::maxSupportPoints][3] ) const
{
  GEOS_UNUSED_VAR( k );
  return LEAF::calcGradN( q, X, stack, gradN );
}

template< typename LEAF >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 FiniteElementBase::getGradN( localIndex const k,
                                    localIndex const q,
                                    int const X,
                                    real64 (& gradN)[LEAF::maxSupportPoints][3] ) const
{
  GEOS_UNUSED_VAR( X );

  LvArray::tensorOps::copy< LEAF::maxSupportPoints, 3 >( gradN, m_viewGradN[ k ][ q ] );

  return m_viewDetJ( k, q );
}

template< typename LEAF >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 FiniteElementBase::getGradN( localIndex const k,
                                    localIndex const q,
                                    int const X,
                                    typename LEAF::StackVariables const & stack,
                                    real64 (& gradN)[LEAF::maxSupportPoints][3] ) const
{
  GEOS_UNUSED_VAR( X );
  GEOS_UNUSED_VAR( stack );

  LvArray::tensorOps::copy< LEAF::maxSupportPoints, 3 >( gradN, m_viewGradN[ k ][ q ] );

  return m_viewDetJ( k, q );
}
/// @endcond


} // namespace geos::finiteElement
} // namespace geos

#endif //GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_FINITEELEMENTBASE_HPP_
