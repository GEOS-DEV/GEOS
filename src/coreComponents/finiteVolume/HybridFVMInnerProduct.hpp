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
 * @file HybridFVMInnerProduct.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FINITEVOLUME_HYBRIDFVMINNERPRODUCT_HPP
#define GEOSX_PHYSICSSOLVERS_FINITEVOLUME_HYBRIDFVMINNERPRODUCT_HPP

#include "common/DataTypes.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

namespace geosx
{

namespace HybridFVMInnerProduct
{

/// Maximum number of faces per cell
static constexpr localIndex MAX_NUM_FACES = 15;

/**
 * @struct HybridFVMInnerProductType
 * @brief Struct describing the possible types of inner product.
 */
struct HybridFVMInnerProductType
{
  static constexpr integer TPFA = 0;       ///< classic TPFA
  static constexpr integer QUASI_TPFA = 1; ///< quasi TPFA
};

/******************************** Helpers ********************************/

/**
 * @struct HybridFVMInnerProductHelper
 * @brief Helper struct handling inner product for hybrid finite volume schemes.
 */
struct HybridFVMInnerProductHelper
{

  /**
   * @brief Create a full tensor from an array.
   * @param[in] values the input array
   * @param[out] result the full tensor
   */
  static
  void makeFullTensor( R1Tensor const & values,
                       stackArray2d< real64, 9 > & result );

};

/******************************** TPFA Kernels ********************************/

/**
 * @struct TPFAFaceInnerProductKernel
 * @brief Struct handling face inner product in the quasi TPFA scheme.
 */
struct TPFAFaceInnerProductKernel
{
  /**
   * @brief In a given element, recompute the transmissibility matrix in a face using TPFA.
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] transMatrix the output
   */
  inline static void
  Compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
           ArrayOfArraysView< localIndex const > const & faceToNodes,
           arraySlice1d< localIndex const > const elemToFaces,
           R1Tensor const & elemCenter,
           R1Tensor const & elemPerm,
           real64 const & lengthTolerance,
           stackArray2d< real64, MAX_NUM_FACES *MAX_NUM_FACES > const & transMatrix )
  {
    GEOSX_UNUSED_VAR( nodePosition );
    GEOSX_UNUSED_VAR( faceToNodes );
    GEOSX_UNUSED_VAR( elemToFaces );
    GEOSX_UNUSED_VAR( elemCenter );
    GEOSX_UNUSED_VAR( elemPerm );
    GEOSX_UNUSED_VAR( lengthTolerance );
    GEOSX_UNUSED_VAR( transMatrix );
    GEOSX_LOG_RANK( "Support for FaceElementSubRegion is not implemented in the Hybrid FVM scheme yet" );
  }
};

/**
 * @struct TPFACellInnerProductKernel
 * @brief Struct handling cell inner product in the TPFA scheme.
 */
struct TPFACellInnerProductKernel
{

  /**
   * @brief In a given element, recompute the transmissibility matrix in a cell using TPFA.
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] transMatrix the output
   */
  static void
  Compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
           ArrayOfArraysView< localIndex const > const & faceToNodes,
           arraySlice1d< localIndex const > const elemToFaces,
           R1Tensor const & elemCenter,
           R1Tensor const & elemPerm,
           real64 const & lengthTolerance,
           stackArray2d< real64, MAX_NUM_FACES *MAX_NUM_FACES > const & transMatrix );

};

/******************************** Quasi TPFA Kernel ********************************/

/**
 * @struct QTPFACellInnerProductKernel
 * @brief Struct handling cell inner product in the quasi TPFA scheme.
 */
struct QTPFACellInnerProductKernel
{

  /**
   * @brief In a given element, recompute the transmissibility matrix using a consistent inner product.
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemVolume the volume of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] tParam parameter used in the transmissibility matrix computations
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[in] orthonormalizeWithSVD flag to choose is use orthonormalizing through SVD
   * @param[inout] transMatrix the output
   *
   * When tParam = 2, we obtain a scheme that reduces to TPFA
   * on orthogonal meshes, but remains consistent on non-orthogonal meshes
   */
  static void
  Compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
           ArrayOfArraysView< localIndex const > const & faceToNodes,
           arraySlice1d< localIndex const > const elemToFaces,
           R1Tensor const & elemCenter,
           real64 const & elemVolume,
           R1Tensor const & elemPerm,
           real64 const & tParam,
           real64 const & lengthTolerance,
           bool const & orthonormalizeWithSVD,
           stackArray2d< real64, MAX_NUM_FACES *MAX_NUM_FACES > const & transMatrix );

};


} // namespace HybridFVMInnerProduct

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FINITEVOLUME_HYBRIDFVMINNERPRODUCT_HPP
