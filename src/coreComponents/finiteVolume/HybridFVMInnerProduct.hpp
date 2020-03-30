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

static constexpr localIndex MAX_NUM_FACES = 15;

struct HybridFVMInnerProductType
{
  static constexpr integer TPFA = 0;
  static constexpr integer QUASI_TPFA = 1;
};

/******************************** Helpers ********************************/

struct HybridFVMInnerProductHelper
{

  // for now, I just copy-pasted this function from TwoPointFluxApproximation
  static
  void makeFullTensor( R1Tensor const & values,
                       stackArray2d< real64, 9 > & result );

};

/******************************** TPFA Kernels ********************************/

struct TPFAFaceInnerProductKernel
{
  /**
   * @brief In a given element, recompute the transmissibility matrix in a face using TPFA
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] transMatrix
   *
   */
  inline static void
  Compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & GEOSX_UNUSED_PARAM( nodePosition ),
           ArrayOfArraysView< localIndex const > const & GEOSX_UNUSED_PARAM( faceToNodes ),
           arraySlice1d< localIndex const > const GEOSX_UNUSED_PARAM( elemToFaces ),
           R1Tensor const & GEOSX_UNUSED_PARAM( elemCenter ),
           R1Tensor const & GEOSX_UNUSED_PARAM( elemPerm ),
           real64 const & GEOSX_UNUSED_PARAM( lengthTolerance ),
           stackArray2d< real64, MAX_NUM_FACES *MAX_NUM_FACES > const & GEOSX_UNUSED_PARAM( transMatrix ) )
  {
    GEOSX_LOG_RANK( "Support for FaceElementSubRegion is not implemented in the Hybrid FVM scheme yet" );
  }
};


struct TPFACellInnerProductKernel
{

  /**
   * @brief In a given element, recompute the transmissibility matrix in a cell using TPFA
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] transMatrix
   *
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

struct QTPFACellInnerProductKernel
{

  /**
   * @brief In a given element, recompute the transmissibility matrix using a consistent inner product
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] tParam parameter used in the transmissibility matrix computations
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] transMatrix
   *
   * When tParam = 2, we obtain a scheme that reduces to TPFA
   * on orthogonal meshes, but remains consistent on non-orthogonal meshes
   *
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
