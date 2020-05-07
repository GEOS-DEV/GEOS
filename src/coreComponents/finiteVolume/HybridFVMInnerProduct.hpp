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

/******************************** Helpers ********************************/

struct HybridFVMInnerProductHelper
{

  GEOSX_HOST_DEVICE
  static
  void MakeFullTensor( R1Tensor const & values,
                       arraySlice2d< real64 > const & result );

  template< localIndex NF >
  GEOSX_HOST_DEVICE
  static
  void Orthonormalize( arraySlice1d< real64 > const & q0,
                       arraySlice1d< real64 > const & q1,
                       arraySlice1d< real64 > const & q2,
                       arraySlice2d< real64 > const & cellToFaceMat );


};

/******************************** TPFA Kernel ********************************/

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
  template< localIndex NF >
  static void
  Compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
           ArrayOfArraysView< localIndex const > const & faceToNodes,
           arraySlice1d< localIndex const > const elemToFaces,
           R1Tensor const & elemCenter,
           R1Tensor const & elemPerm,
           real64 const & lengthTolerance,
           arraySlice2d< real64 > const & transMatrix );

};

#define INST_TPFACellInnerProduct( NF ) \
  extern template \
  void TPFACellInnerProductKernel::Compute< NF >( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition, \
                                                  ArrayOfArraysView< localIndex const > const & faceToNodes, \
                                                  arraySlice1d< localIndex const > const elemToFaces, \
                                                  R1Tensor const & elemCenter, \
                                                  R1Tensor const & elemPerm, \
                                                  real64 const & lengthTolerance, \
                                                  arraySlice2d< real64 > const & transMatrix )

INST_TPFACellInnerProduct( 4 );
INST_TPFACellInnerProduct( 5 );
INST_TPFACellInnerProduct( 6 );

#undef INST_TPFACellInnerProduct


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
  template< localIndex NF >
  static void
  Compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
           ArrayOfArraysView< localIndex const > const & faceToNodes,
           arraySlice1d< localIndex const > const elemToFaces,
           R1Tensor const & elemCenter,
           real64 const & elemVolume,
           R1Tensor const & elemPerm,
           real64 const & tParam,
           real64 const & lengthTolerance,
           arraySlice2d< real64 > const & transMatrix );

};

#define INST_QTPFACellInnerProduct( NF ) \
  extern template \
  void QTPFACellInnerProductKernel::Compute< NF >( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition, \
                                                   ArrayOfArraysView< localIndex const > const & faceToNodes, \
                                                   arraySlice1d< localIndex const > const elemToFaces, \
                                                   R1Tensor const & elemCenter, \
                                                   real64 const & elemVolume, \
                                                   R1Tensor const & elemPerm, \
                                                   real64 const & tParam, \
                                                   real64 const & lengthTolerance, \
                                                   arraySlice2d< real64 > const & transMatrix )

INST_QTPFACellInnerProduct( 4 );
INST_QTPFACellInnerProduct( 5 );
INST_QTPFACellInnerProduct( 6 );

#undef INST_QTPFACellInnerProduct

} // namespace HybridFVMInnerProduct

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FINITEVOLUME_HYBRIDFVMINNERPRODUCT_HPP
