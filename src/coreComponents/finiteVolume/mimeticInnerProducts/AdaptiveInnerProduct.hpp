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
 * @file AdaptiveInnerProduct.hpp
 */

#ifndef GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_ADAPTIVEINNERPRODUCT_HPP_
#define GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_ADAPTIVEINNERPRODUCT_HPP_

#include "finiteVolume/mimeticInnerProducts/TPFAInnerProduct.hpp"

namespace geos
{
namespace mimeticInnerProduct
{

/**
 * @class AdaptiveInnerProduct
 * @tparam MFD_IP the consistent (MFD) inner product used in the MFD-compatible cells
 *
 * @brief Adaptive mixed-form inner product blending the MFD mass matrix of @p MFD_IP with the
 *        diagonal TPFA mass matrix according to a cell-wise weight:
 *
 *          M_adapt = eta * M_mfd + (1 - eta) * M_tpfa,   eta in [0, 1].
 *
 * With the binary Global Adaptation stencil flag (eta = chi in {0, 1}) this reproduces the
 * adaptive TPFA/MFD discretization; a fractional weight yields a convex combination whose
 * spectral bounds follow from those of the two operands (uniform coercivity of the blend).
 *
 * The blend is branch-free by design and lives at the inner-product level so that any
 * mixed mimetic assembly kernel (single-phase today, compositional in the future) can adopt
 * the adaptation by wrapping its inner product type: AdaptiveInnerProduct< IP >::computeM.
 */
template< typename MFD_IP >
struct AdaptiveInnerProduct
{
  /**
   * @brief Compute the blended mixed-form mass matrix M_adapt in a given element.
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemVolume the volume of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[in] weight the cell-wise blending weight eta in [0, 1] (1 = MFD, 0 = TPFA)
   * @param[inout] M the output blended mass matrix
   */
  template< localIndex NF >
  GEOS_HOST_DEVICE
  static void
  computeM( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
            ArrayOfArraysView< localIndex const > const & faceToNodes,
            arraySlice1d< localIndex const > const & elemToFaces,
            arraySlice1d< real64 const > const & elemCenter,
            real64 const & elemVolume,
            real64 const (&elemPerm)[ 3 ],
            real64 const & lengthTolerance,
            real64 const & weight,
            arraySlice2d< real64 > const & M )
  {
    stackArray2d< real64, NF *NF > mfdMassMatrix( NF, NF );
    stackArray2d< real64, NF *NF > tpfaMassMatrix( NF, NF );

    MFD_IP::template computeM< NF >( nodePosition,
                                     faceToNodes,
                                     elemToFaces,
                                     elemCenter,
                                     elemVolume,
                                     elemPerm,
                                     lengthTolerance,
                                     mfdMassMatrix );

    TPFAInnerProduct::computeM< NF >( nodePosition,
                                      faceToNodes,
                                      elemToFaces,
                                      elemCenter,
                                      elemVolume,
                                      elemPerm,
                                      lengthTolerance,
                                      tpfaMassMatrix );

    for( localIndex i = 0; i < NF; ++i )
    {
      for( localIndex j = 0; j < NF; ++j )
      {
        M[i][j] = weight * mfdMassMatrix( i, j ) + ( 1.0 - weight ) * tpfaMassMatrix( i, j );
      }
    }
  }
};

} // end namespace mimeticInnerProduct

} // end namespace geos

#endif //GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_ADAPTIVEINNERPRODUCT_HPP_
