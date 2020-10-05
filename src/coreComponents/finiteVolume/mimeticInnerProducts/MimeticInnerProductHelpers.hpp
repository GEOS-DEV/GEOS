/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MimeticInnerProductHelpers.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FINITEVOLUME_MIMETICINNERPRODUCTS_MIMETICINNERPRODUCTHELPERS_HPP
#define GEOSX_PHYSICSSOLVERS_FINITEVOLUME_MIMETICINNERPRODUCTS_MIMETICINNERPRODUCTHELPERS_HPP

namespace geosx
{
namespace mimeticInnerProduct
{

/**
 * @struct MimeticInnerProductHelpers
 * @brief Helper struct handling inner product for hybrid finite volume schemes.
 */
struct MimeticInnerProductHelpers
{

  /**
   * @brief Create a full tensor from an array.
   * @param[in] values the input array
   * @param[out] result the full tensor
   */
  GEOSX_HOST_DEVICE
  static
  void MakeFullTensor( real64 const (&values)[ 3 ],
                       real64 (& result)[ 3 ][ 3 ] )
  {
    LvArray::tensorOps::fill< 3, 3 >( result, 0.0 );
    result[ 0 ][ 0 ] = values[ 0 ];
    result[ 1 ][ 1 ] = values[ 1 ];
    result[ 2 ][ 2 ] = values[ 2 ];
  }

  /**
   * @brief Orthonormalize a set of three vectors
   * @tparam NF vector space dimensionality
   * @param[in,out] q0 first vector
   * @param[in,out] q1 second vector
   * @param[in,out] q2 third vector
   * @param[out] cellToFaceMat a copy of in/out vectors stacked into a matrix
   */
  template< localIndex NF >
  GEOSX_HOST_DEVICE
  static
  void Orthonormalize( real64 (& q0)[ NF ],
                       real64 (& q1)[ NF ],
                       real64 (& q2)[ NF ],
                       real64 (& cellToFaceMat)[ NF ][ 3 ] )
  {
    // modified Gram-Schmidt algorithm

    // q0
    LvArray::tensorOps::scale< NF >( q0, 1.0/LvArray::tensorOps::l2Norm< NF >( q0 ) );

    // q1
    real64 const q0Dotq1 = LvArray::tensorOps::AiBi< NF >( q0, q1 );
    LvArray::tensorOps::scaledAdd< NF >( q1, q0, -q0Dotq1 );
    LvArray::tensorOps::scale< NF >( q1, 1.0/LvArray::tensorOps::l2Norm< NF >( q1 ) );

    // q2
    real64 const q0Dotq2 = LvArray::tensorOps::AiBi< NF >( q0, q2 );
    LvArray::tensorOps::scaledAdd< NF >( q2, q0, -q0Dotq2 );
    real64 const q1Dotq2 = LvArray::tensorOps::AiBi< NF >( q1, q2 );
    LvArray::tensorOps::scaledAdd< NF >( q2, q1, -q1Dotq2 );
    LvArray::tensorOps::scale< NF >( q2, 1.0/LvArray::tensorOps::l2Norm< NF >( q2 ) );

    for( localIndex i = 0; i < NF; ++i )
    {
      cellToFaceMat[ i ][ 0 ] = q0[ i ];
      cellToFaceMat[ i ][ 1 ] = q1[ i ];
      cellToFaceMat[ i ][ 2 ] = q2[ i ];
    }
  }

};

} // namespace mimeticInnerProduct

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FINITEVOLUME_MIMETICINNERPRODUCTS_MIMETICINNERPRODUCTHELPERS_HPP
