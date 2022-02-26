/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PreconditionerNull.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERNULL_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERNULL_HPP_

#include "linearAlgebra/common/PreconditionerBase.hpp"

namespace geos
{

/**
 * @brief Common interface for identity preconditioning operator
 * @tparam LAI linear algebra interface providing vectors, matrices and solvers
 */
template< typename LAI >
class PreconditionerNull : public PreconditionerBase< LAI >
{
public:

  /// Alias for base type
  using Base = PreconditionerBase< LAI >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  /**
   * @brief Apply operator to a vector.
   *
   * @param src Input vector (src).
   * @param dst Output vector (dst).
   */
  virtual void apply( Vector const & src,
                      Vector & dst ) const override
  {
    GEOS_LAI_ASSERT_EQ( this->numGlobalRows(), dst.globalSize() );
    GEOS_LAI_ASSERT_EQ( this->numGlobalCols(), src.globalSize() );
    dst.zero();
  }
};

} // namespace geos

#endif //GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERNULL_HPP_
