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

#ifndef GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERIDENTITY_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERIDENTITY_HPP_

#include "linearAlgebra/interfaces/LinearOperator.hpp"
#include "linearAlgebra/solvers/PreconditionerBase.hpp"

namespace geosx
{

/**
 * @brief Common interface for identity preconditioning operator
 * @tparam LAI linear algebra interface providing vectors, matrices and solvers
 */
template< typename LAI >
class PreconditionerIdentity : public PreconditionerBase< LAI >
{
public:

  virtual ~PreconditionerIdentity() = default;

  using Vector = typename LinearOperator< typename LAI::ParallelVector >::Vector;
  using Matrix = typename LAI::ParallelMatrix;

  /**
   * @brief Compute the preconditioner from a matrix
   * @param mat the matrix to precondition
   */
  virtual void apply( Vector const & src,
                      Vector & dst ) const override
  {
    GEOSX_LAI_ASSERT_EQ( this->numGlobalRows(), dst.globalSize() );
    GEOSX_LAI_ASSERT_EQ( this->numGlobalCols(), src.globalSize() );
    dst = src;
  }
};

}

#endif //GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERIDENTITY_HPP_
