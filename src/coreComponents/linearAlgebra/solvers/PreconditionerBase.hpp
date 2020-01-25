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

#ifndef GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERBASE_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERBASE_HPP_

#include "linearAlgebra/interfaces/LinearOperator.hpp"

namespace geosx
{

class DofManager;

/**
 * @brief Common interface for preconditioning operators
 * @tparam LAI linear algebra interface providing vectors, matrices and solvers
 */
template< typename LAI >
class PreconditionerBase : public LinearOperator< typename LAI::ParallelVector >
{
public:

  virtual ~PreconditionerBase() = default;

  using Vector = typename LinearOperator< typename LAI::ParallelVector >::Vector;
  using Matrix = typename LAI::ParallelMatrix;

  /**
   * @brief Compute the preconditioner from a matrix
   * @param mat the matrix to precondition
   */
  virtual void compute( Matrix const & mat,
                        DofManager const & dofManager ) = 0;

};

}

#endif //GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERBASE_HPP_
