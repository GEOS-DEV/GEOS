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
  virtual void compute( Matrix const & mat )
  {
    m_numGlobalRows = mat.numGlobalRows();
    m_numGlobalCols = mat.numGlobalCols();
  }

  /**
   * @brief Compute the preconditioner from a matrix
   * @param mat the matrix to precondition
   * @param dofManager the Degree-of-Freedom manager associated with matrix
   */
  virtual void compute( Matrix const & mat,
                        DofManager const & GEOSX_UNUSED_PARAM( dofManager ) )
  {
    compute( mat );
  }

  virtual globalIndex numGlobalRows() const override
  {
    return m_numGlobalRows;
  }

  virtual globalIndex numGlobalCols() const override
  {
    return m_numGlobalCols;
  }

private:

  globalIndex m_numGlobalRows;
  globalIndex m_numGlobalCols;
};

}

#endif //GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERBASE_HPP_
