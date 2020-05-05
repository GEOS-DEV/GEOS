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
 * @file SeparateComponentPreconditioner.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_SOLVERS_SEPARATECOMPONENTPRECONDITIONER_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_SEPARATECOMPONENTPRECONDITIONER_HPP_

#include "PreconditionerBase.hpp"

#include <memory>

namespace geosx
{

template< typename LAI >
class SeparateComponentPreconditioner : public PreconditionerBase< LAI >
{
public:

  using Base = PreconditionerBase< LAI >;
  using Vector = typename Base::Vector;
  using Matrix = typename Base::Matrix;

  SeparateComponentPreconditioner( localIndex const numComp,
                                   std::unique_ptr< PreconditionerBase< LAI > > precond );

  virtual ~SeparateComponentPreconditioner() override;

  virtual void compute( Matrix const & mat, DofManager const & dofManager ) override;

  virtual void apply( Vector const & src, Vector & dst ) const override;

  virtual void clear() override;

  Matrix const & getPrecondMatrix() const
  {
    return m_matSC;
  }

  PreconditionerBase< LAI > const & getNestedPrecond() const
  {
    return *m_precond;
  }

private:

  /// Number of components in the matrix
  localIndex m_numComp;

  /// Separate component approximation of the original matrix
  Matrix m_matSC;

  /// Actual preconditioner
  std::unique_ptr< PreconditionerBase< LAI > > m_precond;
};

}

#endif //GEOSX_LINEARALGEBRA_SOLVERS_SEPARATECOMPONENTPRECONDITIONER_HPP_
