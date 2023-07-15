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

#ifndef GEOS_LINEARALGEBRA_SOLVERS_PRECONDITIONERIDENTITY_HPP_
#define GEOS_LINEARALGEBRA_SOLVERS_PRECONDITIONERIDENTITY_HPP_

#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/common/LinearOperator.hpp"
#include "linearAlgebra/common/PreconditionerBase.hpp"

namespace geos
{

/**
 * @brief Common interface for identity preconditioning operator
 * @tparam LAI linear algebra interface providing vectors, matrices and solvers
 */
template< typename LAI >
class PreconditionerIdentity : public PreconditionerBase< LAI >
{
public:

  /// Alias for base type
  using Base = PreconditionerBase< LAI >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  virtual ~PreconditionerIdentity() = default;

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
    dst.copy( src );
  }
};

template< typename LAI >
class MatrixFreePreconditionerIdentity : public PreconditionerIdentity< LAI >
{
public:
  MatrixFreePreconditionerIdentity( DofManager & dofManager )
    : m_dofManager( dofManager )
  { }

  virtual globalIndex numGlobalRows() const
  {
    return m_dofManager.numGlobalDofs();
  }

  virtual globalIndex numGlobalCols() const
  {
    return m_dofManager.numGlobalDofs();
  }

  virtual localIndex numLocalRows() const
  {
    return m_dofManager.numLocalDofs();
  }

  virtual localIndex numLocalCols() const
  {
    return m_dofManager.numLocalDofs();
  }

  virtual MPI_Comm comm() const
  {
    return MPI_COMM_GEOSX;
  }

private:
  DofManager & m_dofManager;

};

}

#endif //GEOS_LINEARALGEBRA_SOLVERS_PRECONDITIONERIDENTITY_HPP_
