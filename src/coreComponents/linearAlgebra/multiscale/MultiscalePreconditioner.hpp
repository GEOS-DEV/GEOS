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
 * @file MultiscalePreconditioner.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_MULTISCALEPRECONDITIONER_HPP_
#define GEOSX_LINEARALGEBRA_MULTISCALE_MULTISCALEPRECONDITIONER_HPP_

#include "linearAlgebra/common/PreconditionerBase.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "mesh/MeshLevel.hpp"

#include <memory>

namespace geos
{

namespace multiscale
{
template< typename LAI >
class LevelBuilderBase;
}

/**
 * @brief Multiscale preconditioner for near-elliptic and coupled problems.
 * @tparam LAI linear algebra interface type
 */
template< typename LAI >
class MultiscalePreconditioner : public PreconditionerBase< LAI >
{
public:

  /// Alias for base type
  using Base = PreconditionerBase< LAI >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  /// Alias for operator type
  using Operator = LinearOperator< Vector >;

  /**
   * @brief Constructor.
   * @param params linear solver parameters
   * @param domain physical domain
   */
  MultiscalePreconditioner( LinearSolverParameters params, DomainPartition & domain );

  ~MultiscalePreconditioner();

  virtual void setup( Matrix const & mat ) override;

  /**
   * @brief Apply operator to a vector, <tt>dst = this(src)</tt>.
   * @param src input vector
   * @param dst output vector
   *
   * @warning @p src and @p dst cannot alias the same vector (some implementations may allow this).
   */
  virtual void apply( Vector const & src, Vector & dst ) const override;

  virtual void clear() override;

private:

  void createLevels( Matrix const & mat );

  void computeLevel( integer const level ) const;

  void printLevelInfo() const;

  void logMessage( integer const minLevel, string const & msg ) const;

  struct LevelData
  {
    std::unique_ptr< multiscale::LevelBuilderBase< LAI > > builder;
    Matrix const * matrix; ///< level matrix (operator)
    mutable Vector rhs;    ///< level right-hand side vector
    mutable Vector sol;    ///< level solution
    mutable Vector tmp;    ///< level temporary vector used to hold intermediate solutions
  };

  LinearSolverParameters m_params;

  DomainPartition & m_domain;

  std::vector< LevelData > m_levels;

  std::unique_ptr< PreconditionerBase< LAI > > m_coarse_solver;

  bool m_initialized{ false };
};

} // namespace geos

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_MULTISCALEPRECONDITIONER_HPP_
