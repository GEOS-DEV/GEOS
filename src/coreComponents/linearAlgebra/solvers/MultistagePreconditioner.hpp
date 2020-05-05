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
 * @file MultistagePreconditioner.hpp
 */
#ifndef GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERMULTISTAGE_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERMULTISTAGE_HPP_

#include "linearAlgebra/solvers/PreconditionerBase.hpp"

#include <vector>
#include <memory>

namespace geosx
{

template< typename LAI >
class MultistagePreconditioner : public PreconditionerBase< LAI >
{
public:

  using Base = PreconditionerBase< LAI >;
  using Vector = typename Base::Vector;
  using Matrix = typename Base::Matrix;

  MultistagePreconditioner();

  virtual ~MultistagePreconditioner() override;

  virtual void compute( Matrix const & mat,
                        DofManager const & dofManager ) override;

  virtual void apply( Vector const & src, Vector & dst ) const override;

  virtual void clear() override;

  /**
   * @brief Add a stage to the preconditioner
   * @param[in] stage owning pointer to the operator to add (transfers ownership)
   */
  void addStage( std::unique_ptr< PreconditionerBase< LAI > > stage, localIndex numSteps );

  /**
   * @brief Get the number of stages
   * @return the number of stages added
   */
  localIndex numStages() const { return m_stages.size(); }

  PreconditionerBase< LAI > & stage( localIndex k ) { return *m_stages[k]; }

  PreconditionerBase< LAI > const & stage( localIndex k ) const { return *m_stages[k]; }

private:

  /// Individual stage preconditioners
  std::vector< std::unique_ptr< PreconditionerBase< LAI > > > m_stages;

  /// Number of applications of each stage
  array1d< localIndex > m_numSteps;

  /// Internal vector used to keep track of current residual
  mutable Vector m_residual;

  /// Internal vector used as a temporary buffer for stage solution
  mutable Vector m_solution;
};

}

#endif //GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERMULTISTAGE_HPP_
