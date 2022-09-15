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
 * @file LevelBuilderBase.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_LEVELBUILDERBASE_HPP_
#define GEOSX_LINEARALGEBRA_MULTISCALE_LEVELBUILDERBASE_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/common/LinearOperator.hpp"
#include "linearAlgebra/common/PreconditionerBase.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "mesh/MeshLevel.hpp"

#include <memory>

namespace geosx
{
namespace multiscale
{

template< typename LAI >
class LevelBuilderBase
{
public:

  /// Alias for vector type
  using Vector = typename LAI::ParallelVector;

  /// Alias for matrix type
  using Matrix = typename LAI::ParallelMatrix;

  /// Alias for operator type
  using Operator = LinearOperator< Vector >;

  static std::unique_ptr< LevelBuilderBase< LAI > >
  create( string name, LinearSolverParameters params );

  virtual ~LevelBuilderBase() = default;

  virtual Operator const & prolongation() const = 0;

  virtual Operator const & restriction() const = 0;

  virtual Matrix const & matrix() const = 0;

  virtual Operator const * presmoother() const = 0;

  virtual Operator const * postsmoother() const = 0;

  virtual void initializeFineLevel( DomainPartition & domain,
                                    geosx::DofManager const & dofManager,
                                    MPI_Comm const & comm ) = 0;

  virtual void initializeCoarseLevel( LevelBuilderBase< LAI > & fine,
                                      Matrix const & fineMat ) = 0;

  virtual void compute( Matrix const & fineMatrix ) = 0;

  virtual std::unique_ptr< PreconditionerBase< LAI > > makeCoarseSolver() const = 0;

protected:

  explicit LevelBuilderBase( string name, LinearSolverParameters params );

  string m_name;

  LinearSolverParameters m_params;
};

} // namespace multiscale
} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_LEVELBUILDERBASE_HPP_
