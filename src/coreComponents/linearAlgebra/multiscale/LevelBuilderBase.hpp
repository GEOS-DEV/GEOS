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

namespace geos
{
namespace multiscale
{

/**
 * @brief Base class for level builder implementations
 * @tparam LAI linear algebra interface type
 */
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

  /**
   * @brief Factory interface to create level instances.
   * @param name level name
   * @param params solver parameters
   * @return owning pointer to the level instance
   */
  static std::unique_ptr< LevelBuilderBase< LAI > >
  create( string name, LinearSolverParameters params );

  /**
   * @brief Constructor.
   * @param name level name
   * @param params solver parameters (levels may need access to the full set of them)
   */
  explicit LevelBuilderBase( string name, LinearSolverParameters params );

  /**
   * @brief Destructor.
   */
  virtual ~LevelBuilderBase() = default;

  /**
   * @brief @return current level's prolongation operator
   */
  virtual Operator const & prolongation() const = 0;

  /**
   * @brief @return current level's restriction operator
   */
  virtual Operator const & restriction() const = 0;

  /**
   * @brief @return current level's system matrix
   */
  virtual Matrix const & matrix() const = 0;

  /**
   * @brief @return current level's pre-smoothing operator
   */
  virtual Operator const * presmoother() const = 0;

  /**
   * @brief @return current level's post-smoothing operator
   */
  virtual Operator const * postsmoother() const = 0;

  /**
   * @brief Initialize the finest level (level 0).
   * @param domain the physical domain object
   * @param dofManager the source DofManager
   * @param comm MPI communicator
   */
  virtual void initializeFineLevel( DomainPartition & domain,
                                    geos::DofManager const & dofManager,
                                    MPI_Comm const & comm ) = 0;

  /**
   * @brief Initialize a coarse level (levels 1 and above).
   * @param fineLevel the previous (fine) level
   * @param fineMatrix the previous (fine) level system matrix
   */
  virtual void initializeCoarseLevel( LevelBuilderBase< LAI > & fineLevel,
                                      Matrix const & fineMatrix ) = 0;

  /**
   * @brief Compute the current level.
   * @param fineMatrix the previous (fine) level system matrix
   */
  virtual void compute( Matrix const & fineMatrix ) = 0;

  /**
   * @brief Instantiate coarsest level solver.
   * @return owning pointer to the solver
   */
  virtual std::unique_ptr< PreconditionerBase< LAI > > makeCoarseSolver() const = 0;

protected:

  /// Level name (label)
  string m_name;

  /// Linear solver top-level parameters
  LinearSolverParameters m_params;
};

} // namespace multiscale
} // namespace geos

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_LEVELBUILDERBASE_HPP_
