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
 * @file MsrsbLevelBuilderBase.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_MSRSBLEVELBUILDERBASE_HPP
#define GEOSX_LINEARALGEBRA_MULTISCALE_MSRSBLEVELBUILDERBASE_HPP

#include "linearAlgebra/common/PreconditionerBase.hpp"
#include "linearAlgebra/multiscale/mesh/DofManager.hpp"
#include "linearAlgebra/multiscale/mesh/MeshLevel.hpp"
#include "linearAlgebra/multiscale/LevelBuilderBase.hpp"

namespace geos
{
namespace multiscale
{

/**
 * @brief Base class for MsRSB level builders
 * @tparam LAI linear algebra interface type
 */
template< typename LAI >
class MsrsbLevelBuilderBase : public LevelBuilderBase< LAI >
{
public:

  /// Alias for base type
  using Base = LevelBuilderBase< LAI >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  /// Alias for operator type
  using Operator = typename Base::Operator;

  /**
   * @brief Constructor.
   * @param name level name
   * @param params linear solver parameters
   */
  explicit MsrsbLevelBuilderBase( string name, LinearSolverParameters params );

  virtual Matrix const & prolongation() const override
  {
    return m_prolongation;
  }

  virtual Operator const & restriction() const override
  {
    return *m_restriction;
  }

  virtual Matrix const & matrix() const override
  {
    return m_matrix;
  }

  virtual PreconditionerBase< LAI > const * presmoother() const override
  {
    return m_preSmoother.get();
  }

  /**
   * @brief @return current level's pre-smoothing operator (mutable access)
   */
  PreconditionerBase< LAI > * presmoother()
  {
    return m_preSmoother.get();
  }

  virtual PreconditionerBase< LAI > const * postsmoother() const override
  {
    return m_postSmoother.get();
  }

  /**
   * @brief @return current level's post-smoothing operator (mutable access)
   */
  PreconditionerBase< LAI > * postsmoother()
  {
    return m_postSmoother.get();
  }

  /**
   * @brief @return the level's multiscale DofManager
   */
  multiscale::DofManager & dofManager() { return m_dofManager; }

  /**
   * @brief @return the level's multiscale DofManager
   */
  multiscale::DofManager const & dofManager() const { return m_dofManager; }

  virtual void compute( Matrix const & fineMatrix ) override;

  /**
   * @brief Update current level's prolongation using a new previous level matrix
   * @param fineMatrix the previous level matrix
   * @return whether an update was actually performed (it may not be when the change is below threshold)
   */
  virtual bool updateProlongation( Matrix const & fineMatrix ) = 0;

protected:

  using Base::m_params;
  using Base::m_name;

  /// Pointer to the fine level
  MsrsbLevelBuilderBase const * m_fineLevel{};

  /// Prolongation matrix P
  Matrix m_prolongation;

  /// Restriction (kept as abstract operator to allow for memory efficiency, e.g. when R = P^T)
  std::unique_ptr< Operator > m_restriction;

  /// Level operator matrix
  Matrix m_matrix;

  /// DofManager for the matrix
  multiscale::DofManager m_dofManager;

  /// Pre-smoothing operator
  std::unique_ptr< PreconditionerBase< LAI > > m_preSmoother;

  /// Post-smoothing operator
  std::unique_ptr< PreconditionerBase< LAI > > m_postSmoother;
};

} // geosx
} // multiscale

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_MSRSBLEVELBUILDERBASE_HPP
