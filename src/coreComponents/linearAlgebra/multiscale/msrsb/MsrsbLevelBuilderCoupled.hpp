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
 * @file MsrsbLevelBuilderCoupled.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_MSRSBLEVELBUILDERCOUPLED_HPP
#define GEOSX_LINEARALGEBRA_MULTISCALE_MSRSBLEVELBUILDERCOUPLED_HPP

#include "MsrsbLevelBuilder.hpp"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"

namespace geosx
{

namespace multiscale
{

template< typename LAI >
class MsrsbLevelBuilderCoupled : public MsrsbLevelBuilderBase< LAI >
{

public:

  /// Alias for base type
  using Base = MsrsbLevelBuilderBase< LAI >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  explicit MsrsbLevelBuilderCoupled( string name, LinearSolverParameters params );

  virtual void initializeFineLevel( DomainPartition & domain,
                                    geosx::DofManager const & dofManager,
                                    MPI_Comm const & comm ) override;

  virtual void initializeCoarseLevel( LevelBuilderBase< LAI > & fineLevel,
                                      Matrix const & fineMat ) override;

  virtual bool updateProlongation( Matrix const & fineMatrix ) override;

  virtual std::unique_ptr< PreconditionerBase< LAI > > makeCoarseSolver() const override;

private:

  using Base::m_params;
  using Base::m_name;
  using Base::m_prolongation;
  using Base::m_restriction;
  using Base::m_matrix;
  using Base::m_dofManager;
  using Base::m_preSmoother;
  using Base::m_postSmoother;
  using Base::m_fineLevel;

  void initializeCommon( DomainPartition & domain, MPI_Comm const & comm );

  void createSmoothers();

  void buildProlongationStructure( DofManager const & fineDofManager );

  /// A field description for each sub-block
  std::vector< string > m_fields;

  /// Subproblem selector matrices at the current level
  std::vector< Matrix > m_selectors;

  /// Levels for each sub-problem
  std::vector< std::unique_ptr< MsrsbLevelBuilder< LAI > > > m_builders;

  /// Sub-problem prolongators in extracted (local) form
  std::vector< CRSMatrix< real64, globalIndex > > m_prolongationBlocks;

  /// Temporary storage for combined prolongation operator (stored to avoid recreating the structure)
  CRSMatrix< real64, globalIndex > m_localProlongation;
};

} // namespace multiscale

} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_MSRSBLEVELBUILDERCOUPLED_HPP
