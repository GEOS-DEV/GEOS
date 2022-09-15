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

/**
 * @file SinglePhasePoromechanicsSolver.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSSOLVER_HPP_

#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"

namespace geosx
{

class SinglePhasePoromechanicsSolver : public CoupledSolver< SolidMechanicsLagrangianFEM,
                                                             SinglePhaseBase >
{
public:

  using Base = CoupledSolver< SolidMechanicsLagrangianFEM, SinglePhaseBase >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  enum class SolverType : integer
  {
    SolidMechanics = 0,
    Flow = 1
  };

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "poromechanics"; }

  /**
   * @brief main constructor for SinglePhasePoromechanicsSolver objects
   * @param name the name of this instantiation of SinglePhasePoromechanicsSolver in the repository
   * @param parent the parent group of this instantiation of SinglePhasePoromechanicsSolver
   */
  SinglePhasePoromechanicsSolver( const string & name,
                                  Group * const parent );

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new SinglePhasePoromechanicsSolver object through the object catalog.
   */
  static string catalogName() { return "SinglePhasePoromechanics"; }

  /**
   * @brief accessor for the pointer to the solid mechanics solver
   * @return a pointer to the solid mechanics solver
   */
  SolidMechanicsLagrangianFEM * solidMechanicsSolver() const
  {
    return std::get< toUnderlying( SolverType::SolidMechanics ) >( m_solvers );
  }

  /**
   * @brief accessor for the pointer to the flow solver
   * @return a pointer to the flow solver
   */
  SinglePhaseBase * flowSolver() const
  {
    return std::get< toUnderlying( SolverType::Flow ) >( m_solvers );
  }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void registerDataOnMesh( Group & MeshBodies ) override;

  virtual void setupCoupling( DomainPartition const & domain,
                              DofManager & dofManager ) const override;

  virtual void setupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            ParallelVector & rhs,
                            ParallelVector & solution,
                            bool const setSparsity = true ) override;

  virtual std::unique_ptr< PreconditionerBase< LAInterface > >
  createPreconditioner( DomainPartition & domain ) const override;

  virtual void assembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             int const cycleNumber,
                             DomainPartition & domain ) override;

  virtual void updateState( DomainPartition & domain ) override;

  virtual void solveLinearSystem( DofManager const & dofManager,
                                  ParallelMatrix & matrix,
                                  ParallelVector & rhs,
                                  ParallelVector & solution ) override;

  /**@}*/

protected:

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * porousMaterialNamesString() { return "porousMaterialNames"; }
    constexpr static char const * linearSystemScalingString() { return "linearSystemScaling"; }
  };

  virtual void initializePostInitialConditionsPreSubGroups() override;

  virtual void initializePreSubGroups() override;

private:

  // input flag indicating whether system scaling should be performed
  integer m_systemScaling = 0;

  // A vector to perform row/block-wise linear system scaling
  ParallelVector m_scalingVector;
};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSSOLVER_HPP_ */
