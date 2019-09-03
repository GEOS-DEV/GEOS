/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file ReservoirSolver.hpp
 *
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_COUPLEDSOLVERS_RESERVOIRSOLVER_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_COUPLEDSOLVERS_RESERVOIRSOLVER_HPP_

#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{

class FlowSolverBase;
class WellSolverBase;

class ReservoirSolver : public SolverBase
{
public:

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  ReservoirSolver( const std::string& name,
                   Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~ReservoirSolver() override;

  /// deleted copy constructor
  ReservoirSolver( ReservoirSolver const & ) = delete;

  /// default move constructor
  ReservoirSolver( ReservoirSolver && ) = default;

  /// deleted assignment operator
  ReservoirSolver & operator=( ReservoirSolver const & ) = delete;

  /// deleted move operator
  ReservoirSolver & operator=( ReservoirSolver && ) = delete;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return "Reservoir"; }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/
  
  virtual void ImplicitStepSetup( real64 const & time_n,
                                  real64 const & dt,
                                  DomainPartition * const domain,
                                  DofManager & dofManager,
                                  ParallelMatrix & matrix,
                                  ParallelVector & rhs,
                                  ParallelVector & solution ) override final;

  virtual void SetupSystem( DomainPartition * const domain,
                            DofManager & dofManager,
                            ParallelMatrix & matrix,
                            ParallelVector & rhs,
                            ParallelVector & solution ) override;

  virtual void SetupDofs( DomainPartition const * const domain,
                          DofManager & dofManager ) const override;

  virtual void AssembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition * const domain,
                               DofManager const & dofManager,
                               ParallelMatrix & matrix,
                               ParallelVector & rhs ) override;

  virtual void ApplyBoundaryConditions( real64 const time,
                                        real64 const dt,
                                        DomainPartition * const domain,
                                        DofManager const & dofManager,
                                        ParallelMatrix & matrix,
                                        ParallelVector & rhs ) override;

  virtual real64
  CalculateResidualNorm( DomainPartition const * const domain,
                         DofManager const & dofManager,
                         ParallelVector const & rhs ) override;

  virtual void SolveSystem( DofManager const & dofManager,
                            ParallelMatrix & matrix,
                            ParallelVector & rhs,
                            ParallelVector & solution ) override;

  virtual bool
  CheckSystemSolution( DomainPartition const * const domain,
                       DofManager const & dofManager,
                       ParallelVector const & solution,
                       real64 const scalingFactor ) override;

  virtual void
  ApplySystemSolution( DofManager const & dofManager,
                       ParallelVector const & solution,
                       real64 const scalingFactor,
                       DomainPartition * const domain ) override;

  virtual void ResetStateToBeginningOfStep( DomainPartition * const domain ) override;

  
  virtual void ImplicitStepComplete( real64 const & time,
                                     real64 const & dt,
                                     DomainPartition * const domain ) override;

  virtual real64 SolverStep( real64 const & time_n,
                             real64 const & dt,
                             int const cycleNumber,
                             DomainPartition * const domain ) override;

  FlowSolverBase * GetFlowSolver() const { return m_flowSolver; }

  WellSolverBase * GetWellSolver() const { return m_wellSolver; }
  
  struct viewKeyStruct : SolverBase::viewKeyStruct
  {

    // solver that assembles the reservoir equations
    constexpr static auto flowSolverNameString = "flowSolverName";

    // solver that assembles the well
    constexpr static auto wellSolverNameString = "wellSolverName";

  } reservoirWellsSolverViewKeys;


protected:

  virtual void PostProcessInput() override;

private:

  // solver that assembles the reservoir equations
  string m_flowSolverName;

  // solver that assembles the well equations and compute perforation rates
  string m_wellSolverName;

  // pointer to the flow sub-solver
  FlowSolverBase * m_flowSolver;

  // pointer to the well sub-solver
  WellSolverBase * m_wellSolver;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_COUPLEDSOLVERS_RESERVOIRSOLVER_HPP_ */
