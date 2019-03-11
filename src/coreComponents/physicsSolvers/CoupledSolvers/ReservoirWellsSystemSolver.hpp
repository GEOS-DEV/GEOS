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
 * @file ReservoirWellsSystemSolver.hpp
 *
 */

#ifndef RESERVOIRWELLSYSTEMSOLVER_HPP_
#define RESERVOIRWELLSYSTEMSOLVER_HPP_

#include "../SolverBase.hpp"

namespace geosx
{

class ReservoirWellsSystemSolver : public SolverBase
{
public:
  ReservoirWellsSystemSolver( const std::string& name,
                              ManagedGroup * const parent );
  ~ReservoirWellsSystemSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return "ReservoirWellsSystem"; }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/
  
  virtual void ImplicitStepSetup( real64 const& time_n,
                                  real64 const& dt,
                                  DomainPartition * const domain,
                                  systemSolverInterface::EpetraBlockSystem * const blockSystem) override final;

  virtual void AssembleSystem( DomainPartition * const domain,
                               systemSolverInterface::EpetraBlockSystem * const blockSystem,
                               real64 const time,
                               real64 const dt ) override;

  virtual void ApplyBoundaryConditions( DomainPartition * const domain,
                                        systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                        real64 const time,
                                        real64 const dt ) override;

  virtual real64
  CalculateResidualNorm( systemSolverInterface::EpetraBlockSystem const *const blockSystem,
                         DomainPartition * const domain ) override;

  virtual void SolveSystem( systemSolverInterface::EpetraBlockSystem * const blockSystem,
                            SystemSolverParameters const * const params ) override;

  virtual bool
  CheckSystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                       real64 const scalingFactor,
                       DomainPartition * const domain ) override;

  virtual void
  ApplySystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                       real64 const scalingFactor,
                       DomainPartition * const domain ) override;

  virtual void ResetStateToBeginningOfStep( DomainPartition * const domain ) override;

  
  virtual void ImplicitStepComplete( real64 const& time_n,
                                     real64 const& dt,
                                     DomainPartition * const domain) override;

  virtual real64 SolverStep( real64 const & time_n,
                             real64 const & dt,
                             int const cycleNumber,
                             DomainPartition * domain ) override;
  
  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    // solver that assembles the reservoir equations
    constexpr static auto flowSolverNameString = "flowSolverName";
    // solver that assembles the well
    constexpr static auto wellSolverNameString = "wellSolverName";
  } reservoirWellsSystemSolverViewKeys;


protected:
  virtual void InitializePostInitialConditions_PreSubGroups(dataRepository::ManagedGroup * const problemManager) override final;


private:
  
  /**
   * @brief Set up the linear system (DOF indices and sparsity patterns)
   * @param domain the domain containing the mesh and fields
   * @param blockSystem the linear system object
   */
  void SetupSystem ( DomainPartition * const domain,
                     systemSolverInterface::EpetraBlockSystem * const blockSystem );

  // solver that assembles the reservoir equations
  string m_flowSolverName;
  // solver that assembles the well
  string m_wellSolverName;

};

} /* namespace geosx */

#endif /* RESERVOIRWELLSSYSTEMSOLVER_HPP_ */
