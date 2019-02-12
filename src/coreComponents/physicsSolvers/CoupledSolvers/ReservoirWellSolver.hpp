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
 * @file ReservoirWellSolver.hpp
 *
 */

#ifndef RESERVOIRWELLSOLVER_HPP_
#define RESERVOIRWELLSOLVER_HPP_

#include "../SolverBase.hpp"

namespace geosx
{

class ReservoirWellSolver : public SolverBase
{
public:
  ReservoirWellSolver( const std::string& name,
                       ManagedGroup * const parent );
  ~ReservoirWellSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return "ReservoirWell"; }

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
    constexpr static auto flowSolverNameString = "flowSolverName";
    constexpr static auto wellSolverNameString = "wellSolverName";
  } reservoirWellSolverViewKeys;


protected:
  virtual void PostProcessInput() override final;

  virtual void InitializePostInitialConditions_PreSubGroups(dataRepository::ManagedGroup * const problemManager) override final;


private:

  string m_flowSolverName;
  string m_wellSolverName;

};

} /* namespace geosx */

#endif /* RESERVOIRWELLSOLVER_HPP_ */
