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

#ifndef SOLID_MECHANICS_LAGRANGIAN_FEM_HPP_
#define SOLID_MECHANICS_LAGRANGIAN_FEM_HPP_

#include "physicsSolvers/SolverBase.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"

#include "linearAlgebraInterface/src/DofManager.hpp"
#include "linearAlgebraInterface/src/InterfaceTypes.hpp"

struct stabledt
{
  double m_maxdt;
};

namespace geosx
{
namespace dataRepository
{
class ManagedGroup;
}
class FieldSpecificationBase;
class FiniteElementBase;
class DomainPartition;


class LaplaceFEM : public SolverBase
{
public:

  LaplaceFEM( const std::string& name,
              ManagedGroup * const parent );

  virtual ~LaplaceFEM() override;

  static string CatalogName() { return "LaplaceFEM"; }

  virtual void RegisterDataOnMesh( ManagedGroup * const MeshBodies ) override final;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual real64 SolverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition * domain ) override;

  virtual real64 ExplicitStep( real64 const & time_n,
                               real64 const & dt,
                               integer const cycleNumber,
                               DomainPartition * const domain ) override;

  virtual void
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition * const domain,
                     DofManager & dofManager,
                     ParallelMatrix & matrix,
                     ParallelVector & rhs,
                     ParallelVector & solution ) override;


  virtual void
  AssembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition * const domain,
                  DofManager const & dofManager,
                  ParallelMatrix & matrix,
                  ParallelVector & rhs ) override;

  virtual void
  ApplyBoundaryConditions( real64 const time,
                           real64 const dt,
                           DomainPartition * const domain,
                           DofManager const & dofManager,
                           ParallelMatrix & matrix,
                           ParallelVector & rhs ) override;

  virtual void
  SolveSystem( DofManager const & dofManager,
               ParallelMatrix & matrix,
               ParallelVector & rhs,
               ParallelVector & solution ) override;

  virtual void
  ApplySystemSolution( DofManager const & dofManager,
                       ParallelVector const & solution,
                       real64 const scalingFactor,
                       DomainPartition * const domain ) override;

  virtual void
  ResetStateToBeginningOfStep( DomainPartition * const domain ) override
  {}

  virtual void
  ImplicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition * const domain ) override;
  /**@}*/

  void SetupSystem( DomainPartition * const domain,
                    DofManager & dofManager,
                    ParallelMatrix & matrix,
                    ParallelVector & rhs,
                    ParallelVector & solution );

  void ApplyDirichletBC_implicit( real64 const time,
                                  DofManager const & dofManager,
                                  DomainPartition & domain,
                                  ParallelMatrix & matrix,
                                  ParallelVector & rhs );

  enum class timeIntegrationOption
  {
    SteadyState,
    ImplicitTransient,
    ExplicitTransient
  };

  struct viewKeyStruct : public SolverBase::viewKeyStruct
  {
    dataRepository::ViewKey timeIntegrationOption = { "timeIntegrationOption" };
    dataRepository::ViewKey fieldVarName = { "fieldName" };

  } laplaceFEMViewKeys;

  inline ParallelVector const * getSolution() const {
    return & m_solution;
  }

  inline globalIndex getSize() const {
    return m_matrix.globalRows();
  }

protected:
  virtual void PostProcessInput() override final;

private:

  string m_fieldName;
  stabledt m_stabledt;
  timeIntegrationOption m_timeIntegrationOption;
  LaplaceFEM();

};

} /* namespace geosx */

#endif /* SOLID_MECHANICS_LAGRANGIAN_FEM_HPP_ */
