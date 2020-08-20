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

#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "physicsSolvers/SolverBase.hpp"

struct stabledt
{
  double m_maxdt;
};

namespace geosx
{
namespace dataRepository
{
class Group;
}
class FieldSpecificationBase;
class FiniteElementBase;
class DomainPartition;

class ReactionDiffusionFEM : public SolverBase
{
public:
  ReactionDiffusionFEM( const std::string & name, Group * const parent );

  virtual ~ReactionDiffusionFEM() override;

  static string CatalogName()
  {
    return "ReactionDiffusionFEM";
  }

  virtual void RegisterDataOnMesh( Group * const MeshBodies ) override final;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived
   * classes
   */
  /**@{*/

  virtual real64 SolverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override;

  virtual real64 ExplicitStep( real64 const & time_n,
                               real64 const & dt,
                               integer const cycleNumber,
                               DomainPartition & domain ) override;

  virtual void ImplicitStepSetup( real64 const & time_n,
                                  real64 const & dt,
                                  DomainPartition & domain ) override;

  virtual void SetupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override;

  virtual void AssembleSystem( real64 const time, real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

  virtual void ApplyBoundaryConditions( real64 const time, real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs ) override;

  virtual void SolveSystem( DofManager const & dofManager, ParallelMatrix & matrix,
                            ParallelVector & rhs,
                            ParallelVector & solution ) override;

  virtual void ApplySystemSolution( DofManager const & dofManager,
                                    arrayView1d< real64 const > const & localSolution,
                                    real64 const scalingFactor,
                                    DomainPartition & domain ) override;

  virtual void ResetStateToBeginningOfStep(
    DomainPartition & GEOSX_UNUSED_PARAM( domain ) ) override
  {}

  virtual void ImplicitStepComplete( real64 const & time, real64 const & dt,
                                     DomainPartition & domain ) override;
  /**@}*/

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
    // static constexpr auto coeffFieldName = "coeffFieldName";
    // static constexpr auto coeffName = "coeffField";

    dataRepository::ViewKey timeIntegrationOption =
    { "timeIntegrationOption" };
    dataRepository::ViewKey fieldVarName =
    { "fieldName" };

  } reactionDiffusionFEMViewKeys;

  inline ParallelVector const * getSolution() const
  {
    return &m_solution;
  }

  inline globalIndex getSize() const
  {
    return m_matrix.numGlobalRows();
  }

protected:
  virtual void PostProcessInput() override final;

private:
  string m_fieldName;
  stabledt m_stabledt;
  timeIntegrationOption m_timeIntegrationOption;

//  string m_coeffFieldName;

  ReactionDiffusionFEM();
};

} /* namespace geosx */

#endif /* SOLID_MECHANICS_LAGRANGIAN_FEM_HPP_ */
