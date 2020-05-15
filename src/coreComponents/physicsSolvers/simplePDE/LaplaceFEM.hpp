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

#ifndef GEOSX_PHYSICSSOLVERS_LAGRANGIAN_FEM_HPP_
#define GEOSX_PHYSICSSOLVERS_LAGRANGIAN_FEM_HPP_

#include "physicsSolvers/SolverBase.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
//START_SPHINX_INCLUDE_00
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
//END_SPHINX_INCLUDE_00

//START_SPHINX_INCLUDE_02
class LaplaceFEM : public SolverBase
{
public:

  LaplaceFEM() = delete;

  LaplaceFEM( const std::string & name,
              Group * const parent );

  virtual ~LaplaceFEM() override;

  static string CatalogName() { return "LaplaceFEM"; }

  virtual void RegisterDataOnMesh( Group * const MeshBodies ) override final;

//END_SPHINX_INCLUDE_02
/**
 * @defgroup Solver Interface Functions
 *
 * These functions provide the primary interface that is required for derived classes
 */
/**@{*/

  //START_SPHINX_INCLUDE_03
  virtual real64 SolverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition * domain ) override;

  virtual real64 ExplicitStep( real64 const & time_n,
                               real64 const & dt,
                               integer const cycleNumber,
                               DomainPartition * const domain ) override;

  void sparsityGeneration( DomainPartition const & domain );

  virtual void
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition * const domain,
                     DofManager & dofManager,
                     ParallelMatrix & matrix,
                     ParallelVector & rhs,
                     ParallelVector & solution ) override;

  virtual void
  SetupDofs( DomainPartition const * const domain,
             DofManager & dofManager ) const override;

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
  ResetStateToBeginningOfStep( DomainPartition * const GEOSX_UNUSED_PARAM( domain ) ) override
  {}

  virtual void
  ImplicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition * const domain ) override;

  //END_SPHINX_INCLUDE_03
  /**@}*/

  void ApplyDirichletBC_implicit( real64 const time,
                                  DofManager const & dofManager,
                                  DomainPartition & domain,
                                  ParallelMatrix & matrix,
                                  ParallelVector & rhs );

  //START_SPHINX_INCLUDE_01
  enum class timeIntegrationOption
  {
    SteadyState,
    ImplicitTransient,
    ExplicitTransient
  };
  //END_SPHINX_INCLUDE_01

  //START_SPHINX_INCLUDE_04
  struct viewKeyStruct : public SolverBase::viewKeyStruct
  {
    dataRepository::ViewKey timeIntegrationOption = { "timeIntegrationOption" };
    dataRepository::ViewKey fieldVarName = { "fieldName" };

  } laplaceFEMViewKeys;
  //END_SPHINX_INCLUDE_04

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
  timeIntegrationOption m_timeIntegrationOption;

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_LAGRANGIAN_FEM_HPP_ */
