/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_FEM_HPP_
#define GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_FEM_HPP_

#include "common/EnumStrings.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

namespace geosx
{

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
                             DomainPartition & domain ) override;

  virtual void
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  SetupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  SetupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               array1d< real64 > & localRhs,
               array1d< real64 > & localSolution,
               bool const setSparsity = false ) override;

  virtual void
  AssembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void
  ApplyBoundaryConditions( real64 const time,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  virtual void
  SolveSystem( DofManager const & dofManager,
               ParallelMatrix & matrix,
               ParallelVector & rhs,
               ParallelVector & solution ) override;

  virtual void
  ApplySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void
    ResetStateToBeginningOfStep( DomainPartition & GEOSX_UNUSED_PARAM( domain ) ) override;

  virtual void
  ImplicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  //END_SPHINX_INCLUDE_03
  /**@}*/

  void ApplyDirichletBC_implicit( real64 const time,
                                  DofManager const & dofManager,
                                  DomainPartition & domain,
                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                  arrayView1d< real64 > const & localRhs );

  //START_SPHINX_INCLUDE_01
  enum class TimeIntegrationOption
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

private:

  string m_fieldName;
  TimeIntegrationOption m_timeIntegrationOption;

};

//START_SPHINX_INCLUDE_05
ENUM_STRINGS( LaplaceFEM::TimeIntegrationOption, "SteadyState", "ImplicitTransient", "ExplicitTransient" )
//END_SPHINX_INCLUDE_05

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_FEM_HPP_ */
