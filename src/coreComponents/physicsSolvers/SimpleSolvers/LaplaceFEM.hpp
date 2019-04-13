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

/*
 * NewtonianMechanics.hpp
 *
 *  Created on: Dec 4, 2014
 *      Author: rrsettgast
 */

#ifndef SOLID_MECHANICS_LAGRANGIAN_FEM_HPP_
#define SOLID_MECHANICS_LAGRANGIAN_FEM_HPP_

#include "physicsSolvers/SolverBase.hpp"
#include "systemSolverInterface/LinearSolverWrapper.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"

#include "DofManager.hpp"
#include "DofManager.hpp"

// Just to print the matrix
#include "TrilinosInterface.hpp"

struct stabledt
{
  double m_maxdt;
};

namespace ML_Epetra
{ class MultiLevelPreconditioner; }

namespace geosx
{
namespace dataRepository
{
class ManagedGroup;
}
class FieldSpecificationBase;
class FiniteElementBase;
class DomainPartition;

using ParallelMatrix = typename TrilinosInterface::ParallelMatrix;
using ParallelVector = typename TrilinosInterface::ParallelVector;

class LaplaceFEM : public SolverBase
{
public:
  LaplaceFEM( const std::string& name,
              ManagedGroup * const parent );


  virtual ~LaplaceFEM() override;

  static string CatalogName() { return "LaplaceFEM"; }

  virtual void RegisterDataOnMesh( ManagedGroup * const MeshBodies ) override final;

  virtual void InitializePreSubGroups(ManagedGroup * const rootGroup) override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual real64 SolverStep( real64 const& time_n,
                         real64 const& dt,
                         integer const cycleNumber,
                         DomainPartition * domain ) override;

  virtual real64 ExplicitStep( real64 const & time_n,
                                 real64 const & dt,
                                 integer const cycleNumber,
                                 DomainPartition * const domain ) override;

  virtual void ImplicitStepSetup( real64 const& time_n,
                              real64 const& dt,
                              DomainPartition * const domain,
                              systemSolverInterface::EpetraBlockSystem * const blockSystem ) override;


  virtual void AssembleSystem( DomainPartition * const domain,
                               systemSolverInterface::EpetraBlockSystem * const blockSystem,
                               real64 const time,
                               real64 const dt ) override;

  virtual void ApplyBoundaryConditions( DomainPartition * const domain,
                                        systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                        real64 const time,
                                        real64 const dt ) override;

//  virtual real64
//  CalculateResidualNorm( systemSolverInterface::EpetraBlockSystem const * const blockSystem ) override;

  virtual void SolveSystem( systemSolverInterface::EpetraBlockSystem * const blockSystem,
                            SystemSolverParameters const * const params ) override;

  virtual void
  ApplySystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                       real64 const scalingFactor,
                       DomainPartition * const domain ) override;

  virtual void ResetStateToBeginningOfStep( DomainPartition * const domain ) override {}

  virtual  void ImplicitStepComplete( real64 const & time,
                                      real64 const & dt,
                                      DomainPartition * const domain ) override;
  /**@}*/


  void TimeStepQuasiStatic( real64 const& time_n,
                            real64 const& dt,
                            integer const cycleNumber,
                            DomainPartition& domain );

//  real64 TimeStepImplicit( real64 const & time_n,
//                           real64 const & dt,
//                           integer const cycleNumber,
//                           DomainPartition * const domain );

  void SetupSystem( DomainPartition * const domain,
                    systemSolverInterface::EpetraBlockSystem * const blockSystem );

  void SetupMLPreconditioner( DomainPartition const & domain,
                              ML_Epetra::MultiLevelPreconditioner* MLPrec );

  void ApplyDirichletBC_implicit( real64 const time,
                                  DomainPartition & domain,
                                  systemSolverInterface::EpetraBlockSystem & blockSystem );

  void ApplyDirichletBC_implicit( real64 const time,
                                  DomainPartition & domain,
                                  ParallelMatrix & matrix,
                                  ParallelVector & rhs );

  /////////////////////////////////////////////////////////////////////////////////////////
  void ApplyBoundaryConditionToSystem( FieldSpecificationManager const & fsManager,
                                       string const & functionName,
                                       set<localIndex> const & targetSet,
                                       real64 const time,
                                       dataRepository::ManagedGroup * dataGroup,
                                       string const & fieldName,
                                       string const & dofMapName,
                                       integer const & dofDim,
                                       integer const & component,
                                       real64 const & scale,
                                       ParallelMatrix & matrix,
                                       ParallelVector & rhs );

  template< typename LAMBDA >
  void ApplyBoundaryConditionToSystem( FieldSpecificationManager const & fsManager,
                                       string const & functionName,
                                       set<localIndex> const & targetSet,
                                       real64 const time,
                                       dataRepository::ManagedGroup * dataGroup,
                                       arrayView1d<globalIndex const> const & dofMap,
                                       integer const & dofDim,
                                       integer const & component,
                                       real64 const & scale,
                                       ParallelMatrix & matrix,
                                       ParallelVector & rhs,
                                       LAMBDA && lambda )
  {
    NewFunctionManager * functionManager = NewFunctionManager::Instance();

    globalIndex_array dof( targetSet.size() );
    real64_array rhsContribution( targetSet.size() );

    if( functionName.empty() )
    {

      integer counter=0;
      for( auto a : targetSet )
      {
        dof( counter ) = dofDim*dofMap[a]+component;
        SpecifyFieldValue( dof( counter ),
                           matrix,
                           rhsContribution( counter ),
                           scale,
                           lambda( a ) );

        ++counter;
      }
      ReplaceGlobalValues( rhs, counter, dof.data(), rhsContribution.data() );
    }
    else
    {
      FunctionBase const * const function  = functionManager->GetGroup<FunctionBase>( functionName );

      GEOS_ERROR_IF( function == nullptr, "Function '" << functionName << "' not found" );

      if( function->isFunctionOfTime()==2 )
      {
        real64 value = scale * function->Evaluate( &time );
        integer counter=0;
        for( auto a : targetSet )
        {
          dof( counter ) = dofDim*integer_conversion<int>( dofMap[a] )+component;
          SpecifyFieldValue( dof( counter ),
                             matrix,
                             rhsContribution( counter ),
                             value,
                             lambda( a ) );
          ++counter;
        }
        ReplaceGlobalValues( rhs, counter, dof.data(), rhsContribution.data() );
      }
      else
      {
        real64_array result;
        result.resize( integer_conversion<localIndex>( targetSet.size()));
        function->Evaluate( dataGroup, time, targetSet, result );
        integer counter=0;
        for( auto a : targetSet )
        {
          dof( counter ) = dofDim*integer_conversion<int>( dofMap[a] )+component;
          SpecifyFieldValue( dof( counter ),
                             matrix,
                             rhsContribution( counter ),
                             scale*result[counter],
                             lambda( a ) );
          ++counter;
        }
        ReplaceGlobalValues( rhs, counter, dof.data(), rhsContribution.data() );
      }
    }
  }

  static inline void SpecifyFieldValue( globalIndex const dof,
                                        ParallelMatrix & matrix,
                                        real64 & rhs,
                                        real64 const & bcValue,
                                        real64 const fieldValue )
  {
    if( true )//node_is_ghost[*nd] < 0 )
    {
      //real64 LARGE = blockSystem->ClearSystemRow( blockID, static_cast< int >( dof ), 1.0 );
      //rhs = -LARGE*( bcValue - fieldValue );
      if( matrix.ParallelMatrixGetLocalRowID( dof ) >= 0 )
      {
        matrix.clearRow( dof, 1.0 );
        rhs = bcValue;
      }
      else
      {
        rhs = 0.0;
      }
    }
  }

  static inline void ReplaceGlobalValues( ParallelVector & rhs,
                                          localIndex const num,
                                          globalIndex const * const dof,
                                          real64 const * const values )
  {
    rhs.set( dof, values, num );
  }

  void solve( ParallelMatrix const * const matrix,
              ParallelVector * const rhs,
              ParallelVector * const solution,
              SystemSolverParameters const * const params );
  /////////////////////////////////////////////////////////////////////////////////////////

  enum class timeIntegrationOption
  {
    SteadyState,
    ImplicitTransient,
    ExplicitTransient
  };

  struct viewKeyStruct : public SolverBase::viewKeyStruct
  {
    constexpr static auto blockLocalDofNumberString = "blockLocalDofNumber_Laplace";

    dataRepository::ViewKey timeIntegrationOption = { "timeIntegrationOption" };
    dataRepository::ViewKey fieldVarName = { "fieldName" };

    dataRepository::ViewKey blockLocalDofNumber = { blockLocalDofNumberString };

  } laplaceFEMViewKeys;

  inline ParallelVector const * getSolution() const {
    return & m_solution;
  }

protected:
  virtual void PostProcessInput() override final;

private:
  string m_fieldName;
  stabledt m_stabledt;
  timeIntegrationOption m_timeIntegrationOption;
  LaplaceFEM();
  DofManager dofManager;
  bool firstTime = true;

  int mpiRank = CommunicationTools::MPI_Rank( MPI_COMM_GEOSX );

  ParallelMatrix m_matrix;
  ParallelVector m_rhs;
  ParallelVector m_solution;
};

} /* namespace geosx */

#endif /* SOLID_MECHANICS_LAGRANGIAN_FEM_HPP_ */
