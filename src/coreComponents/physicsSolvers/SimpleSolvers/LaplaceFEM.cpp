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
 * @file LaplaceFEM.cpp
 *
 */

#include "LaplaceFEM.hpp"

#include <vector>
#include <math.h>

#include "common/TimingMacros.hpp"

#include "dataRepository/ManagedGroup.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteElement/ElementLibrary/FiniteElement.h"
#include "finiteElement/Kinematics.h"
#include "managers/NumericalMethodsManager.hpp"
#include "codingUtilities/Utilities.hpp"

#include "managers/DomainPartition.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{}
}

using namespace dataRepository;
using namespace constitutive;

LaplaceFEM::LaplaceFEM( const std::string& name,
                        ManagedGroup * const parent ):
  SolverBase( name, parent ),
  m_fieldName("primaryField")
{
  RegisterViewWrapper<string>(laplaceFEMViewKeys.timeIntegrationOption.Key())->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("option for default time integration method");

  RegisterViewWrapper<string>(laplaceFEMViewKeys.fieldVarName.Key(), &m_fieldName, false)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("name of field variable");
}

LaplaceFEM::~LaplaceFEM()
{
  // TODO Auto-generated destructor stub
}

void LaplaceFEM::RegisterDataOnMesh( ManagedGroup * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    NodeManager * const nodes = mesh.second->group_cast<MeshBody*>()->getMeshLevel(0)->getNodeManager();

    nodes->RegisterViewWrapper<real64_array >( m_fieldName )->
      setApplyDefaultValue(0.0)->
      setPlotLevel(PlotLevel::LEVEL_0)->
      setDescription("Primary field variable");
  }
}

void LaplaceFEM::PostProcessInput()
{
  SolverBase::PostProcessInput();

  string tiOption = this->getReference<string>(laplaceFEMViewKeys.timeIntegrationOption);

  if( tiOption == "SteadyState" )
  {
    this->m_timeIntegrationOption = timeIntegrationOption::SteadyState;
  }
  else if( tiOption == "ImplicitTransient" )
  {
    this->m_timeIntegrationOption = timeIntegrationOption::ImplicitTransient;
  }
  else if ( tiOption == "ExplicitTransient" )
  {
    this->m_timeIntegrationOption = timeIntegrationOption::ExplicitTransient;
  }
  else
  {
    GEOS_ERROR("invalid time integration option");
  }

  // Set basic parameters for solver
  m_linearSolverParameters.verbosity = 0;
  m_linearSolverParameters.solverType = "gmres";
  m_linearSolverParameters.krylov.tolerance = 1e-8;
  m_linearSolverParameters.krylov.maxIterations = 250;
  m_linearSolverParameters.krylov.maxRestart = 250;
  m_linearSolverParameters.preconditionerType = "amg";
  m_linearSolverParameters.amg.smootherType = "gaussSeidel";
  m_linearSolverParameters.amg.coarseType = "direct";
}

real64 LaplaceFEM::SolverStep( real64 const& time_n,
                               real64 const& dt,
                               const int cycleNumber,
                               DomainPartition * domain )
{
  real64 dtReturn = dt;
  if( m_timeIntegrationOption == timeIntegrationOption::ExplicitTransient )
  {
    dtReturn = ExplicitStep( time_n, dt, cycleNumber, domain );
  }
  else if( m_timeIntegrationOption == timeIntegrationOption::ImplicitTransient ||
           m_timeIntegrationOption == timeIntegrationOption::SteadyState )
  {
    dtReturn = this->LinearImplicitStep( time_n, dt, cycleNumber, domain, m_dofManager, m_matrix, m_rhs, m_solution );
  }
  return dtReturn;
}

real64 LaplaceFEM::ExplicitStep( real64 const& time_n,
                                 real64 const& dt,
                                 const int cycleNumber,
                                 DomainPartition * const domain )
{
  return dt;
}

void LaplaceFEM::ImplicitStepSetup( real64 const & time_n,
                                    real64 const & dt,
                                    DomainPartition * const domain,
                                    DofManager & dofManager,
                                    ParallelMatrix & matrix,
                                    ParallelVector & rhs,
                                    ParallelVector & solution )
{
  // Computation of the sparsity pattern
  SetupSystem( domain, dofManager, matrix, rhs, solution );
}

void LaplaceFEM::ImplicitStepComplete( real64 const & time_n,
                                       real64 const & dt,
                                       DomainPartition * const domain)
{
}

void LaplaceFEM::SetupSystem( DomainPartition * const domain,
                              DofManager & dofManager,
                              ParallelMatrix & matrix,
                              ParallelVector & rhs,
                              ParallelVector & solution )
{
  m_dofManager.setMesh( domain, 0, 0 );
  m_dofManager.addField( m_fieldName,
                         DofManager::Location::Node,
                         DofManager::Connectivity::Elem );

  m_dofManager.setSparsityPattern( matrix, m_fieldName, m_fieldName );
  m_dofManager.setVector( rhs, m_fieldName, m_fieldName );
  m_dofManager.setVector( solution, m_fieldName, m_fieldName );
}

void LaplaceFEM::AssembleSystem( real64 const time_n,
                                 real64 const dt,
                                 DomainPartition * const domain,
                                 DofManager const & dofManager,
                                 ParallelMatrix & matrix,
                                 ParallelVector & rhs )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  //ManagedGroup * const nodeManager = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();
  NumericalMethodsManager const *
  numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);
  FiniteElementDiscretizationManager const *
  feDiscretizationManager = numericalMethodManager->
    GetGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);

  //globalIndex_array const & indexArray = nodeManager->getReference<globalIndex_array>( dofManager.getKey( m_fieldName ) );

  // Initialize all entries to zero
  matrix.zero();
  rhs.zero();

  matrix.open();
  rhs.open();

  // begin region loop
  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegion * const elementRegion = elemManager->GetRegion(er);

    FiniteElementDiscretization const *
    feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

    elementRegion->forElementSubRegionsIndex<CellElementSubRegion>([&]( localIndex const esr,
                                                                        CellElementSubRegion const * const elementSubRegion )
    {
      array3d<R1Tensor> const &
      dNdX = elementSubRegion->getReference< array3d< R1Tensor > >(keys::dNdX);

      arrayView2d<real64> const &
      detJ = elementSubRegion->getReference< array2d<real64> >(keys::detJ);

      arrayView2d<localIndex> const & elemsToNodes = elementSubRegion->nodeList();
      const int numNodesPerElement = integer_conversion<int>(elemsToNodes.size(1));

      globalIndex_array element_index( numNodesPerElement );
      real64_array element_rhs( numNodesPerElement );
      real64_array2d element_matrix( numNodesPerElement, numNodesPerElement );

      integer_array const & elemGhostRank = elementSubRegion->m_ghostRank;
      const int n_q_points = feDiscretization->m_finiteElement->n_quadrature_points();

      // begin element loop, skipping ghost elements
      for( localIndex k=0 ; k<elementSubRegion->size() ; ++k )
      {
        if(elemGhostRank[k] < 0)
        {
          dofManager.getIndices( element_index, DofManager::Connectivity::Elem, er, esr, k, m_fieldName );

          element_rhs = 0.0;
          element_matrix = 0.0;
          for( localIndex q=0 ; q<n_q_points ; ++q)
          {
            for( localIndex a=0 ; a<numNodesPerElement ; ++a)
            {
              real64 diffusion = 1.0;
              for( localIndex b=0 ; b<numNodesPerElement ; ++b)
              {
                element_matrix(a,b) += detJ[k][q] *
                                       diffusion *
                                     + Dot( dNdX[k][q][a], dNdX[k][q][b] );
              }

            }
          }
          matrix.add( element_index, element_index, element_matrix );
          rhs.add( element_index, element_rhs );
        }
      }
    });
  }
  matrix.close();
  rhs.close();

  if( verboseLevel() == 2 )
  {
    GEOS_LOG_RANK_0( "After LaplaceFEM::AssembleSystem" );
    GEOS_LOG_RANK_0("\nJacobian:\n" << matrix);
    GEOS_LOG_RANK_0("\nResidual:\n" << rhs);
  }

  if( verboseLevel() >= 3 )
  {
    SystemSolverParameters * const solverParams = getSystemSolverParameters();
    integer newtonIter = solverParams->numNewtonIterations();

    string filename_mat = "matrix_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    matrix.write( filename_mat, true );

    string filename_rhs = "rhs_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    rhs.write( filename_rhs, true );

    GEOS_LOG_RANK_0( "After LaplaceFEM::AssembleSystem" );
    GEOS_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOS_LOG_RANK_0( "Residual: written to " << filename_rhs );
  }
}

void LaplaceFEM::ApplySystemSolution( DofManager const & dofManager,
                                      ParallelVector const & solution,
                                      real64 const scalingFactor,
                                      DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();

  dofManager.copyVectorToField( solution, m_fieldName, nodeManager );

  // Syncronize ghost nodes
  std::map<string, string_array> fieldNames;
  fieldNames["node"].push_back( m_fieldName );

  CommunicationTools::
  SynchronizeFields( fieldNames, mesh,
                     domain->getReference<array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  if( dofManager.needDoubleSync() )
  {
    CommunicationTools::
    SynchronizeFields( fieldNames, mesh,
                       domain->getReference<array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );
  }
}

void LaplaceFEM::ApplyBoundaryConditions( real64 const time_n,
                                          real64 const dt,
                                          DomainPartition * const domain,
                                          DofManager const & dofManager,
                                          ParallelMatrix & matrix,
                                          ParallelVector & rhs )
{
  ApplyDirichletBC_implicit( time_n + dt, dofManager, *domain, m_matrix, m_rhs );

  if( verboseLevel() == 2 )
  {
    GEOS_LOG_RANK_0( "After LaplaceFEM::ApplyBoundaryConditions" );
    GEOS_LOG_RANK_0("\nJacobian:\n" << matrix);
    GEOS_LOG_RANK_0("\nResidual:\n" << rhs);
  }

  if( verboseLevel() >= 3 )
  {
    SystemSolverParameters * const solverParams = getSystemSolverParameters();
    integer newtonIter = solverParams->numNewtonIterations();

    string filename_mat = "matrix_bc_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    matrix.write( filename_mat, true );

    string filename_rhs = "rhs_bc_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    rhs.write( filename_rhs, true );

    GEOS_LOG_RANK_0( "After LaplaceFEM::ApplyBoundaryConditions" );
    GEOS_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOS_LOG_RANK_0( "Residual: written to " << filename_rhs );
  }
}

void LaplaceFEM::SolveSystem( DofManager const & dofManager,
                              ParallelMatrix & matrix,
                              ParallelVector & rhs,
                              ParallelVector & solution )
{
  rhs.scale( -1.0 ); // TODO decide if we want this here
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );

  if( verboseLevel() == 2 )
  {
    GEOS_LOG_RANK_0("After LaplaceFEM::SolveSystem");
    GEOS_LOG_RANK_0("\nSolution\n" << solution);
  }
}

void LaplaceFEM::ApplyDirichletBC_implicit( real64 const time,
                                            DofManager const & dofManager,
                                            DomainPartition & domain,
                                            ParallelMatrix & matrix,
                                            ParallelVector & rhs )
{
  FieldSpecificationManager const * const fsManager = FieldSpecificationManager::get();

  fsManager->Apply( time,
                    &domain,
                    "nodeManager",
                    m_fieldName,
                    [&]( FieldSpecificationBase const * const bc,
                    string const &,
                    set<localIndex> const & targetSet,
                    ManagedGroup * const targetGroup,
                    string const fieldName )->void
  {
    bc->ApplyBoundaryConditionToSystem<FieldSpecificationEqual, LAInterface>( targetSet,
                                                                              false,
                                                                              time,
                                                                              targetGroup,
                                                                              m_fieldName,
                                                                              dofManager.getKey( m_fieldName ),
                                                                              1,
                                                                              matrix,
                                                                              rhs );
  });
}

REGISTER_CATALOG_ENTRY( SolverBase, LaplaceFEM, std::string const &, ManagedGroup * const )
} /* namespace ANST */
