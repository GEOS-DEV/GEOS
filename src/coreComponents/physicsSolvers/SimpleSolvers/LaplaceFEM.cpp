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

#include "RAJA/RAJA.hpp"
#include "RAJA/util/defines.hpp"

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
using namespace systemSolverInterface;

LaplaceFEM::LaplaceFEM( const std::string& name,
                        ManagedGroup * const parent ):
  SolverBase( name, parent )
{
//  this->RegisterGroup<SystemSolverParameters>( groupKeys.systemSolverParameters.Key() );
  // To generate the schema, multiple solvers of that use this command are constructed
  // Doing this can cause an error in the block setup, so move it to InitializePreSubGroups
  // getLinearSystemRepository()->SetBlockID( BlockIDs::dummyScalarBlock, this->getName() );

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

    nodes->RegisterViewWrapper<array1d<globalIndex> >( viewKeyStruct::blockLocalDofNumberString )->
      setApplyDefaultValue(-1)->
      setPlotLevel(PlotLevel::LEVEL_1)->
      setDescription("Global DOF numbers for the primary field variable");
  }
}

void LaplaceFEM::PostProcessInput()
{
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
}

void LaplaceFEM::InitializePreSubGroups( ManagedGroup * const problemManager )
{
  SolverBase::InitializePreSubGroups(problemManager);

  // set the blockID for the block system interface
  getLinearSystemRepository()->SetBlockID( BlockIDs::dummyScalarBlock, this->getName() );
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
    dtReturn = this->LinearImplicitStep( time_n, dt, cycleNumber, domain, getLinearSystemRepository() );
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

void LaplaceFEM::ImplicitStepSetup( real64 const& time_n,
                                    real64 const& dt,
                                    DomainPartition * const domain,
                                    systemSolverInterface::EpetraBlockSystem * const blockSystem )
{
  // Computation of the sparsity pattern
  SetupSystem( domain, blockSystem );
}

void LaplaceFEM::ImplicitStepComplete( real64 const & time_n,
                                       real64 const & dt,
                                       DomainPartition * const domain)
{
}

void LaplaceFEM::SetupSystem( DomainPartition * const domain,
                              EpetraBlockSystem * const blockSystem )
{
  // Set basic parameters for solver
  m_parameters.verbosity = 0;
  m_parameters.solverType = "gmres";
  m_parameters.krylov.tolerance = 1e-8;
  m_parameters.krylov.maxIterations = 250;
  m_parameters.krylov.maxRestart = 250;
  m_parameters.preconditionerType = "amg";
  m_parameters.amg.smootherType = "gaussSeidel";
  m_parameters.amg.coarseType = "direct";

  dofManager.setMesh( domain, 0, 0 );
  dofManager.addField( m_fieldName, DofManager::Location::Node, DofManager::Connectivity::Elem );

  ParallelMatrix & sparsity = m_matrix;
  dofManager.setSparsityPattern( sparsity, m_fieldName, m_fieldName );
  dofManager.setVector( m_rhs, m_fieldName, m_fieldName );
  dofManager.setVector( m_solution, m_fieldName, m_fieldName );
}

void LaplaceFEM::AssembleSystem ( DomainPartition * const  domain,
                                  EpetraBlockSystem * const blockSystem,
                                  real64 const time_n,
                                  real64 const dt )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ManagedGroup * const nodeManager = mesh->getNodeManager();
  ConstitutiveManager  * const constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  NumericalMethodsManager const *
  numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);
  FiniteElementDiscretizationManager const *
  feDiscretizationManager = numericalMethodManager->
    GetGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);

  globalIndex_array const & indexArray = nodeManager->getReference<globalIndex_array>( dofManager.getKey( m_fieldName ) );

  // Initialize all entries to zero
  m_matrix.zero();
  m_rhs.zero();

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
          m_matrix.add( element_index, element_index, element_matrix );
          m_rhs.add( element_index, element_rhs );
        }
      }
    });
  }
  m_matrix.close();
  m_rhs.close();

  if( verboseLevel() >= 2 )
  {
    string name = "matrix_" + std::to_string( time_n ) + ".mtx";
    m_matrix.write( name.c_str() );
    name = "rhs_" + std::to_string( time_n ) + ".mtx";
    m_rhs.write( name.c_str() );
  }
}

void LaplaceFEM::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                      real64 const scalingFactor,
                                      DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  dofManager.copyVectorToField( m_solution, m_fieldName, nodeManager );

  // Syncronize ghost nodes
  std::map<string, string_array> fieldNames;
  fieldNames["node"].push_back( m_fieldName );

  CommunicationTools::
  SynchronizeFields( fieldNames, mesh,
                     domain->getReference<array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );
}

void LaplaceFEM::ApplyBoundaryConditions( DomainPartition * const domain,
                                          systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                          real64 const time_n,
                                          real64 const dt )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ManagedGroup * const nodeManager = mesh->getNodeManager();
  FieldSpecificationManager * fsManager = FieldSpecificationManager::get();

  ApplyDirichletBC_implicit( time_n + dt, *domain, m_matrix, m_rhs );
  if( verboseLevel() >= 2 )
  {
    string name = "matrixDir_" + std::to_string( time_n+dt ) + ".mtx";
    m_matrix.write( name.c_str() );
    name = "rhsDir_" + std::to_string( time_n+dt ) + ".mtx";
    m_rhs.write( name.c_str() );
  }
}

void LaplaceFEM::SolveSystem( systemSolverInterface::EpetraBlockSystem * const blockSystem,
                              SystemSolverParameters const * const params )
{
  // Now create a solver from the parameter list
  LinearSolver solver( m_parameters );

  // Solve using the iterative solver and compare norms with true solution
  solver.solve( m_matrix, m_solution, m_rhs );

  if( verboseLevel() >= 2 )
  {
    string name = "sol.mtx";
    m_solution.write( name.c_str() );
  }
}

void LaplaceFEM::ApplyDirichletBC_implicit( real64 const time,
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
    bc->ApplyBoundaryConditionToSystem<FieldSpecificationEqual, TrilinosInterface>( targetSet,
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
