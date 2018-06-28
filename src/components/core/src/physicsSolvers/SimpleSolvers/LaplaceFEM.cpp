/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * NewtonianMechanics.cpp
 *
 *  Created on: Dec 4, 2014
 *      Author: rrsettgast
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
#include "constitutive/LinearElasticIsotropic.hpp"
#include "finiteElement/FiniteElementManager.hpp"
#include "finiteElement/FiniteElementSpaceManager.hpp"
#include "finiteElement/ElementLibrary/FiniteElement.h"
#include "finiteElement/Kinematics.h"
//#include "finiteElement/ElementLibrary/FiniteElementUtilities.h"
#include "physicsSolvers/BoundaryConditions/BoundaryConditionManager.hpp"

#include "codingUtilities/Utilities.hpp"

#include "managers/DomainPartition.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#define verbose 1


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
  this->RegisterGroup<SystemSolverParameters>( groupKeys.systemSolverParameters.Key() );
  m_linearSystem.SetBlockID( EpetraBlockSystem::BlockIDs::displacementBlock, this->getName() );
}



LaplaceFEM::~LaplaceFEM()
{
  // TODO Auto-generated destructor stub
}


void LaplaceFEM::FillDocumentationNode(  )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  SolverBase::FillDocumentationNode();

  docNode->setName(this->CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("An example solid mechanics solver");


  docNode->AllocateChildNode( viewKeys.timeIntegrationOption.Key(),
                              viewKeys.timeIntegrationOption.Key(),
                              -1,
                              "string",
                              "string",
                              "option for default time integration method",
                              "option for default time integration method",
                              "ExplicitDynamic",
                              "",
                              0,
                              1,
                              0 );


  docNode->AllocateChildNode( viewKeys.fieldVarName.Key(),
                              viewKeys.fieldVarName.Key(),
                              -1,
                              "string",
                              "string",
                              "name of field variable",
                              "name of field variable",
                              "Pressure",
                              "",
                              0,
                              1,
                              0 );
}

void LaplaceFEM::FillOtherDocumentationNodes( dataRepository::ManagedGroup * const rootGroup )
{
  DomainPartition * domain  = rootGroup->GetGroup<DomainPartition>(keys::domain);

  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    NodeManager * const nodes = ManagedGroup::group_cast<MeshBody*>(mesh.second)->getMeshLevel(0)->getNodeManager();
    cxx_utilities::DocumentationNode * const docNode = nodes->getDocumentationNode();


    docNode->AllocateChildNode( viewKeys.trilinosIndex.Key(),
                                viewKeys.trilinosIndex.Key(),
                                -1,
                                "localIndex_array",
                                "localIndex_array",
                                "",
                                "",
                                "",
                                keys::nodeManager,
                                1,
                                0,
                                0 );

    docNode->AllocateChildNode( viewKeys.ghostRank.Key(),
                                viewKeys.ghostRank.Key(),
                                -1,
                                "integer_array",
                                "integer_array",
                                "",
                                "",
                                "",
                                keys::nodeManager,
                                1,
                                0,
                                0 );

    docNode->AllocateChildNode( "Temperature",
                                "Temperature",
                                -1,
                                "real64_array",
                                "real64_array",
                                "",
                                "",
                                "",
                                keys::nodeManager,
                                1,
                                0,
                                0 );
  }
}


void LaplaceFEM::ReadXML_PostProcess()
{
  string tiOption = this->getReference<string>(viewKeys.timeIntegrationOption);

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

//void LaplaceFEM::BuildDataStructure( ManagedGroup * const domain )
//{
//  SolverBase::BuildDataStructure( domain );
//
//  // Test auto-registration:
//  RegisterDocumentationNodes();
//
//}


void LaplaceFEM::InitializePreSubGroups( ManagedGroup * const problemManager )
{

}

void LaplaceFEM::TimeStep( real64 const& time_n,
                                             real64 const& dt,
                                             const int cycleNumber,
                                             ManagedGroup * domain )
{

  if( m_timeIntegrationOption == timeIntegrationOption::ExplicitTransient )
  {
    TimeStepExplicit( time_n, dt, cycleNumber, ManagedGroup::group_cast<DomainPartition*>(domain) );
  }
  else if( m_timeIntegrationOption == timeIntegrationOption::ImplicitTransient ||
           m_timeIntegrationOption == timeIntegrationOption::SteadyState )
  {
    this->TimeStepImplicit( time_n, dt, cycleNumber, domain->group_cast<DomainPartition*>() );
  }
}

void LaplaceFEM::TimeStepExplicit( real64 const& time_n,
                                                     real64 const& dt,
                                                     const int cycleNumber,
                                                     DomainPartition * const domain )
{


}


void LaplaceFEM::ApplyDirichletBC_implicit( ManagedGroup * object,
                                            BoundaryConditionBase const * const bc,
                                            lSet const & set,
                                            real64 const time_n,
                                            EpetraBlockSystem & blockSystem )
{
  bc->ApplyDirichletBounaryConditionDefaultMethod<0>( set,
                                                      time_n,
                                                      object,
                                                      "Temperature",
                                                      viewKeys.trilinosIndex.Key(),
                                                      1,
                                                      &blockSystem,
                                                      EpetraBlockSystem::BlockIDs::displacementBlock );
}


void LaplaceFEM::TimeStepImplicitSetup( real64 const& time_n,
                                                          real64 const& dt,
                                                          DomainPartition * const domain )
{
}

void LaplaceFEM::TimeStepImplicitComplete( real64 const & time_n,
                                                             real64 const & dt,
                                                             DomainPartition * const domain)
{
}


real64 LaplaceFEM::TimeStepImplicit( real64 const & time_n,
                                                       real64 const & dt,
                                                       integer const cycleNumber,
                                                       DomainPartition * const domain )
{
  real64 dt_return = dt;
////  view_rtype<r1_array> uhat  =
// domain.m_feNodeManager.getData<r1_array>(keys::IncrementalDisplacement);

#if 1
  TimeStepImplicitSetup( time_n, dt, domain );


  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


      SetupSystem( domain, &m_linearSystem );

      Epetra_FECrsMatrix * matrix = m_linearSystem.GetMatrix( EpetraBlockSystem::BlockIDs::displacementBlock,
                                                              EpetraBlockSystem::BlockIDs::displacementBlock );
      matrix->Scale(0.0);

      Epetra_FEVector * rhs = m_linearSystem.GetResidualVector( EpetraBlockSystem::BlockIDs::displacementBlock );
      rhs->Scale(0.0);

      Epetra_FEVector * solution = m_linearSystem.GetSolutionVector( EpetraBlockSystem::BlockIDs::displacementBlock );
      solution->Scale(0.0);

      Assemble( domain, &m_linearSystem, time_n+dt, dt );

//      matrix->Print(std::cout);
//      rhs->Print(std::cout);

      MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
      ManagedGroup * const nodeManager = mesh->getNodeManager();

      BoundaryConditionManager * bcManager = BoundaryConditionManager::get();

      bcManager->ApplyBoundaryCondition( this,
                                         &LaplaceFEM::ApplyDirichletBC_implicit,
                                         nodeManager,
                                         "Temperature",
                                         time_n + dt,
                                         m_linearSystem );

      matrix->Print(std::cout);
      rhs->Print(std::cout);

      rhs->Scale(-1.0);

      if(verbose)
        std::cout<<"Solving system"<<std::endl;




      this->m_linearSolverWrapper.SolveSingleBlockSystem( &m_linearSystem,
                                                          getSystemSolverParameters(),
                                                          systemSolverInterface::EpetraBlockSystem::BlockIDs::displacementBlock );

      solution->Print(std::cout);

      ApplySystemSolution( &m_linearSystem, 1.0, 0, nodeManager );

  TimeStepImplicitComplete( time_n, dt,  domain );
#endif

  return dt_return;
}


void LaplaceFEM::SetNumRowsAndTrilinosIndices( ManagedGroup * const nodeManager,
                                               localIndex & numLocalRows,
                                               localIndex & numGlobalRows,
                                               localIndex_array& localIndices,
                                               localIndex offset )
{
//  dim =
// domain.m_feElementManager.m_ElementRegions.begin()->second.m_ElementDimension;
  int dim = 1;


  int n_mpi_processes;
  MPI_Comm_size( MPI_COMM_WORLD, &n_mpi_processes );

  int this_mpi_process = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &this_mpi_process );

  std::vector<int> gather(n_mpi_processes);

  int intNumLocalRows = static_cast<int>(numLocalRows);
  m_linearSolverWrapper.m_epetraComm.GatherAll( &intNumLocalRows,
                                                &gather.front(),
                                                1 );
  numLocalRows = intNumLocalRows;

  localIndex first_local_row = 0;
  numGlobalRows = 0;

  for( integer p=0 ; p<n_mpi_processes ; ++p)
  {
    numGlobalRows += gather[p];
    if(p<this_mpi_process)
      first_local_row += gather[p];
  }

  // create trilinos dof indexing

  localIndex_array& trilinos_index = nodeManager->getReference<localIndex_array>(viewKeys.trilinosIndex);
  integer_array const & is_ghost       = nodeManager->getReference<integer_array>(viewKeys.ghostRank);


  trilinos_index = -1;

  integer local_count = 0;
  for(integer r=0 ; r<trilinos_index.size() ; ++r )
  {
    if(is_ghost[r] < 0)
    {
      trilinos_index[r] = first_local_row+local_count+offset;
      localIndices.push_back(trilinos_index[r]);
      local_count++;
    }
    else
    {
      trilinos_index[r] = -INT_MAX;
    }
  }

  assert(local_count == numLocalRows );


}


void LaplaceFEM :: SetupSystem ( DomainPartition * const domain,
                                                   EpetraBlockSystem * const blockSystem )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();

  localIndex dim = 1;//domain.m_feElementManager.m_ElementRegions.begin()->second.m_ElementDimension;
  localIndex n_ghost_rows  = nodeManager->GetNumberOfGhosts();
  localIndex n_local_rows  = nodeManager->size()-n_ghost_rows;
  localIndex n_global_rows = 0;

  localIndex_array displacementIndices;
  SetNumRowsAndTrilinosIndices( nodeManager,
                                n_local_rows,
                                n_global_rows,
                                displacementIndices,
                                0 );

  std::map<string, array<string> > fieldNames;
  fieldNames["node"].push_back("trilinosIndex_LaplaceFEM");

  CommunicationTools::SynchronizeFields(fieldNames,
                              mesh,
                              domain->getReference< array<NeighborCommunicator> >( domain->viewKeys.neighbors ) );



  // create epetra map

  Epetra_Map * junk =  new Epetra_Map( static_cast<int>(dim*n_global_rows),
                                       static_cast<int>(dim*n_local_rows),
                                       0,
                                       m_linearSolverWrapper.m_epetraComm );

  Epetra_Map * const rowMap = blockSystem->SetRowMap( EpetraBlockSystem::BlockIDs::displacementBlock,
                                                      std::make_unique<Epetra_Map>( static_cast<int>(dim*n_global_rows),
                                                                                    static_cast<int>(dim*n_local_rows),
                                                                                    0,
                                                                                    m_linearSolverWrapper.m_epetraComm ) );

  Epetra_FECrsGraph * const sparsity = blockSystem->SetSparsity( EpetraBlockSystem::BlockIDs::displacementBlock,
                                                                 EpetraBlockSystem::BlockIDs::displacementBlock,
                                                                 std::make_unique<Epetra_FECrsGraph>(Copy,*rowMap,0) );



  SetSparsityPattern( domain, sparsity );

  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();

  blockSystem->SetMatrix( EpetraBlockSystem::BlockIDs::displacementBlock,
                          EpetraBlockSystem::BlockIDs::displacementBlock,
                          std::make_unique<Epetra_FECrsMatrix>(Copy,*sparsity) );

  blockSystem->SetSolutionVector( EpetraBlockSystem::BlockIDs::displacementBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

  blockSystem->SetResidualVector( EpetraBlockSystem::BlockIDs::displacementBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

}

void LaplaceFEM::SetSparsityPattern( DomainPartition const * const domain,
                                                       Epetra_FECrsGraph * const sparsity )
{
  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ManagedGroup const * const nodeManager = mesh->getNodeManager();

  localIndex_array const & trilinos_index = nodeManager->getReference<localIndex_array>(viewKeys.trilinosIndex);
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  for( localIndex elemRegIndex=0 ; elemRegIndex<elemManager->numRegions() ; ++elemRegIndex )
  {
    ElementRegion const * const elementRegion = elemManager->GetRegion( elemRegIndex );
      auto const & numMethodName = elementRegion->getData<string>(keys::numericalMethod);

      for( localIndex subRegionIndex=0 ; subRegionIndex<elementRegion->numSubRegions() ; ++subRegionIndex )
      {
        CellBlockSubRegion const * const cellBlock = elementRegion->GetSubRegion(subRegionIndex);
        localIndex const numElems = cellBlock->size();
        lArray2d const & elemsToNodes = cellBlock->getWrapper<FixedOneToManyRelation>(cellBlock->viewKeys().nodeList)->reference();// getData<lArray2d>(keys::nodeList);
        localIndex const numNodesPerElement = elemsToNodes.size(1);

        integer_array elementLocalDofIndex (numNodesPerElement);

        array<integer> const & elemGhostRank = cellBlock->m_ghostRank;

        for( localIndex k=0 ; k<numElems ; ++k )
        {
          if( elemGhostRank[k] < 0 )
          {
            arrayView1d<localIndex const> const localNodeIndices = elemsToNodes[k];

            for( localIndex a=0 ; a<numNodesPerElement ; ++a )
            {
              for(localIndex i=0 ; i<numNodesPerElement ; ++i)
              {
                elementLocalDofIndex[i] = integer_conversion<int>(trilinos_index[localNodeIndices[i]]);
              }

              sparsity->InsertGlobalIndices(integer_conversion<int>(elementLocalDofIndex.size()),
                                            elementLocalDofIndex.data(),
                                            integer_conversion<int>(elementLocalDofIndex.size()),
                                            elementLocalDofIndex.data());
            }
          }

        }
      }
    }
}



real64 LaplaceFEM::Assemble ( DomainPartition * const  domain,
                                                EpetraBlockSystem * const blockSystem,
                                                real64 const time_n,
                                                real64 const dt )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ManagedGroup * const nodeManager = mesh->getNodeManager();
  ConstitutiveManager  * const constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  FiniteElementManager const * numericalMethodManager = domain->getParent()->GetGroup<FiniteElementManager>(keys::finiteElementManager);
  FiniteElementSpaceManager const * feSpaceManager = numericalMethodManager->GetGroup<FiniteElementSpaceManager>(keys::finiteElementSpaces);


  Epetra_FECrsMatrix * const matrix = blockSystem->GetMatrix( EpetraBlockSystem::BlockIDs::displacementBlock,
                                                              EpetraBlockSystem::BlockIDs::displacementBlock );
  Epetra_FEVector * const rhs = blockSystem->GetResidualVector( EpetraBlockSystem::BlockIDs::displacementBlock );
  Epetra_FEVector * const solution = blockSystem->GetSolutionVector( EpetraBlockSystem::BlockIDs::displacementBlock );

  localIndex_array const & trilinos_index = nodeManager->getReference<localIndex_array>(viewKeys.trilinosIndex);

  for( auto & region : elemManager->GetGroup(dataRepository::keys::elementRegions)->GetSubGroups() )
  {
    ElementRegion * const elementRegion = ManagedGroup::group_cast<ElementRegion *>(region.second);
    auto const & numMethodName = elementRegion->getData<string>(keys::numericalMethod);
    FiniteElementSpace const * feSpace = feSpaceManager->GetGroup<FiniteElementSpace>(numMethodName);

    for( auto & cellBlock : elementRegion->GetGroup(dataRepository::keys::cellBlockSubRegions)->GetSubGroups() )
    {
      CellBlockSubRegion * const cellBlockSubRegion = ManagedGroup::group_cast<CellBlockSubRegion*>(cellBlock.second );

      multidimensionalArray::ManagedArray< R1Tensor, 3 > & dNdX = cellBlockSubRegion->getReference< multidimensionalArray::ManagedArray< R1Tensor, 3 > >(keys::dNdX);

      Array2dT<real64> const & detJ            = cellBlockSubRegion->getReference< Array2dT<real64> >(keys::detJ);

      lArray2d const & elemsToNodes = cellBlockSubRegion->getWrapper<FixedOneToManyRelation>(cellBlockSubRegion->viewKeys().nodeList)->reference();
      const integer numNodesPerElement = integer_conversion<int>(elemsToNodes.size(1));

      Epetra_IntSerialDenseVector  element_index   (numNodesPerElement);
      Epetra_SerialDenseVector     element_rhs     (numNodesPerElement);
      Epetra_SerialDenseMatrix     element_matrix  (numNodesPerElement,
                                                    numNodesPerElement);

      array<integer> const & elemGhostRank = cellBlockSubRegion->m_ghostRank;
      const int n_q_points = feSpace->m_finiteElement->n_quadrature_points();

      // begin element loop, skipping ghost elements
      for( localIndex k=0 ; k<cellBlockSubRegion->size() ; ++k )
      {
        if(elemGhostRank[k] < 0)
        {
          arrayView1d<localIndex const> const local_index = elemsToNodes[k];

          for( int a=0 ; a<numNodesPerElement ; ++a)
          {
            const localIndex n = local_index[a];
            element_index[a] = integer_conversion<int>(trilinos_index[n]);
          }

          element_rhs.Scale(0);
          element_matrix.Scale(0);

          for( int q=0 ; q<n_q_points ; ++q)
          {
            for( int a=0 ; a<numNodesPerElement ; ++a)
            {
//              element_rhs(a) += detJ[k][q] *
//                                equation_data.source *
//                                fe.value(a,q);

              double diffusion = 1;
              for( int b=0 ; b<numNodesPerElement ; ++b)
              {
                element_matrix(a,b) += detJ[k][q] *
                                       diffusion *
                                     + Dot( dNdX[k][q][a], dNdX[k][q][b] );
              }

            }
//            element_matrix.Print(std::cout);
          }

//          element_matrix.Print(std::cout);

          matrix->SumIntoGlobalValues( element_index,
                                       element_matrix);


          rhs->SumIntoGlobalValues( element_index,
                                    element_rhs);

        }
      }
    }
  }

  matrix->GlobalAssemble(true);
  rhs->GlobalAssemble();

  return 0;
}


void LaplaceFEM::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                      real64 const scalingFactor,
                                      localIndex const dofOffset,
                                      dataRepository::ManagedGroup * const nodeManager )
{
  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( EpetraBlockSystem::BlockIDs::displacementBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( EpetraBlockSystem::BlockIDs::displacementBlock );
  localIndex_array const & trilinos_index = nodeManager->getReference<localIndex_array>(viewKeys.trilinosIndex);

  string const & fieldName = getReference<string>(viewKeys.fieldVarName);
  real64_array & fieldVar = nodeManager->getReference<real64_array>(string("Temperature"));

//  integer_array & ghostRank = nodeManager->getReference<integer_array>(NodeManager::viewKeyStruct::ghostRankString);

  int dummy;
  double* local_solution = nullptr;
  solution->ExtractView(&local_solution,&dummy);

  for( localIndex r=0 ; r<fieldVar.size() ; ++r)
  {
//    if(ghostRank[r] < 0)
    {
      int lid = rowMap->LID(integer_conversion<int>(trilinos_index[r]));
      fieldVar[r] = local_solution[lid];
    }
  }


}


REGISTER_CATALOG_ENTRY( SolverBase, LaplaceFEM, std::string const &, ManagedGroup * const )
} /* namespace ANST */
