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

#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
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
  SetupSystem( domain, blockSystem );
}

void LaplaceFEM::ImplicitStepComplete( real64 const & time_n,
                                                             real64 const & dt,
                                                             DomainPartition * const domain)
{
}



void LaplaceFEM::SetNumRowsAndTrilinosIndices( ManagedGroup * const nodeManager,
                                               localIndex & numLocalRows,
                                               globalIndex & numGlobalRows,
                                               localIndex_array& localIndices,
                                               localIndex offset )
{
  int n_mpi_processes;
  MPI_Comm_size( MPI_COMM_GEOSX, &n_mpi_processes );

  int this_mpi_process = 0;
  MPI_Comm_rank( MPI_COMM_GEOSX, &this_mpi_process );

  std::vector<int> gather(n_mpi_processes);

  int intNumLocalRows = integer_conversion<int>(numLocalRows);
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
  globalIndex_array& trilinos_index = nodeManager->getReference<globalIndex_array>(laplaceFEMViewKeys.blockLocalDofNumber);
  integer_array const & is_ghost       = nodeManager->getReference<integer_array>(NodeManager::viewKeyStruct::ghostRankString);

  trilinos_index = -1;

  integer local_count = 0;
  for(integer r=0 ; r<trilinos_index.size() ; ++r )
  {
    if(is_ghost[r] < 0)
    {
      trilinos_index[r] = first_local_row+local_count+offset;
      local_count++;
    }
    else
    {
      trilinos_index[r] = -INT_MAX;
    }
  }

  GEOS_ASSERT(local_count == numLocalRows);
}


void LaplaceFEM :: SetupSystem ( DomainPartition * const domain,
                                                   EpetraBlockSystem * const blockSystem )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();

  localIndex dim = 1;//domain.m_feElementManager.m_ElementRegions.begin()->second.m_ElementDimension;
  localIndex n_ghost_rows  = nodeManager->GetNumberOfGhosts();
  localIndex n_local_rows  = nodeManager->size()-n_ghost_rows;
  globalIndex n_global_rows = 0;

  localIndex_array displacementIndices;
  SetNumRowsAndTrilinosIndices( nodeManager,
                                n_local_rows,
                                n_global_rows,
                                displacementIndices,
                                0 );

  std::map<string, string_array > fieldNames;
  fieldNames["node"].push_back(viewKeyStruct::blockLocalDofNumberString);

  CommunicationTools::SynchronizeFields(fieldNames,
                              mesh,
                              domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );



  // create epetra map

  Epetra_Map * const
  rowMap = blockSystem->SetRowMap( BlockIDs::dummyScalarBlock,
                                   std::make_unique<Epetra_Map>( dim*n_global_rows,
                                                                 dim*n_local_rows,
                                                                 0,
                                                                 m_linearSolverWrapper.m_epetraComm ) );

  Epetra_FECrsGraph * const
  sparsity = blockSystem->SetSparsity( BlockIDs::dummyScalarBlock,
                                       BlockIDs::dummyScalarBlock,
                                       std::make_unique<Epetra_FECrsGraph>(Copy,*rowMap,0) );


  SetSparsityPattern( domain, sparsity );

  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();

  blockSystem->SetMatrix( BlockIDs::dummyScalarBlock,
                          BlockIDs::dummyScalarBlock,
                          std::make_unique<Epetra_FECrsMatrix>(Copy,*sparsity) );

  blockSystem->SetSolutionVector( BlockIDs::dummyScalarBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

  blockSystem->SetResidualVector( BlockIDs::dummyScalarBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

}

void LaplaceFEM::SetSparsityPattern( DomainPartition const * const domain,
                                     Epetra_FECrsGraph * const sparsity )
{
  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ManagedGroup const * const nodeManager = mesh->getNodeManager();

  globalIndex_array const & trilinos_index = nodeManager->getReference<globalIndex_array>(laplaceFEMViewKeys.blockLocalDofNumber);
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  for( localIndex elemRegIndex=0 ; elemRegIndex<elemManager->numRegions() ; ++elemRegIndex )
  {
    ElementRegion const * const elementRegion = elemManager->GetRegion( elemRegIndex );

    elementRegion->forElementSubRegions([&]( auto const * const elementSubRegion )
    {
      localIndex const numElems = elementSubRegion->size();
      TYPEOFPTR(elementSubRegion)::NodeMapType const & elemsToNodes = elementSubRegion->nodeList();
      localIndex const numNodesPerElement = elemsToNodes.size(1);

      globalIndex_array elementLocalDofIndex (numNodesPerElement);

      array1d<integer> const & elemGhostRank = elementSubRegion->m_ghostRank;

      for( localIndex k=0 ; k<numElems ; ++k )
      {
        //if( elemGhostRank[k] < 0 )
        {
          for( localIndex a=0 ; a<numNodesPerElement ; ++a )
          {
            for(localIndex i=0 ; i<numNodesPerElement ; ++i)
            {
              elementLocalDofIndex[i] = trilinos_index[elemsToNodes[k][i]];
            }

            sparsity->InsertGlobalIndices(integer_conversion<int>(elementLocalDofIndex.size()),
                                          elementLocalDofIndex.data(),
                                          integer_conversion<int>(elementLocalDofIndex.size()),
                                          elementLocalDofIndex.data());
          }
        }
      }
    });
  }
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
  NumericalMethodsManager const * numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);
  FiniteElementDiscretizationManager const * feDiscretizationManager = numericalMethodManager->GetGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);


  Epetra_FECrsMatrix * const matrix = blockSystem->GetMatrix( BlockIDs::dummyScalarBlock,
                                                              BlockIDs::dummyScalarBlock );
  Epetra_FEVector * const rhs = blockSystem->GetResidualVector( BlockIDs::dummyScalarBlock );
  Epetra_FEVector * const solution = blockSystem->GetSolutionVector( BlockIDs::dummyScalarBlock );

  globalIndex_array const & trilinos_index = nodeManager->getReference<globalIndex_array>(laplaceFEMViewKeys.blockLocalDofNumber);

  for( auto & region : elemManager->GetGroup(dataRepository::keys::elementRegions)->GetSubGroups() )
  {
    ElementRegion * const elementRegion = ManagedGroup::group_cast<ElementRegion *>(region.second);

    FiniteElementDiscretization const * feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

    for( auto & elementSubRegion : elementRegion->GetGroup(ElementRegion::viewKeyStruct::elementSubRegions)->GetSubGroups() )
    {
      CellElementSubRegion * const cellElementSubRegion = ManagedGroup::group_cast<CellElementSubRegion*>(elementSubRegion.second );

      array3d< R1Tensor > & dNdX = cellElementSubRegion->getReference< array3d< R1Tensor > >(keys::dNdX);

      array2d<real64> const & detJ            = cellElementSubRegion->getReference< array2d<real64> >(keys::detJ);

      array2d<localIndex> const & elemsToNodes = cellElementSubRegion->nodeList();
      const integer numNodesPerElement = integer_conversion<int>(elemsToNodes.size(1));

      Epetra_LongLongSerialDenseVector  element_index   (numNodesPerElement);
      Epetra_SerialDenseVector     element_rhs     (numNodesPerElement);
      Epetra_SerialDenseMatrix     element_matrix  (numNodesPerElement,
                                                    numNodesPerElement);

      array1d<integer> const & elemGhostRank = cellElementSubRegion->m_ghostRank;
      const int n_q_points = feDiscretization->m_finiteElement->n_quadrature_points();

      // begin element loop, skipping ghost elements
      for( localIndex k=0 ; k<cellElementSubRegion->size() ; ++k )
      {
        if(elemGhostRank[k] < 0)
        {
          for( int a=0 ; a<numNodesPerElement ; ++a)
          {
            const localIndex n = elemsToNodes[k][a];
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
          }

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

  if( verboseLevel() >= 2 )
  {
    matrix->Print(std::cout);
    rhs->Print(std::cout);
  }

}


void LaplaceFEM::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                      real64 const scalingFactor,
                                      DomainPartition * const domain )
{
  NodeManager * const nodeManager = domain->getMeshBody(0)->getMeshLevel(0)->getNodeManager();
  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( BlockIDs::dummyScalarBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( BlockIDs::dummyScalarBlock );
  globalIndex_array const & trilinos_index = nodeManager->getReference<globalIndex_array>(laplaceFEMViewKeys.blockLocalDofNumber);

  string const & fieldName = getReference<string>(laplaceFEMViewKeys.fieldVarName);
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
void LaplaceFEM::ApplyBoundaryConditions( DomainPartition * const domain,
                                          systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                          real64 const time_n,
                                          real64 const dt )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ManagedGroup * const nodeManager = mesh->getNodeManager();
  FieldSpecificationManager * fsManager = FieldSpecificationManager::get();

  ApplyDirichletBC_implicit( time_n + dt, *domain, *blockSystem );

  if( verboseLevel() >= 2 )
  {
    Epetra_FECrsMatrix * const matrix = blockSystem->GetMatrix( BlockIDs::dummyScalarBlock,
                                                                BlockIDs::dummyScalarBlock );
    Epetra_FEVector * const rhs = blockSystem->GetResidualVector( BlockIDs::dummyScalarBlock );
    matrix->Print(std::cout);
    rhs->Print(std::cout);
  }

}

void LaplaceFEM::SolveSystem( systemSolverInterface::EpetraBlockSystem * const blockSystem,
                              SystemSolverParameters const * const params )
{
  SolverBase::SolveSystem( blockSystem, params, BlockIDs::dummyScalarBlock );
}

void LaplaceFEM::ApplyDirichletBC_implicit( real64 const time,
                                            DomainPartition & domain,
                                            EpetraBlockSystem & blockSystem )
{

  FieldSpecificationManager const * const fsManager = FieldSpecificationManager::get();

  fsManager->Apply( time,
                    &domain,
                    "nodeManager",
                    "Temperature",
                    [&]( FieldSpecificationBase const * const bc,
                    string const &,
                    set<localIndex> const & targetSet,
                    ManagedGroup * const targetGroup,
                    string const fieldName )->void
  {
    bc->ApplyBoundaryConditionToSystem<FieldSpecificationEqual>( targetSet,
                                                                 false,
                                                                 time,
                                                                 targetGroup,
                                                                 "Temperature",
                                                                 laplaceFEMViewKeys.blockLocalDofNumber.Key(),
                                                                 1,
                                                                 &blockSystem,
                                                                 BlockIDs::dummyScalarBlock );
  });
}



REGISTER_CATALOG_ENTRY( SolverBase, LaplaceFEM, std::string const &, ManagedGroup * const )
} /* namespace ANST */
