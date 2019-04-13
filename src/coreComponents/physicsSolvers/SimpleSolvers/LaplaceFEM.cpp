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
#include "constitutive/LinearElasticIsotropic.hpp"
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
  // Just one computation of the sparsity pattern!!!
  if( firstTime )
  {
    firstTime = false;
    SetupSystem( domain, blockSystem );
  }
}

void LaplaceFEM::ImplicitStepComplete( real64 const & time_n,
                                       real64 const & dt,
                                       DomainPartition * const domain)
{
}

void LaplaceFEM::SetupSystem( DomainPartition * const domain,
                              EpetraBlockSystem * const blockSystem )
{
  dofManager.setMesh( domain, 0, 0 );
  dofManager.addField( "Temperature", DofManager::Location::Node, DofManager::Connectivity::Elem );

  ParallelMatrix & sparsity = m_matrix;
  dofManager.getSparsityPattern( sparsity, "Temperature", "Temperature", true, &m_rhs );

  /////////////////////////////////////////////////////////////////////////////////////
  // TRILINOS CODE
  Epetra_FECrsMatrix const * const sparsityPtr = sparsity.unwrappedPointer();

  Epetra_Map const * const rowMap0 = &(sparsityPtr->RowMap());

  // create epetra map
  Epetra_Map * const
  rowMap = blockSystem->SetRowMap( BlockIDs::dummyScalarBlock,
                                   std::make_unique<Epetra_Map>( *rowMap0 ));

  blockSystem->SetMatrix( BlockIDs::dummyScalarBlock,
                          BlockIDs::dummyScalarBlock,
                          std::make_unique<Epetra_FECrsMatrix>(*sparsityPtr) );

  blockSystem->SetSolutionVector( BlockIDs::dummyScalarBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

  blockSystem->SetResidualVector( BlockIDs::dummyScalarBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );
  /////////////////////////////////////////////////////////////////////////////////////
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

  globalIndex_array const & indexArray = nodeManager->getReference<globalIndex_array>( "Temperature_dof_indices" );

  {
  // Initialize all entries to zero
  m_matrix.scale( 0.0 );
  m_rhs.scale( 0.0 );

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
          dofManager.getIndices( element_index, DofManager::Connectivity::Elem, er, esr, k, "Temperature" );

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
    string name = "NEWmatrix_" + std::to_string( time_n ) + ".mtx";
    m_matrix.write( name.c_str() );
    name = "NEWrhs_" + std::to_string( time_n ) + ".mtx";
    EpetraExt::MultiVectorToMatrixMarketFile( name.c_str(), *m_rhs.unwrappedPointer() );
  }
  }

  /////////////////////////////////////////////////////////////////////////////////////
  // TRILINOS CODE
  Epetra_FECrsMatrix * const matrix = blockSystem->GetMatrix( BlockIDs::dummyScalarBlock,
                                                              BlockIDs::dummyScalarBlock );
  Epetra_FEVector * const rhs = blockSystem->GetResidualVector( BlockIDs::dummyScalarBlock );
  Epetra_FEVector * const solution = blockSystem->GetSolutionVector( BlockIDs::dummyScalarBlock );

  matrix->Scale( 0.0 );

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

      Epetra_LongLongSerialDenseVector element_index(numNodesPerElement);
      Epetra_SerialDenseVector element_rhs(numNodesPerElement);
      Epetra_SerialDenseMatrix element_matrix(numNodesPerElement, numNodesPerElement);

      array1d<integer> const & elemGhostRank = elementSubRegion->m_ghostRank;
      const int n_q_points = feDiscretization->m_finiteElement->n_quadrature_points();

      // begin element loop, skipping ghost elements
      for( localIndex k=0 ; k<elementSubRegion->size() ; ++k )
      {
        if(elemGhostRank[k] < 0)
        {
          globalIndex_array indices;
          dofManager.getIndices( indices, DofManager::Connectivity::Elem, er, esr, k, "Temperature" );
          for( int i = 0 ; i < indices.size() ; ++i )
          {
            element_index[i] = integer_conversion<long long>( indices[i] );
          }

          element_rhs.Scale(0);
          element_matrix.Scale(0);

          for( int q=0 ; q<n_q_points ; ++q)
          {
            for( int a=0 ; a<numNodesPerElement ; ++a)
            {
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
    });
  }

  matrix->GlobalAssemble(true);
  rhs->GlobalAssemble();

  if( verboseLevel() >= 2 )
  {
    string name = "matrix_" + std::to_string( time_n ) + ".mtx";
    EpetraExt::RowMatrixToMatrixMarketFile( name.c_str(), *matrix );
    name = "rhs_" + std::to_string( time_n ) + ".mtx";
    EpetraExt::MultiVectorToMatrixMarketFile( name.c_str(), *rhs );
  }
  /////////////////////////////////////////////////////////////////////////////////////
}

void LaplaceFEM::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                      real64 const scalingFactor,
                                      DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  // Retrieve fieldVar
  real64_array & fieldVar = nodeManager->getReference<real64_array>(string("Temperature"));

  /////////////////////////////////////////////////////////////////////////////////////
  // TRILINOS CODE
  {
  // Retrieve original map and solution
  Epetra_Map const * const rowMap = blockSystem->GetRowMap( BlockIDs::dummyScalarBlock );
  Epetra_FEVector const * solution = blockSystem->GetSolutionVector( BlockIDs::dummyScalarBlock );

  int dummy;
  double* local_solution = nullptr;
  solution->ExtractView(&local_solution,&dummy);

  globalIndex_array const & indexArray = nodeManager->getReference<globalIndex_array>( "Temperature_dof_indices" );

  // Map values from local_solution to fieldVar
  for( localIndex r = 0 ; r < indexArray.size() ; ++r )
  {
    int lid = rowMap->LID(integer_conversion<int>(indexArray[r]));
    // Check if it is available
    if( lid >= 0 )
    {
      fieldVar[r] = local_solution[lid];
    }
  }
  }
  /////////////////////////////////////////////////////////////////////////////////////

  {
  int dummy;
  double* local_solution = nullptr;
  m_solution.unwrappedPointer()->ExtractView(&local_solution,&dummy);

  globalIndex_array const & indexArray = nodeManager->getReference<globalIndex_array>( "Temperature_dof_indices" );

  // Map values from local_solution to fieldVar
  for( localIndex r = 0 ; r < indexArray.size() ; ++r )
  {
    localIndex lid = m_matrix.ParallelMatrixGetLocalRowID( indexArray[r] );
    // Check if it is available
    if( lid >= 0 )
    {
      fieldVar[r] = local_solution[lid];
    }
  }
  }

  // Syncronize ghost nodes
  std::map<string, string_array> fieldNames;
  fieldNames["node"].push_back( "Temperature" );

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
    string name = "NEWmatrixDir_" + std::to_string( time_n ) + ".mtx";
    m_matrix.write( name.c_str() );
    name = "NEWrhsDir_" + std::to_string( time_n ) + ".mtx";
    EpetraExt::MultiVectorToMatrixMarketFile( name.c_str(), *m_rhs.unwrappedPointer() );
  }

  ///////////////////////////////////////////////////////////////////////////////
  // TRILINOS CODE
  ApplyDirichletBC_implicit( time_n + dt, *domain, *blockSystem );

  if( verboseLevel() >= 2 )
  {
    Epetra_FECrsMatrix * const matrix = blockSystem->GetMatrix( BlockIDs::dummyScalarBlock,
                                                                BlockIDs::dummyScalarBlock );
    Epetra_FEVector * const rhs = blockSystem->GetResidualVector( BlockIDs::dummyScalarBlock );

    string name = "matrixDir_" + std::to_string( time_n ) + ".mtx";
    EpetraExt::RowMatrixToMatrixMarketFile( name.c_str(), *matrix );
    name = "rhsDir_" + std::to_string( time_n ) + ".mtx";
    EpetraExt::MultiVectorToMatrixMarketFile( name.c_str(), *rhs );
  }
  ///////////////////////////////////////////////////////////////////////////////
}

void LaplaceFEM::SolveSystem( systemSolverInterface::EpetraBlockSystem * const blockSystem,
                              SystemSolverParameters const * const params )
{
  m_solution.create( m_rhs );
  solve( &m_matrix, &m_rhs, &m_solution, params );

  if( verboseLevel() >= 2 )
  {
    string name = "NEWsol.mtx";
    EpetraExt::MultiVectorToMatrixMarketFile( name.c_str(), *m_solution.unwrappedPointer() );
  }

  ///////////////////////////////////////////////////////////////////////////////
  // TRILINOS CODE
  SolverBase::SolveSystem( blockSystem, params, BlockIDs::dummyScalarBlock );

  if( verboseLevel() >= 2 )
  {
    Epetra_FEVector const * solution = blockSystem->GetSolutionVector( BlockIDs::dummyScalarBlock );
    string name = "sol.mtx";
    EpetraExt::MultiVectorToMatrixMarketFile( name.c_str(), *solution );
  }
  ///////////////////////////////////////////////////////////////////////////////
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
                                                                 time,
                                                                 targetGroup,
                                                                 "Temperature",
                                                                 "Temperature_dof_indices",
                                                                 1,
                                                                 &blockSystem,
                                                                 BlockIDs::dummyScalarBlock );
  });
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
                    "Temperature",
                    [&]( FieldSpecificationBase const * const bc,
                    string const &,
                    set<localIndex> const & targetSet,
                    ManagedGroup * const targetGroup,
                    string const fieldName )->void
  {
    integer const component = bc->GetComponent();
    real64 const scale = bc->GetScale();
    string const & functionName = bc->getReference<string>( FieldSpecificationBase::viewKeyStruct::functionNameString );
    ApplyBoundaryConditionToSystem( *fsManager,
                                    functionName,
                                    targetSet,
                                    time,
                                    targetGroup,
                                    "Temperature",
                                    "Temperature_dof_indices",
                                    1,
                                    component,
                                    scale,
                                    matrix,
                                    rhs );
  });
}

/////////////////////////////////////////////////////////////////////////////////////////////
void LaplaceFEM::ApplyBoundaryConditionToSystem( FieldSpecificationManager const & fsManager,
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
                                                 ParallelVector & rhs )
{
  dataRepository::ViewWrapperBase * vw = dataGroup->getWrapperBase( fieldName );
  std::type_index typeIndex = std::type_index( vw->get_typeid());
  arrayView1d<globalIndex> const & dofMap = dataGroup->getReference<array1d<globalIndex>>( dofMapName );

  rtTypes::ApplyArrayTypeLambda1( rtTypes::typeID( typeIndex ),
    [&]( auto type ) -> void
    {
      using fieldType = decltype(type);
      dataRepository::ViewWrapper<fieldType> & view = dynamic_cast< dataRepository::ViewWrapper<fieldType> & >(*vw);
      fieldType & field = view.reference();

      ApplyBoundaryConditionToSystem( fsManager, functionName, targetSet, time, dataGroup, dofMap, dofDim, component, scale, matrix, rhs,
        [&]( localIndex const a )->real64
        {
          return static_cast<real64>(rtTypes::value( field[a], component ));
        }
      );
    }
  );
}

/*
template< typename LAMBDA >
void LaplaceFEM::ApplyBoundaryConditionToSystem( FieldSpecificationManager const & fsManager,
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
*/

void LaplaceFEM::solve( ParallelMatrix const * const matrix,
                        ParallelVector * const rhs,
                        ParallelVector * const solution,
                        SystemSolverParameters const * const params )
{
  solution->set( 0.0 );


  //if(params->scalingOption())
  if(true)
  {
    Epetra_Vector scaling(matrix->unwrappedPointer()->RowMap());
    matrix->unwrappedPointer()->InvRowSums(scaling);
    matrix->unwrappedPointer()->LeftScale(scaling);

    Epetra_MultiVector tmp(*rhs->unwrappedPointer());
    rhs->unwrappedPointer()->Multiply(1.0,scaling,tmp,0.0);
  }

  Epetra_LinearProblem problem( matrix->unwrappedPointer(),
                                solution->unwrappedPointer(),
                                rhs->unwrappedPointer() );
  AztecOO solver(problem);
  solver.SetAztecOption(AZ_precond,AZ_dom_decomp);
  solver.SetAztecOption(AZ_subdomain_solve,AZ_ilut);
  solver.SetAztecParam(AZ_ilut_fill,params->ilut_fill());
  solver.SetAztecParam(AZ_drop,params->ilut_drop());
  solver.SetAztecOption(AZ_output,params->verbose());
  solver.Iterate(params->numKrylovIter(),
                 params->krylovTol() );
}

/////////////////////////////////////////////////////////////////////////////////////////////

REGISTER_CATALOG_ENTRY( SolverBase, LaplaceFEM, std::string const &, ManagedGroup * const )
} /* namespace ANST */
