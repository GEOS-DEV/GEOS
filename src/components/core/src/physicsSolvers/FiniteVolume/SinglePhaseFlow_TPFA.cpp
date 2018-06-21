/*
 * NewtonianMechanics.cpp
 *
 *  Created on: Dec 4, 2014
 *      Author: rrsettgast
 */

#include "SinglePhaseFlow_TPFA.hpp"

#include <vector>
#include <math.h>

#include "common/TimingMacros.hpp"

#include "dataRepository/ManagedGroup.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "finiteElement/FiniteElementManager.hpp"
#include "finiteElement/FiniteElementSpaceManager.hpp"
#include "finiteElement/ElementLibrary/FiniteElement.h"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "physicsSolvers/BoundaryConditions/BoundaryConditionManager.hpp"




#define verbose 1

#define local_dim 3 //number of local dimensions

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace systemSolverInterface;


SinglePhaseFlow_TPFA::SinglePhaseFlow_TPFA( const std::string& name,
                                                            ManagedGroup * const parent ):
  SolverBase( name, parent )
{
  this->RegisterGroup<SystemSolverParameters>( groupKeys.systemSolverParameters.Key() );
  m_linearSystem.SetBlockID( EpetraBlockSystem::BlockIDs::fluidPressureBlock, this->getName() );
}



SinglePhaseFlow_TPFA::~SinglePhaseFlow_TPFA()
{}


void SinglePhaseFlow_TPFA::FillDocumentationNode(  )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  SolverBase::FillDocumentationNode();

  docNode->setName(this->CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("An example solid mechanics solver");


  docNode->AllocateChildNode( viewKeys.functionalSpace.Key(),
                              viewKeys.functionalSpace.Key(),
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

void SinglePhaseFlow_TPFA::FillOtherDocumentationNodes( dataRepository::ManagedGroup * const rootGroup )
{
  DomainPartition * domain  = rootGroup->GetGroup<DomainPartition>(keys::domain);

  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel * meshLevel = ManagedGroup::group_cast<MeshBody*>(mesh.second)->getMeshLevel(0);
    ElementRegionManager * const elemManager = meshLevel->getElemManager();
//    cxx_utilities::DocumentationNode * const docNode = elemManager->getDocumentationNode();

    elemManager->forCellBlocks( [&]( CellBlockSubRegion * const cellBlock )->void
    {
      cxx_utilities::DocumentationNode * const docNode = cellBlock->getDocumentationNode();
      docNode->AllocateChildNode( viewKeys.trilinosIndex.Key(),
                                  viewKeys.trilinosIndex.Key(),
                                  -1,
                                  "localIndex_array",
                                  "localIndex_array",
                                  "",
                                  "",
                                  "",
                                  elemManager->getName(),
                                  1,
                                  0,
                                  0 );

      docNode->AllocateChildNode( "ghostRank",
                                  viewKeys.ghostRank.Key(),
                                  -1,
                                  "integer_array",
                                  "integer_array",
                                  "",
                                  "",
                                  "",
                                  elemManager->getName(),
                                  1,
                                  0,
                                  0 );

      docNode->AllocateChildNode( "Pressure",
                                  "Pressure",
                                  -1,
                                  "real64_array",
                                  "real64_array",
                                  "",
                                  "",
                                  "",
                                  elemManager->getName(),
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
                                  elemManager->getName(),
                                  1,
                                  0,
                                  0 );

    });
  }
}


//void SinglePhaseFlow_TPFA::ReadXML_PostProcess()
//{
//
//}

void SinglePhaseFlow_TPFA::InitializePreSubGroups( ManagedGroup * const problemManager )
{

}

void SinglePhaseFlow_TPFA::TimeStep( real64 const& time_n,
                                             real64 const& dt,
                                             const int cycleNumber,
                                             ManagedGroup * domain )
{
  this->TimeStepImplicit( time_n, dt, cycleNumber, domain->group_cast<DomainPartition*>() );
}


void SinglePhaseFlow_TPFA::ApplyDirichletBC_implicit( ManagedGroup * object,
                                            BoundaryConditionBase const * const bc,
                                            lSet const & set,
                                            real64 const time_n,
                                            EpetraBlockSystem & blockSystem )
{
//  bc->ApplyDirichletBounaryConditionDefaultMethod<0>( set,
//                                                      time_n,
//                                                      object,
//                                                      "Temperature",
//                                                      viewKeys.trilinosIndex.Key(),
//                                                      1,
//                                                      &blockSystem,
//                                                      EpetraBlockSystem::BlockIDs::fluidPressureBlock );
}


void SinglePhaseFlow_TPFA::TimeStepImplicitSetup( real64 const& time_n,
                                                          real64 const& dt,
                                                          DomainPartition * const domain )
{
}

void SinglePhaseFlow_TPFA::TimeStepImplicitComplete( real64 const & time_n,
                                                             real64 const & dt,
                                                             DomainPartition * const domain)
{
}


real64 SinglePhaseFlow_TPFA::TimeStepImplicit( real64 const & time_n,
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

      Epetra_FECrsMatrix * matrix = m_linearSystem.GetMatrix( EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                                                              EpetraBlockSystem::BlockIDs::fluidPressureBlock );
      matrix->Scale(0.0);

      Epetra_FEVector * rhs = m_linearSystem.GetResidualVector( EpetraBlockSystem::BlockIDs::fluidPressureBlock );
      rhs->Scale(0.0);

      Epetra_FEVector * solution = m_linearSystem.GetSolutionVector( EpetraBlockSystem::BlockIDs::fluidPressureBlock );
      solution->Scale(0.0);

      Assemble( domain, &m_linearSystem, time_n+dt, dt );

//      matrix->Print(std::cout);
//      rhs->Print(std::cout);

      MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
      ManagedGroup * const nodeManager = mesh->getNodeManager();

      BoundaryConditionManager * bcManager = BoundaryConditionManager::get();

      bcManager->ApplyBoundaryCondition( this,
                                         &SinglePhaseFlow_TPFA::ApplyDirichletBC_implicit,
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
                                                          systemSolverInterface::EpetraBlockSystem::BlockIDs::fluidPressureBlock );

      solution->Print(std::cout);

      ApplySystemSolution( &m_linearSystem, 1.0, 0, nodeManager );

  TimeStepImplicitComplete( time_n, dt,  domain );
#endif

  return dt_return;
}


void SinglePhaseFlow_TPFA::SetNumRowsAndTrilinosIndices( ManagedGroup * const nodeManager,
                                                                 localIndex & numLocalRows,
                                                                 localIndex & numGlobalRows,
                                                                 localIndex_array& localIndices,
                                                                 localIndex offset )
{
//  dim =
// domain.m_feElementManager.m_ElementRegions.begin()->second.m_ElementDimension;
  int dim = 1;

  int n_mpi_processes = 1;
  int this_mpi_process = 0;
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
  integer_array& is_ghost       = nodeManager->getReference<integer_array>(viewKeys.ghostRank);


  trilinos_index = -1;
  is_ghost = -2;

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

//  partition.SynchronizeFields(m_syncedFields, CommRegistry::lagrangeSolver02);

}


void SinglePhaseFlow_TPFA :: SetupSystem ( DomainPartition * const domain,
                                                   EpetraBlockSystem * const blockSystem )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ManagedGroup * const nodeManager = mesh->getNodeManager();

  localIndex dim = 1;//domain.m_feElementManager.m_ElementRegions.begin()->second.m_ElementDimension;
  localIndex n_ghost_rows  = 0;//domain.m_feNodeManager.GetNumGhosts();
  localIndex n_local_rows  = nodeManager->size()-n_ghost_rows;
  localIndex n_global_rows = 0;

  localIndex_array displacementIndices;
  SetNumRowsAndTrilinosIndices( nodeManager,
                                n_local_rows,
                                n_global_rows,
                                displacementIndices,
                                0 );



  // create epetra map


  Epetra_Map * const rowMap = blockSystem->SetRowMap( EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                                                      std::make_unique<Epetra_Map>( static_cast<int>(dim*n_global_rows),
                                                                                    static_cast<int>(dim*n_local_rows),
                                                                                    0,
                                                                                    m_linearSolverWrapper.m_epetraComm ) );

  Epetra_FECrsGraph * const sparsity = blockSystem->SetSparsity( EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                                                                 EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                                                                 std::make_unique<Epetra_FECrsGraph>(Copy,*rowMap,0) );



  SetSparsityPattern( domain, sparsity );

  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();

  blockSystem->SetMatrix( EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                          EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                          std::make_unique<Epetra_FECrsMatrix>(Copy,*sparsity) );

  blockSystem->SetSolutionVector( EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

  blockSystem->SetResidualVector( EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

}

void SinglePhaseFlow_TPFA::SetSparsityPattern( DomainPartition const * const domain,
                                                       Epetra_FECrsGraph * const sparsity )
{
  int dim=3;
  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ManagedGroup const * const nodeManager = mesh->getNodeManager();

  localIndex_array const & trilinos_index = nodeManager->getReference<localIndex_array>(viewKeys.trilinosIndex);
  ElementRegionManager const * const elemManager = mesh->getElemManager();


  elemManager->forElementRegions([&](ElementRegion const * const elementRegion)
    {
      auto const & numMethodName = elementRegion->getData<string>(keys::numericalMethod);

      elementRegion->forCellBlocks([&](CellBlockSubRegion const * const cellBlock)
      {
        localIndex const numElems = cellBlock->size();
        lArray2d const & elemsToNodes = cellBlock->getWrapper<lArray2d>(cellBlock->m_CellBlockSubRegionViewKeys.nodeList)->reference();// getData<lArray2d>(keys::nodeList);
        localIndex const numNodesPerElement = elemsToNodes.size(1);

        integer_array elementLocalDofIndex (numNodesPerElement);

        for( localIndex k=0 ; k<numElems ; ++k )
        {
          const localIndex* const localNodeIndices = elemsToNodes[k];

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
      });
    });
}



real64 SinglePhaseFlow_TPFA::Assemble ( DomainPartition * const  domain,
                                                EpetraBlockSystem * const blockSystem,
                                                real64 const time_n,
                                                real64 const dt )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ManagedGroup * const nodeManager = mesh->getNodeManager();
  ConstitutiveManager const * const constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);
//  ElementRegionManager * const elemManager = mesh->getElemManager();
//  FiniteElementManager const * numericalMethodManager = domain->getParent()->GetGroup<FiniteElementManager>(keys::finiteElementManager);
//  FiniteElementSpaceManager const * feSpaceManager = numericalMethodManager->GetGroup<FiniteElementSpaceManager>(keys::finiteElementSpaces);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  FaceManager * const faceManager = mesh->getFaceManager();

  Array2dT<localIndex> const facesToElements;
  view_rtype_const<integer_array> faceGhostRank = faceManager->getData<integer_array>(viewKeys.ghostRank);



  Array2dT<localIndex> & elemRegionList     = faceManager->elementRegionList();
  Array2dT<localIndex> & elemSubRegionList  = faceManager->elementSubRegionList();
  Array2dT<localIndex> & elemList           = faceManager->elementList();


  Epetra_FECrsMatrix * const matrix = blockSystem->GetMatrix( EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                                                              EpetraBlockSystem::BlockIDs::fluidPressureBlock );
  Epetra_FEVector * const rhs = blockSystem->GetResidualVector( EpetraBlockSystem::BlockIDs::fluidPressureBlock );
  Epetra_FEVector * const solution = blockSystem->GetSolutionVector( EpetraBlockSystem::BlockIDs::fluidPressureBlock );

  Epetra_IntSerialDenseVector elem_index(1);
  Epetra_SerialDenseVector elem_rhs(1);
  Epetra_SerialDenseMatrix elem_matrix(1, 1);

  CellBlockSubRegion * const activeCellBlock = elemManager->GetRegion( elemRegionList[0][0] )->GetGroup<CellBlockSubRegion>(elemSubRegionList[0][0]);

  auto const & constitutiveMap = activeCellBlock->getReference< std::pair< Array2dT<localIndex>,Array2dT<localIndex> > >(activeCellBlock->viewKeys().constitutiveMap);

  view_rtype_const<real64_array> volume = activeCellBlock->getData<real64_array>("Volume");
  view_rtype_const<real64_array> porosity = activeCellBlock->getData<real64_array>("Porosity");

  view_rtype_const<real64_array> dVolume = activeCellBlock->getData<real64_array>("dVolume");
  view_rtype_const<real64_array> dPorosity = activeCellBlock->getData<real64_array>("dPorosity");

  view_rtype_const<real64_array> permeability = activeCellBlock->getData<real64_array>("Permeability");



  view_rtype_const<localIndex_array> trilinos_index = activeCellBlock->getData<localIndex_array>(viewKeys.trilinosIndex);


//  constitutive::ConstitutiveBase * constitutiveModel = constitutiveManager->GetGroup<constitutive::ConstitutiveBase>( constitutiveName );
  void * constitutiveModelData;
//  constitutiveModel->SetParamStatePointers(constitutiveModelData);
//  constitutive::ConstitutiveBase::UpdateFunctionPointer
//  constitutiveUpdate = constitutiveModel->GetStateUpdateFunctionPointer();

  view_rtype_const< real64_array > pressure    = activeCellBlock->getData<real64_array>("Pressure");
  view_rtype_const< real64_array > rho        = activeCellBlock->getData<real64_array>("Density");
  view_rtype_const< real64 >       bulkModulus = activeCellBlock->getData<real64_array>("BulkModulus");

  Array2dT<localIndex> const & consitutiveModelIndex      = constitutiveMap.first;
  Array2dT<localIndex> const & consitutiveModelArrayIndex = constitutiveMap.second;

  FOR_ELEMS_IN_SUBREGION( activeCellBlock, k )
  {
//    localIndex const matIndex1 = consitutiveModelIndex[k][0];
//    localIndex const matIndex2 = consitutiveModelArrayIndex[k][0];
//
//    elem_matrix(0, 0) = 1.0 / bulkModulus[matIndex1][matIndex2] * rho[matIndex1][matIndex2] * volume[k] * porosity[k];
//
//    elem_index(0) = trilinos_index[k];
//
//    matrix->SumIntoGlobalValues(elem_index, elem_matrix);
//
//    elem_rhs(0) = ( (dPorosity[k] + porosity[k]) * (dRho[matIndex1][matIndex2] + rho[matIndex1][matIndex2]) ) * dVolume[k]
//                + (dPorosity[k] * (dRho[matIndex1][matIndex2] + rho[matIndex1][matIndex2]) + porosity[k] * dRho[matIndex1][matIndex2] ) * volume[k];
//
//    rhs->SumIntoGlobalValues(elem_index, elem_rhs);
//  });
  } END_FOR




  for( localIndex kf=0 ; kf<faceManager->size() ; ++kf )
  {
    if( faceGhostRank[kf] < 0 )
    {
      localIndex const elem1 = facesToElements[0][kf];
      localIndex const elem2 = facesToElements[1][kf];
//
//      realT const face_weight = m_faceToElemLOverA[kf][0] / (m_faceToElemLOverA[kf][0] + m_faceToElemLOverA[kf][1]);
//
//      realT const face_trans = 1.0 / (m_faceToElemLOverA[kf][0] / (*m_permeabilityPtr[reg0])[k0] + m_faceToElemLOverA[kf][1] / (*m_permeabilityPtr[reg1])[k1]);
//
//      rhoav = face_weight * (*m_densityPtr[reg0])[k0] + (1.0 - face_weight) * (*m_densityPtr[reg1])[k1];
//
//      dP = (*m_pressurePtr[reg0])[k0] - (*m_pressurePtr[reg1])[k1];
//      if (m_applyGravity) dP -= (*gdz)[kf] * rhoav;
//
//
//      rhoTrans = rhoav * face_trans / m_fluidRockProperty.m_mu * dt;
//
//      dRdP1 = dP * face_weight * (*m_densityPtr[reg0])[k0] * ( m_fluidRockProperty.m_compress + (*m_poreCompressPtr[reg0])[k0]) * face_trans / m_fluidRockProperty.m_mu * dt;
//
//      dRdP2 = dP * (1.0 - face_weight) * (*m_densityPtr[reg1])[k1] * ( m_fluidRockProperty.m_compress + (*m_poreCompressPtr[reg1])[k1]) * face_trans / m_fluidRockProperty.m_mu * dt;
//
//      dRgdP1 = 0.0;
//      dRgdP2 = 0.0;
//
//      dRgdP1 = (*gdz)[kf] * face_weight * (*m_densityPtr[reg0])[k0] * ( m_fluidRockProperty.m_compress + (*m_poreCompressPtr[reg0])[k0]);
//      dRgdP2 = (*gdz)[kf] * ( 1 - face_weight ) * (*m_densityPtr[reg1])[k1] * ( m_fluidRockProperty.m_compress + (*m_poreCompressPtr[reg1])[k1]);
//
//      face_matrix(0, 0) = rhoTrans + dRdP1 - rhoTrans * dRgdP1;
//      face_matrix(1, 1) = rhoTrans - dRdP2 + rhoTrans * dRgdP2;
//
//      face_matrix(0, 1) = -rhoTrans + dRdP2 - rhoTrans * dRgdP2;
//      face_matrix(1, 0) = -rhoTrans - dRdP1 + rhoTrans * dRgdP1;
//
//      face_index[0] = (*m_trilinos_index[reg0])[k0];
//      face_index[1] = (*m_trilinos_index[reg1])[k1];
//
//      m_matrix->SumIntoGlobalValues(face_index, face_matrix);
//
//      face_rhs(0) = rhoTrans * dP;
//      face_rhs(1) = -rhoTrans * dP;
//      m_rhs->SumIntoGlobalValues(face_index, face_rhs);

    }
  }


  matrix->GlobalAssemble(true);
  rhs->GlobalAssemble();

  return 0;
}


void SinglePhaseFlow_TPFA::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                      real64 const scalingFactor,
                                      localIndex const dofOffset,
                                      dataRepository::ManagedGroup * const nodeManager )
{
  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( EpetraBlockSystem::BlockIDs::fluidPressureBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( EpetraBlockSystem::BlockIDs::fluidPressureBlock );
  localIndex_array const & trilinos_index = nodeManager->getReference<localIndex_array>(viewKeys.trilinosIndex);

  string const & fieldName = getReference<string>(viewKeys.fieldVarName);
  real64_array & fieldVar = nodeManager->getReference<real64_array>(string("Temperature"));

  int dummy;
  double* local_solution = nullptr;
  solution->ExtractView(&local_solution,&dummy);

  for( localIndex r=0 ; r<fieldVar.size() ; ++r)
  {
//    if(is_ghost[r] < 0)
    {
      int lid = rowMap->LID(integer_conversion<int>(trilinos_index[r]));
      fieldVar[r] = local_solution[lid];
    }
  }





}


void SinglePhaseFlow_TPFA::MakeGeometryParameters( DomainPartition * const  domain )
{

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  NodeManager * const nodeManager = mesh->getNodeManager();
  FaceManager * const faceManager = mesh->getFaceManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  array<R1Tensor> faceCenter(faceManager->size());

  // matrix elements
  //const iArray1d& face_is_ghost = faceManager.GetFieldData<FieldInfo::ghostRank>();

  view_rtype_const<r1_array> nodePos = nodeManager->getData<r1_array>(nodeManager->viewKeys.referencePosition);
//  rArray1d * gdz = faceManager.GetFieldDataPointer<realT>("gdz");

//  if (m_applyGravity) (*gdz) = 0.0;

  Array2dT<localIndex> const & elemRegionList     = faceManager->elementRegionList();
  Array2dT<localIndex> const & elemSubRegionList  = faceManager->elementSubRegionList();
  Array2dT<localIndex> const & elemList           = faceManager->elementList();


  CellBlockSubRegion * const activeCellBlock = elemManager->GetRegion( elemRegionList[0][0] )->GetGroup<CellBlockSubRegion>(elemSubRegionList[0][0]);


  localIndex faceIndex, nnodes;
  R1Tensor dummy;

//  ElementRegionManager::ElementViewAccessor<real64_array> volume =
//      elemManager->ConstructViewAccessor<real64_array>( "Volume", "" );

  view_rtype<real64_array> volume     = activeCellBlock->getData<real64_array>("Volume");
  view_rtype_const<r1_array> elemCenter = activeCellBlock->getData<r1_array>("ElementCenter");
  view_rtype_const<lArray2d>   elemsToNodes = activeCellBlock->getData<lArray2d>(activeCellBlock->viewKeys().nodeList);
  view_rtype_const<r1_array>  X = nodeManager->getData<r1_array>(nodeManager->viewKeys.referencePosition);

  FOR_ELEMS_IN_SUBREGION( activeCellBlock, k )
  {
//    elemCenter[k] = 0.0;
//    for( localIndex a=0 ; a<activeCellBlock->numNodesPerElement() ; ++a )
//    {
//      elemCenter[k] += X[ elemsToNodes[a] ];
//    }
//    elemCenter[k] /= activeCellBlock->numNodesPerElement();


//    R1Tensor X7_X1  = X[elemsToNodes[7]];
//             X7_X1 -= X[elemsToNodes[1]];
//
//    R1Tensor temp;
//    volume[k] = 1.0;//1.0/12.0 * ( temp.Cross(X7_X1,));


  } END_FOR

  // loop over faces
  view_rtype<real64_array> faceArea = faceManager->getData<real64_array>("faceArea");

  for (localIndex kf = 0; kf < faceManager->size(); ++kf )
  {
    faceArea[kf] = 1.0;
  }

  view_rtype< array<R1Tensor> > faceConnectionVector = faceManager->getData<R1Tensor>("faceConnectionVector");
//  faceConnectionVector = 0.0;


  R1Tensor fCenter, fNormal;
  localIndex m, n;

  m_faceToElemLOverA.resize( faceManager->size(), 2);
  for(localIndex kf = 0; kf < faceManager->size(); ++kf)
  {

//    index = faceManager.GetParentIndex(kf);

    //if(face_is_ghost[kf] < 0) //Need ghost values for postprocessing
    {

      localIndex numElems = elemRegionList.size(0);
//      faceManager.FaceCenterAndNormal(nodeManager, kf, fCenter, fNormal);

      assert(numElems <= 2);

      for (localIndex k = 0; k < numElems; ++numElems)
      {
        m = elemRegionList(k,1);
        n = elemSubRegionList(k,1);

        R1Tensor la = elemCenter[ k ];
                 la -= fCenter;

        m_faceToElemLOverA( kf, k) = ( fabs(Dot(la, fNormal)) / faceArea[kf] );

//        if (m_applyGravity)
//        {
//          R1Tensor dz(la);
//          if (ele == 1) dz *= -1.0;
//          (*gdz)[index] += Dot(dz, m_gravityVector);
//        }
//        faceConnectionVector[index] += (0.5- ele) * 2 * la;
      }
//      faceConnectionVector[index].Normalize();
    }
  }
}



REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseFlow_TPFA, std::string const &, ManagedGroup * const )
} /* namespace ANST */
