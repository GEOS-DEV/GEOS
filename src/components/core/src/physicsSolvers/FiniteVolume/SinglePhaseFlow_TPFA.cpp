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

/**
 * @file SinglePhaseFlow_TPFA.cpp
 */

#include "SinglePhaseFlow_TPFA.hpp"

#include <vector>
#include <math.h>


#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "finiteElement/FiniteElementManager.hpp"
#include "finiteElement/FiniteElementSpaceManager.hpp"
#include "finiteElement/ElementLibrary/FiniteElement.h"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "physicsSolvers/BoundaryConditions/BoundaryConditionManager.hpp"
#include "systemSolverInterface/LinearSolverWrapper.hpp"
#include "systemSolverInterface/EpetraBlockSystem.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace systemSolverInterface;


SinglePhaseFlow_TPFA::SinglePhaseFlow_TPFA( const std::string& name,
                                            ManagedGroup * const parent ):
  SolverBase( name, parent )
{

  // set the blockID for the block system interface
  m_linearSystem->SetBlockID( EpetraBlockSystem::BlockIDs::fluidPressureBlock, this->getName() );

  // register data members with the repository
  this->RegisterViewWrapper( viewKeyStruct::gravityFlagString, &m_gravityFlag, 0 );
}


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

  docNode->AllocateChildNode( viewKeyStruct::gravityFlagString,
                              viewKeyStruct::gravityFlagString,
                              -1,
                              "integer",
                              "integer",
                              "name of field variable",
                              "name of field variable",
                              "0",
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

    {
      FaceManager * const faceManager = meshLevel->getFaceManager();
      cxx_utilities::DocumentationNode * const docNode = faceManager->getDocumentationNode();

      docNode->AllocateChildNode( viewKeyStruct::faceAreaString,
                                  viewKeyStruct::faceAreaString,
                                  -1,
                                  "real64_array",
                                  "real64_array",
                                  "",
                                  "",
                                  "",
                                  faceManager->getName(),
                                  1,
                                  0,
                                  0 );

      docNode->AllocateChildNode( viewKeyStruct::faceCenterString,
                                  viewKeyStruct::faceCenterString,
                                  -1,
                                  "r1_array",
                                  "r1_array",
                                  "",
                                  "",
                                  "",
                                  faceManager->getName(),
                                  1,
                                  0,
                                  0 );
    }



    ElementRegionManager * const elemManager = meshLevel->getElemManager();

    elemManager->forCellBlocks( [&]( CellBlockSubRegion * const cellBlock ) -> void
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

        docNode->AllocateChildNode( viewKeyStruct::deltaFluidPressureString,
                                    viewKeyStruct::deltaFluidPressureString,
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

        docNode->AllocateChildNode( viewKeyStruct::deltaFluidDensityString,
                                    viewKeyStruct::deltaFluidDensityString,
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

        docNode->AllocateChildNode( viewKeyStruct::fluidPressureString,
                                    viewKeyStruct::fluidPressureString,
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

        docNode->AllocateChildNode( viewKeyStruct::volumeString,
                                    viewKeyStruct::volumeString,
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

        docNode->AllocateChildNode( viewKeyStruct::deltaVolumeString,
                                    viewKeyStruct::deltaVolumeString,
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

        docNode->AllocateChildNode( viewKeyStruct::porosityString,
                                    viewKeyStruct::porosityString,
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
        docNode->AllocateChildNode( viewKeyStruct::deltaPorosityString,
                                    viewKeyStruct::deltaPorosityString,
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

        docNode->AllocateChildNode( viewKeyStruct::permeabilityString,
                                    viewKeyStruct::permeabilityString,
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

void SinglePhaseFlow_TPFA::InitializeFinalLeaf( ManagedGroup * const problemManager )
{
  DomainPartition * domain = problemManager->GetGroup<DomainPartition>(keys::domain);

  // Call function to fill geometry parameters for use forming system
  MakeGeometryParameters( domain );
}

void SinglePhaseFlow_TPFA::TimeStep( real64 const& time_n,
                                     real64 const& dt,
                                     const int cycleNumber,
                                     ManagedGroup * domain )
{
  // currently the only method is implcit time integration
  this->NonlinearImplicitStep( time_n, dt, cycleNumber, domain->group_cast<DomainPartition*>() );
}


/**
 * This function currently applies Dirichlet boundary conditions on the elements/zones as they
 * hold the DOF. Futher work will need to be done to apply a Dirichlet bc to the connectors (faces)
 */
void SinglePhaseFlow_TPFA::ApplyDirichletBC_implicit( ManagedGroup * object,
                                                      real64 const time,
                                                      EpetraBlockSystem & blockSystem )
{
  BoundaryConditionManager * bcManager = BoundaryConditionManager::get();
  ElementRegionManager * elemManager = object->group_cast<ElementRegionManager *>();

  ElementRegionManager::ElementViewAccessor<localIndex_array>
  trilinosIndex = elemManager->
                  ConstructViewAccessor<localIndex_array>( viewKeyStruct::trilinosIndexString );

  ElementRegionManager::ElementViewAccessor<real64_array>
  pressure_n = elemManager->
               ConstructViewAccessor<real64_array>( viewKeyStruct::fluidPressureString );

  ElementRegionManager::ElementViewAccessor<real64_array>
  dP = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidPressureString);


  // loop through cell block sub-regions
  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegion * const elemRegion = elemManager->GetRegion(er);
    for( localIndex esr=0 ; esr<elemRegion->numSubRegions() ; ++esr)
    {
      CellBlockSubRegion * const subRegion = elemRegion->GetSubRegion(esr);



      // call the BoundaryConditionManager::ApplyBoundaryCondition function that will check to see
      // if the boundary condition should be applied to this subregion
      bcManager->ApplyBoundaryCondition( time,
                                         subRegion,
                                         viewKeyStruct::fluidPressureString,
                                         [&]( BoundaryConditionBase const * const bc,
                                             lSet const & set ) -> void
        {

      // call the application of the boundray condition to alter the matrix and rhs
      bc->ApplyDirichletBounaryConditionDefaultMethod<0>( set,
                                                          time,
                                                          trilinosIndex[er][esr].get(),
                                                          1,
                                                          m_linearSystem,
                                                          EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                                                          [&](localIndex const a)->real64
        {
        return pressure_n[er][esr][a] + dP[er][esr][a];
        });
      });
    }
  }
}


void SinglePhaseFlow_TPFA::ImplicitStepSetup( real64 const& time_n,
                                                  real64 const& dt,
                                                  DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  ConstitutiveManager * const
  constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);

  auto
  dPressure = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidPressureString);

  auto
  dRho = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidDensityString);

  auto
  constitutiveMap = elemManager->
                    ConstructViewAccessor< std::pair< Array2dT<localIndex>,Array2dT<localIndex> > >( CellBlockSubRegion::viewKeyStruct::constitutiveMapString,                                                                                                          string() );

  //***** loop over all elements and initialize the dRho_dP array *****
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k)->void
  {
    dPressure[er][esr][k] = 0.0;
    dRho[er][esr][k] = 0.0;

    // initialize dRho_dP
    localIndex const matIndex1 = constitutiveMap[er][esr].get().first[k][0];
    localIndex const matIndex2 = constitutiveMap[er][esr].get().second[k][0];
    ConstitutiveBase * const EOS = constitutiveManager->GetGroup<ConstitutiveBase>(matIndex1);
    EOS->EquationOfStateDensityUpdate( dPressure[er][esr][k],
                                       matIndex2,
                                       dRho[er][esr][k],
                                       m_dRho_dP[er][esr][k] );

  });


  // setup dof numbers and linear system
  SetupSystem( domain, m_linearSystem );

  // extract the system out of the block system interface.
  Epetra_FECrsMatrix * const
  matrix = m_linearSystem->GetMatrix( EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                                     EpetraBlockSystem::BlockIDs::fluidPressureBlock );

  Epetra_FEVector * const
  rhs = m_linearSystem->GetResidualVector( EpetraBlockSystem::BlockIDs::fluidPressureBlock );

  Epetra_FEVector * const
  solution = m_linearSystem->GetSolutionVector( EpetraBlockSystem::BlockIDs::fluidPressureBlock );

  // zero out the system
  matrix->Scale(0.0);
  rhs->Scale(0.0);
  solution->Scale(0.0);
}

void SinglePhaseFlow_TPFA::ImplicitStepComplete( real64 const & time_n,
                                                     real64 const & dt,
                                                     DomainPartition * const domain)
{

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  auto
  pressure = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::fluidPressureString);

  auto
  dP = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidPressureString);

  auto
  dFluidDensity =
  elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidDensityString);


  ConstitutiveManager * const
  constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);

  constitutive::ViewAccessor< array<real64> >
  rho = constitutiveManager->GetStateData< array<real64> >("fluidDensity");

  auto constitutiveMap = elemManager->
                         ConstructViewAccessor< std::pair< Array2dT<localIndex>,Array2dT<localIndex> > >( CellBlockSubRegion::viewKeyStruct::constitutiveMapString,                                                                                                          string() );

  //***** Loop over all elements and update the pressure and density *****
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k)->void
  {
    // matIndex1 is the index of the material contained in the element
    localIndex const matIndex1 = constitutiveMap[er][esr].get().first[k][0];

    // matIndex2 is the index of the point within material specified in matIndex1
    localIndex const matIndex2 = constitutiveMap[er][esr].get().second[k][0];

    // update the fluid pressure and density.
    pressure[er][esr][k] += dP[er][esr][k];
    rho[matIndex1][matIndex2] += dFluidDensity[er][esr][k];

    if( verboseLevel() >= 1 )
    std::cout<<"pressure["<<er<<"]["<<esr<<"][k] = "<<pressure[er][esr][k]<<std::endl;
  });

}




void SinglePhaseFlow_TPFA::SetNumRowsAndTrilinosIndices( MeshLevel * const meshLevel,
                                                         localIndex & numLocalRows,
                                                         localIndex & numGlobalRows,
                                                         localIndex_array& localIndices,
                                                         localIndex offset )
{

  ElementRegionManager * const elementRegionManager = meshLevel->getElemManager();
  ElementRegionManager::ElementViewAccessor<localIndex_array>
  trilinosIndex = elementRegionManager->
                  ConstructViewAccessor<localIndex_array>( viewKeys.trilinosIndex.Key(),
                                                           string() );

  ElementRegionManager::ElementViewAccessor< integer_array >
  ghostRank = elementRegionManager->
              ConstructViewAccessor<integer_array>( ObjectManagerBase::viewKeyStruct::ghostRankString );

  int numMpiProcesses;
  MPI_Comm_size( MPI_COMM_WORLD, &numMpiProcesses );

  int thisMpiProcess = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &thisMpiProcess );

  localIndex numLocalRowsToSend = numLocalRows;
  array<localIndex> gather(numMpiProcesses);

  // communicate the number of local rows to each process
  m_linearSolverWrapper->m_epetraComm.GatherAll( &numLocalRowsToSend,
                                                &gather.front(),
                                                1 );

  GEOS_ASSERT( numLocalRows == numLocalRowsToSend, "number of local rows inconsistent" );

  // find the first local row on this partition, and find the number of total global rows.
  localIndex firstLocalRow = 0;
  numGlobalRows = 0;

  for( integer p=0 ; p<numMpiProcesses ; ++p)
  {
    numGlobalRows += gather[p];
    if(p<thisMpiProcess)
      firstLocalRow += gather[p];
  }

  // create trilinos dof indexing, setting initial values to -1 to indicate unset values.
  for( localIndex er=0 ; er<ghostRank.size() ; ++er )
  {
    for( localIndex esr=0 ; esr<ghostRank[er].size() ; ++esr )
    {
      trilinosIndex[er][esr] = -1;
    }
  }

  // loop over all elements and set the dof number if the element is not a ghost
  integer localCount = 0;
  forAllElemsInMesh( meshLevel, [&]( localIndex const er,
                                     localIndex const esr,
                                     localIndex const k)->void
  {
    if( ghostRank[er][esr][k] < 0 )
    {
      trilinosIndex[er][esr][k] = firstLocalRow+localCount+offset;
      ++localCount;
    }
    else
    {
      trilinosIndex[er][esr][k] = -1;
    }
  });

  GEOS_ASSERT(localCount == numLocalRows, "Number of DOF assigned does not match numLocalRows" );
}


void SinglePhaseFlow_TPFA :: SetupSystem ( DomainPartition * const domain,
                                           EpetraBlockSystem * const blockSystem )
{
  // assume that there is only a single MeshLevel for now
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elementRegionManager = mesh->getElemManager();

  // for this solver, the dof are on the cell center, and the row corrosponds to an element
  localIndex numGhostRows  = 0;
  localIndex numLocalRows  = 0;
  localIndex numGlobalRows = 0;

  // get the number of local elements, and ghost elements...i.e. local rows and ghost rows
  elementRegionManager->forCellBlocks( [&]( CellBlockSubRegion * const subRegion )
    {
    localIndex subRegionGhosts = subRegion->GetNumberOfGhosts();
    numGhostRows += subRegionGhosts;
    numLocalRows += subRegion->size() - subRegionGhosts;
    });


  localIndex_array displacementIndices;
  SetNumRowsAndTrilinosIndices( mesh,
                                numLocalRows,
                                numGlobalRows,
                                displacementIndices,
                                0 );

  std::map<string, array<string> > fieldNames;
  fieldNames["node"].push_back(viewKeys.trilinosIndex.Key());

  CommunicationTools::
  SynchronizeFields(fieldNames,
                    mesh,
                    domain->getReference< array<NeighborCommunicator> >( domain->viewKeys.neighbors ) );


  // construct row map, and set a pointer to the row map
  Epetra_Map * const
  rowMap = blockSystem->
           SetRowMap( EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                      std::make_unique<Epetra_Map>( static_cast<int>(m_dim*numGlobalRows),
                                                    static_cast<int>(m_dim*numLocalRows),
                                                    0,
                                                    m_linearSolverWrapper->m_epetraComm ) );

  // construct sparisty matrix, set a pointer to the sparsity pattern matrix
  Epetra_FECrsGraph * const
  sparsity = blockSystem->SetSparsity( EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                                       EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                                       std::make_unique<Epetra_FECrsGraph>(Copy,*rowMap,0) );



  // set the sparsity patter
  SetSparsityPattern( domain, sparsity );

  // assemble the global sparsity matrix
  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();

  // construct system matrix
  blockSystem->SetMatrix( EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                          EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                          std::make_unique<Epetra_FECrsMatrix>(Copy,*sparsity) );

  // construct solution vector
  blockSystem->SetSolutionVector( EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

  // construct residual vector
  blockSystem->SetResidualVector( EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

}

void SinglePhaseFlow_TPFA::SetSparsityPattern( DomainPartition const * const domain,
                                               Epetra_FECrsGraph * const sparsity )
{
  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elementRegionManager = meshLevel->getElemManager();
  ElementRegionManager::ElementViewAccessor<localIndex_array const>
  trilinosIndex = elementRegionManager->
                  ConstructViewAccessor<localIndex_array>( viewKeys.trilinosIndex.Key() );

  ElementRegionManager::ElementViewAccessor< integer_array const >
  elemGhostRank = elementRegionManager->
              ConstructViewAccessor<integer_array>( ObjectManagerBase::viewKeyStruct::ghostRankString );


  FaceManager const * const faceManager = meshLevel->getFaceManager();
  Array2dT<localIndex> const & elementRegionList = faceManager->elementRegionList();
  Array2dT<localIndex> const & elementSubRegionList = faceManager->elementSubRegionList();
  Array2dT<localIndex> const & elementIndexList = faceManager->elementList();
  integer_array const & faceGhostRank = faceManager->GhostRank();
  integer_array elementLocalDofIndex;
  elementLocalDofIndex.resize(2);

  //**** loop over all faces. Fill in sparsity for all pairs of DOF/elem that are connected by face
  for( auto const kf : m_faceConnectors )
  {
    elementLocalDofIndex[0] = trilinosIndex[elementRegionList[kf][0]][elementSubRegionList[kf][0]][elementIndexList[kf][0]];
    elementLocalDofIndex[1] = trilinosIndex[elementRegionList[kf][1]][elementSubRegionList[kf][1]][elementIndexList[kf][1]];

    sparsity->InsertGlobalIndices( 2,
                                   elementLocalDofIndex.data(),
                                   2,
                                   elementLocalDofIndex.data());
  }

  // loop over all elements and add all locals just in case the above connector loop missed some
  forAllElemsInMesh( meshLevel, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k)->void
  {
    if( elemGhostRank[er][esr][k] < 0 )
    {
      elementLocalDofIndex[0] = trilinosIndex[er][esr][k];

      sparsity->InsertGlobalIndices( 1,
                                     elementLocalDofIndex.data(),
                                     1,
                                     elementLocalDofIndex.data());
    }
  });

}



real64 SinglePhaseFlow_TPFA::Assemble ( DomainPartition * const  domain,
                                        EpetraBlockSystem * const blockSystem,
                                        real64 const time_n,
                                        real64 const dt )
{
  //***** extract data required for assembly of system *****
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ManagedGroup * const nodeManager = mesh->getNodeManager();
  ConstitutiveManager * const constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  FaceManager const * const faceManager = mesh->getFaceManager();

  integer_array const & faceGhostRank = faceManager->getReference<integer_array>(ObjectManagerBase::viewKeyStruct::ghostRankString);



  Array2dT<localIndex> const & faceToElemRegionList     = faceManager->elementRegionList();
  Array2dT<localIndex> const & faceToElemSubRegionList  = faceManager->elementSubRegionList();
  Array2dT<localIndex> const & faceToElemList           = faceManager->elementList();


  Epetra_FECrsMatrix * const dRdP = blockSystem->GetMatrix( EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                                                              EpetraBlockSystem::BlockIDs::fluidPressureBlock );
  Epetra_FEVector * const residual = blockSystem->GetResidualVector( EpetraBlockSystem::BlockIDs::fluidPressureBlock );
  Epetra_FEVector * const solution = blockSystem->GetSolutionVector( EpetraBlockSystem::BlockIDs::fluidPressureBlock );


  constitutive::ViewAccessor<real64> fluidBulkModulus = constitutiveManager->GetParameterData< real64 >("BulkModulus");

  constitutive::ViewAccessor< array<real64> >
  rho = constitutiveManager->GetStateData< array<real64> >("fluidDensity");



  auto permeability = elemManager->
                      ConstructViewAccessor<real64_array>( viewKeyStruct::permeabilityString );

  auto constitutiveMap = elemManager->
                         ConstructViewAccessor< std::pair< Array2dT<localIndex>,Array2dT<localIndex> > >( CellBlockSubRegion::viewKeyStruct::constitutiveMapString,                                                                                                          string() );

  auto pressure_n = elemManager->
                  ConstructViewAccessor<real64_array>( viewKeyStruct::fluidPressureString );


  auto trilinosIndex = elemManager->
                       ConstructViewAccessor<localIndex_array>( viewKeyStruct::trilinosIndexString );

  auto
  dRho = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidDensityString);

  auto
  elemGhostRank = elemManager->ConstructViewAccessor<integer_array>( ObjectManagerBase::
                                                                     viewKeyStruct::
                                                                     ghostRankString );


  auto volume = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::volumeString);
  auto porosity_n = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::porosityString);
  auto dP = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidPressureString);
  auto dVolume = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaVolumeString);
  auto dPorosity = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaPorosityString);

  constitutive::ViewAccessor< real64 >
  fluidViscosity = constitutiveManager->GetParameterData< real64 >("fluidViscosity");



  //***** Loop over all elements and assemble the change in volume/density terms *****
  Epetra_IntSerialDenseVector elemDOF(1);
  Epetra_SerialDenseVector localElemResidual(1);
  Epetra_SerialDenseMatrix localElem_dRdP(1, 1);

  localIndex numLocalDOF = 0;
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k)->void
  {
    if( elemGhostRank[er][esr][k]<0 )
    {
      // matIndex1 is the index of the material contained in the element
      localIndex const matIndex1 = constitutiveMap[er][esr].get().first[k][0];

      // matIndex2 is the index of the point within material specified in matIndex1
      localIndex const matIndex2 = constitutiveMap[er][esr].get().second[k][0];

      elemDOF(0) = trilinosIndex[er][esr][k];

      // under the assumption that pore volume change = volume change
      dPorosity[er][esr][k] = ( dVolume[er][esr][k] * ( 1.0 - porosity_n[er][esr][k]) )
                              /( volume[er][esr][k] + dVolume[er][esr][k] );

      // Residual contribution is mass conservation in the cell
      localElemResidual(0) = ( (dPorosity[er][esr][k] + porosity_n[er][esr][k])
                             * (dRho[er][esr][k] + rho[matIndex1][matIndex2])
                             * (volume[er][esr][k] + dVolume[er][esr][k] ) )
                           - porosity_n[er][esr][k] * rho[matIndex1][matIndex2] * volume[er][esr][k] ;

      real64 dVdP = 0.0; // is this always zero, even in a coupled problem?
      real64 dPorositydP = volume[er][esr][k] * ( 1.0 - porosity_n[er][esr][k] )
                         / ( ( volume[er][esr][k] + dVolume[er][esr][k] )
                           * ( volume[er][esr][k] + dVolume[er][esr][k] ) ) * dVdP ;
      // Derivative of residual wrt to pressure in the cell
      localElem_dRdP(0, 0) = ( dPorositydP * (dRho[er][esr][k] + rho[matIndex1][matIndex2])
                                           * (volume[er][esr][k] + dVolume[er][esr][k] ) )
                           + ( m_dRho_dP[er][esr][k] * (dPorosity[er][esr][k] + porosity_n[er][esr][k])
                                                   * (volume[er][esr][k] + dVolume[er][esr][k] ) )
                           + ( dVdP * (dPorosity[er][esr][k] + porosity_n[er][esr][k])
                                    * (dRho[er][esr][k] + rho[matIndex1][matIndex2]) );

      // add contribution to global residual and dRdP
      residual->SumIntoGlobalValues(elemDOF, localElemResidual);
      dRdP->SumIntoGlobalValues(elemDOF, localElem_dRdP);
      ++numLocalDOF;
    }
  });



  Epetra_IntSerialDenseVector faceDOF(2);
  Epetra_SerialDenseVector localFaceResidual(2);
  Epetra_SerialDenseMatrix localFace_dRdP(2, 2);

  //***** Now loop over all faces/connectors to calculate the flux contributions *****
  for( auto const kf : m_faceConnectors )
  {
      // get the maps to the element indices
      localIndex er1 = faceToElemRegionList[kf][0];
      localIndex er2 = faceToElemRegionList[kf][1];
      localIndex esr1 = faceToElemSubRegionList[kf][0];
      localIndex esr2 = faceToElemSubRegionList[kf][1];
      localIndex ei1 = faceToElemList[kf][0];
      localIndex ei2 = faceToElemList[kf][1];


      faceDOF[0] = trilinosIndex[er1][esr1][ei1];
      faceDOF[1] = trilinosIndex[er2][esr2][ei2];

      // get the constitutive indices
      localIndex const consitutiveModelIndex1      = constitutiveMap[er1][esr1].get().first[ei1][0];
      localIndex const consitutiveModelIndex2      = constitutiveMap[er2][esr2].get().first[ei2][0];
      localIndex const consitutiveModelArrayIndex1 = constitutiveMap[er1][esr1].get().second[ei1][0];
      localIndex const consitutiveModelArrayIndex2 = constitutiveMap[er2][esr2].get().second[ei2][0];

      real64 const mu = fluidViscosity[consitutiveModelIndex1];

      real64 const rho1 = rho[consitutiveModelIndex1][consitutiveModelArrayIndex1];
      real64 const rho2 = rho[consitutiveModelIndex2][consitutiveModelArrayIndex2];

      //***** calculation of flux *****
      real64 const face_weight = m_faceToElemLOverA[kf][0] / ( m_faceToElemLOverA[kf][0]
                                                             + m_faceToElemLOverA[kf][1] );

      real64 const face_trans = 1.0 / ( m_faceToElemLOverA[kf][0] / permeability[er1][esr1][ei1]
                                      + m_faceToElemLOverA[kf][1] / permeability[er2][esr2][ei2] );

      real64 const rhoBar = face_weight * rho1 + (1.0 - face_weight) * rho2;

      real64 dRhoBardP1 =         face_weight * m_dRho_dP[er1][esr1][ei1];
      real64 dRhoBardP2 = (1.0 - face_weight) * m_dRho_dP[er2][esr2][ei2];


      real64 const deltaP = pressure_n[er1][esr1][ei1] - pressure_n[er2][esr2][ei2]
                          + dP[er1][esr1][ei1] - dP[er2][esr2][ei2];

      real64 const gravityForce = m_gravityFlag ? m_gravityForce[kf] * rhoBar : 0.0;

      real64 const rhoTrans = rhoBar * face_trans / mu * dt;
      //***** end flux terms *****


      // face residual is conservation of mass flux across face
      localFaceResidual(0) =  rhoTrans * ( deltaP - gravityForce );
      localFaceResidual(1) = -localFaceResidual(0);

      // derivatives of residuals wrt to pressure
      real64 const dR1dP1 = ( face_trans / mu * dt )
                            * ( ( rhoBar * ( 1 - m_gravityForce[kf] * dRhoBardP1) )
                              + dRhoBardP1 * ( deltaP - gravityForce ) );

      real64 const dR2dP1 = -dR1dP1;

      real64 const dR1dP2 = ( face_trans / mu * dt )
                            * ( ( rhoBar * ( -1 - m_gravityForce[kf] * dRhoBardP2) )
                                + dRhoBardP2 * ( deltaP - gravityForce ) );

      real64 const dR2dP2 = -dR1dP2;

      localFace_dRdP(0, 0) = dR1dP1;
      localFace_dRdP(1, 1) = dR2dP2;
      localFace_dRdP(0, 1) = dR1dP2;
      localFace_dRdP(1, 0) = dR2dP1;

      dRdP->SumIntoGlobalValues(faceDOF, localFace_dRdP);
      residual->SumIntoGlobalValues(faceDOF, localFaceResidual);
  }

  dRdP->GlobalAssemble(true);
  residual->GlobalAssemble();




  // apply pressure boundary conditions.
  ApplyDirichletBC_implicit( elemManager,
                             time_n + dt,
                             *m_linearSystem );

  real64 localResidual = 0.0;
//  residual->Norm2(&scalarResidual);

  real64 * residualData = nullptr;
  int dummy;
  residual->ExtractView(&residualData,&dummy);
  for( localIndex i=0 ; i<numLocalDOF ; ++i )
  {
    localResidual += residualData[i]*residualData[i];
  }
  realT globalResidualNorm;
  MPI_Allreduce (&localResidual,&globalResidualNorm,1,MPI_DOUBLE,MPI_SUM ,MPI_COMM_WORLD);


  return sqrt(globalResidualNorm);
}


void SinglePhaseFlow_TPFA::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                                real64 const scalingFactor,
                                                localIndex const dofOffset,
                                                DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( EpetraBlockSystem::BlockIDs::fluidPressureBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( EpetraBlockSystem::BlockIDs::fluidPressureBlock );

  int dummy;
  double* local_solution = nullptr;
  solution->ExtractView(&local_solution,&dummy);

  ElementRegionManager * const elementRegionManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor<localIndex_array>
  trilinosIndex = elementRegionManager->
                  ConstructViewAccessor<localIndex_array>( viewKeys.trilinosIndex.Key() );

  ElementRegionManager::ElementViewAccessor<real64_array>
  dP = elementRegionManager->
       ConstructViewAccessor<real64_array>( viewKeyStruct::deltaFluidPressureString );

  ConstitutiveManager * const constitutiveManager =
      domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);

  auto
  constitutiveMap = elementRegionManager->
                    ConstructViewAccessor< std::pair< Array2dT<localIndex>,
                                                      Array2dT<localIndex> > >( CellBlockSubRegion::
                                                                                viewKeyStruct::
                                                                                constitutiveMapString );

  auto
  dRho = elementRegionManager->
         ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidDensityString);

  auto
  elemGhostRank = elemManager->ConstructViewAccessor<integer_array>( ObjectManagerBase::
                                                                     viewKeyStruct::
                                                                     ghostRankString );

  // loop over all elements to update incremental pressure
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k)->void
  {
    if( elemGhostRank[er][esr][k]<0 )
    {
      localIndex const matIndex1 = constitutiveMap[er][esr].get().first[k][0];
      localIndex const matIndex2 = constitutiveMap[er][esr].get().second[k][0];

      // extract solution and apply to dP
      int const lid = rowMap->LID(integer_conversion<int>(trilinosIndex[er][esr][k]));
      dP[er][esr][k] += local_solution[lid];
    }
  });


  // TODO Sync dP

  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k)->void
  {
    // update dRho though EOS call
    ConstitutiveBase * const EOS = constitutiveManager->GetGroup<ConstitutiveBase>(matIndex1);
    EOS->EquationOfStateDensityUpdate( dP[er][esr][k],
                                       matIndex2,
                                       dRho[er][esr][k],
                                       m_dRho_dP[er][esr][k] );
  });
}


void SinglePhaseFlow_TPFA::MakeGeometryParameters( DomainPartition * const  domain )
{

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();
  FaceManager * const faceManager = mesh->getFaceManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  R1Tensor const & gravityVector = getGravityVector();

  Array2dT<localIndex> const & elemRegionList     = faceManager->elementRegionList();
  Array2dT<localIndex> const & elemSubRegionList  = faceManager->elementSubRegionList();
  Array2dT<localIndex> const & elemList           = faceManager->elementList();
  r1_array const & X = nodeManager->referencePosition();

  auto elemCenter = elemManager->ConstructViewAccessor< r1_array >( CellBlock::
                                                                    viewKeyStruct::
                                                                    elementCenterString );
  auto elemsToNodes = elemManager->
                      ConstructViewAccessor<FixedOneToManyRelation>( CellBlockSubRegion::viewKeyStruct::nodeListString );
  auto volume       = elemManager->ConstructViewAccessor<real64_array>( viewKeyStruct::volumeString );

  integer_array const &
  faceGhostRank = faceManager->getReference<integer_array>(ObjectManagerBase::
                                                           viewKeyStruct::
                                                           ghostRankString);

  // Loop over all the elements and calculate element centers, and element volumes
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k)->void
  {
    arrayView1d<localIndex> nodeList = elemsToNodes[er][esr][k];
    array< R1Tensor > Xlocal;
    Xlocal.resize(nodeList.size());
    elemCenter[er][esr][k] = 0.0;
    for( localIndex a=0 ; a<nodeList.size() ; ++a )
    {
      Xlocal[a] = X[ nodeList[a] ];
      elemCenter[er][esr][k] += Xlocal[a];
    }
    elemCenter[er][esr][k] /= nodeList.size();

    volume[er][esr][k] = computationalGeometry::HexVolume( Xlocal );
  });


  r1_array & faceCenter = faceManager->getReference<r1_array>(viewKeyStruct::faceCenterString);
  real64_array & faceArea = faceManager->getReference<real64_array>(viewKeyStruct::faceAreaString);
  array< array<localIndex> > const & faceToNodes = faceManager->nodeList();
  array<R1Tensor> faceConnectionVector( faceManager->size() );


  R1Tensor fNormal;

  m_gravityForce.resize(faceManager->size());
  m_gravityForce = 0.0;
  m_faceToElemLOverA.resize( faceManager->size(), 2);

  // loop over faces and calculate faceArea, faceNormal and faceCenter
  localIndex numFaceConnectors = 0;
  m_faceConnectors.resize(faceManager->size());
  for(localIndex kf = 0 ; kf < faceManager->size() ; ++kf)
  {
    faceArea[kf] = computationalGeometry::Centroid_3DPolygon( faceToNodes[kf],
                                                              X,
                                                              faceCenter[kf],
                                                              fNormal );

    constexpr localIndex numElems = 2;
    for (localIndex ke = 0 ; ke < numElems ; ++ke)
    {
      if( elemRegionList[kf][ke] != -1 )
      {
        R1Tensor la = elemCenter[ elemRegionList[kf][ke] ]
                                [ elemSubRegionList[kf][ke] ]
                                [ elemList[kf][ke] ];
        la -= faceCenter[kf];

        m_faceToElemLOverA( kf, ke) = ( fabs(Dot(la, fNormal)) / faceArea[kf] );

        if (m_gravityFlag)
        {
          R1Tensor dz(la);
          if (ke == 1)
          {
            dz *= -1.0;
          }
          m_gravityForce[kf] += Dot(dz, getGravityVector());
        }
      }
    }
    if( faceGhostRank[kf]<0 && elemRegionList[kf][0] != -1 && elemRegionList[kf][1] != -1)
    {
      m_faceConnectors[numFaceConnectors] = kf;
      ++numFaceConnectors;
    }
  }
  m_faceConnectors.resize(numFaceConnectors);

  // temporary storage for m_dRho_dP on each element
  m_dRho_dP.resize(elemManager->numRegions());
  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegion const * const elemRegion = elemManager->GetRegion(er);
    m_dRho_dP[er].resize(elemRegion->numSubRegions());
    for( localIndex esr=0 ; esr<elemRegion->numSubRegions() ; ++esr )
    {
      CellBlockSubRegion const * const cellBlockSubRegion = elemRegion->GetSubRegion(esr);
      m_dRho_dP[er][esr].resize(cellBlockSubRegion->size());
    }
  }

}



REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseFlow_TPFA, std::string const &, ManagedGroup * const )
} /* namespace ANST */
