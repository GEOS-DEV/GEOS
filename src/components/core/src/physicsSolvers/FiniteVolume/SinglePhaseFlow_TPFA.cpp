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
  SolverBase(name, parent),
  m_precomputeDone(false)
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
  docNode->setShortDescription("An example single phase flow solver");


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

//  docNode->AllocateChildNode( viewKeyStruct::gravityFlagString,
//                              viewKeyStruct::gravityFlagString,
//                              -1,
//                              "integer",
//                              "integer",
//                              "Flag that enables/disables gravity",
//                              "Flag that enables/disables gravity",
//                              "0",
//                              "",
//                              0,
//                              1,
//                              0 );
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
                                  "Face surface area",
                                  "Face surface area",
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
                                  "Face centroid coordinates",
                                  "Face centroid coordinates",
                                  "",
                                  faceManager->getName(),
                                  1,
                                  0,
                                  0 );
      docNode->AllocateChildNode( viewKeyStruct::transmissibilityString,
                                  viewKeyStruct::transmissibilityString,
                                  -1,
                                  "real64_array",
                                  "real64_array",
                                  "Face TPFA transmissibility",
                                  "Face TPFA transmissibility",
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
        docNode->AllocateChildNode( viewKeys.cellLocalIndex.Key(),
                                    viewKeys.cellLocalIndex.Key(),
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
                                    "Change in fluid pressure",
                                    "Change in fluid pressure",
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
                                    "Change in fluid density",
                                    "Change in fluid density",
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
                                    "Fluid pressure",
                                    "Fluid pressure",
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
                                    "Cell volume",
                                    "Cell volume",
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
                                    "Change in cell volume",
                                    "Change in cell volume",
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
                                    "Porosity",
                                    "Porosity",
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
                                    "Change in porosity",
                                    "Change in porosity",
                                    "",
                                    elemManager->getName(),
                                    1,
                                    0,
                                    0 );

        docNode->AllocateChildNode( viewKeyStruct::permeabilityString,
                                    viewKeyStruct::permeabilityString,
                                    -1,
                                    "r1_array",
                                    "r1_array",
                                    "Permeability (principal values)",
                                    "Permeability (principal values)",
                                    "",
                                    elemManager->getName(),
                                    1,
                                    0,
                                    0 );
        docNode->AllocateChildNode( viewKeyStruct::gravityDepthString,
                                    viewKeyStruct::gravityDepthString,
                                    -1,
                                    "real64_array",
                                    "real64_array",
                                    "Precomputed (gravity dot depth)",
                                    "Precomputed (gravity dot depth)",
                                    "",
                                    elemManager->getName(),
                                    1,
                                    0,
                                    0 );

      });
  }
}

void SinglePhaseFlow_TPFA::FinalInitialization( ManagedGroup * const problemManager )
{
  DomainPartition * domain = problemManager->GetGroup<DomainPartition>(keys::domain);

  // Allocate additional storage for derivatives
  AllocateAuxStorage(domain);
}

void SinglePhaseFlow_TPFA::SolverStep( real64 const& time_n,
                                       real64 const& dt,
                                       const int cycleNumber,
                                       ManagedGroup * domain )
{
  // Call function to fill geometry parameters for use forming system
  // Can't call this in FinalInitialization() as field data has not been loaded there yet
  PrecomputeData(ManagedGroup::group_cast<DomainPartition *>(domain));

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
  cellLocalIndex = elemManager->
                   ConstructViewAccessor<localIndex_array>( viewKeyStruct::cellLocalIndexString );

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
                                                          cellLocalIndex[er][esr].get(),
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
                                       m_dDens_dPres[er][esr][k] );

  });


  // setup dof numbers and linear system
  SetupSystem( domain, m_linearSystem );
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
    std::cout<<"pressure["<<er<<"]["<<esr<<"]["<<k<<"] = "<<pressure[er][esr][k]<<std::endl;
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
  cellLocalIndex = elementRegionManager->
                   ConstructViewAccessor<localIndex_array>( viewKeys.cellLocalIndex.Key(),
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
      cellLocalIndex[er][esr] = -1;
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
      cellLocalIndex[er][esr][k] = firstLocalRow+localCount+offset;
      ++localCount;
    }
    else
    {
      cellLocalIndex[er][esr][k] = -1;
    }
  });

  GEOS_ASSERT(localCount == numLocalRows, "Number of DOF assigned does not match numLocalRows" );
}


void SinglePhaseFlow_TPFA::SetupSystem ( DomainPartition * const domain,
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
  fieldNames["node"].push_back(viewKeys.cellLocalIndex.Key());

  CommunicationTools::
  SynchronizeFields(fieldNames,
                    mesh,
                    domain->getReference< array<NeighborCommunicator> >( domain->viewKeys.neighbors ) );


  // construct row map, and set a pointer to the row map
  Epetra_Map * const
  rowMap = blockSystem->
           SetRowMap( EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                      std::make_unique<Epetra_Map>( static_cast<int>(numGlobalRows),
                                                    static_cast<int>(numLocalRows),
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
  cellLocalIndex = elementRegionManager->
                   ConstructViewAccessor<localIndex_array>( viewKeys.cellLocalIndex.Key() );

  ElementRegionManager::ElementViewAccessor< integer_array const >
  elemGhostRank = elementRegionManager->
              ConstructViewAccessor<integer_array>( ObjectManagerBase::viewKeyStruct::ghostRankString );


  FaceManager const * const faceManager = meshLevel->getFaceManager();
  Array2dT<localIndex> const & elementRegionList = faceManager->elementRegionList();
  Array2dT<localIndex> const & elementSubRegionList = faceManager->elementSubRegionList();
  Array2dT<localIndex> const & elementIndexList = faceManager->elementList();
  integer_array const & faceGhostRank = faceManager->GhostRank();

  integer_array elementLocalDofIndexRow;
  integer_array elementLocalDofIndexCol;

  //**** loop over all faces. Fill in sparsity for all pairs of DOF/elem that are connected by face
  for (auto const & conn : m_faceConnectors)
  {
    constexpr localIndex numElems = 2;
    elementLocalDofIndexRow.resize(numElems);
    for (localIndex ke = 0; ke < numElems; ++ke)
    {
      // DOF index is equal to cell index since only one DOF per cell
      elementLocalDofIndexRow[ke] = cellLocalIndex[conn.connectedCellIndices[ke].region   ]
                                                  [conn.connectedCellIndices[ke].subRegion]
                                                  [conn.connectedCellIndices[ke].index    ];
    }

    const localIndex stencilSize = conn.stencilCellIndices.size();
    elementLocalDofIndexCol.resize(stencilSize);
    for (localIndex ke = 0; ke < stencilSize; ++ke)
    {
      // DOF index is equal to cell index since only one DOF per cell
      elementLocalDofIndexCol[ke] = cellLocalIndex[conn.stencilCellIndices[ke].region   ]
                                                  [conn.stencilCellIndices[ke].subRegion]
                                                  [conn.stencilCellIndices[ke].index    ];

    }

    sparsity->InsertGlobalIndices(numElems,
                                  elementLocalDofIndexRow.data(),
                                  stencilSize,
                                  elementLocalDofIndexCol.data());
  }

  // loop over all elements and add all locals just in case the above connector loop missed some
  forAllElemsInMesh(meshLevel, [&] (localIndex const er,
                                    localIndex const esr,
                                    localIndex const k) -> void
  {
    if (elemGhostRank[er][esr][k] < 0)
    {
      elementLocalDofIndexRow[0] = cellLocalIndex[er][esr][k];

      sparsity->InsertGlobalIndices( 1,
                                     elementLocalDofIndexRow.data(),
                                     1,
                                     elementLocalDofIndexRow.data());
    }
  });

}



real64 SinglePhaseFlow_TPFA::AssembleSystem ( DomainPartition * const  domain,
                                        EpetraBlockSystem * const blockSystem,
                                        real64 const time_n,
                                        real64 const dt )
{
  //***** extract data required for assembly of system *****
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  ConstitutiveManager * const  constitutiveManager = domain->GetGroup<ConstitutiveManager>(keys::ConstitutiveManager);
  ElementRegionManager * const elemManager         = mesh->getElemManager();
  FaceManager const * const    faceManager         = mesh->getFaceManager();

  integer_array const & faceGhostRank = faceManager->getReference<integer_array>(ObjectManagerBase::viewKeyStruct::ghostRankString);

  Array2dT<localIndex> const & faceToElemRegionList     = faceManager->elementRegionList();
  Array2dT<localIndex> const & faceToElemSubRegionList  = faceManager->elementSubRegionList();
  Array2dT<localIndex> const & faceToElemList           = faceManager->elementList();

  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix(EpetraBlockSystem::BlockIDs::fluidPressureBlock,
                                                               EpetraBlockSystem::BlockIDs::fluidPressureBlock);
  Epetra_FEVector * const residual = blockSystem->GetResidualVector(EpetraBlockSystem::BlockIDs::fluidPressureBlock);

  jacobian->Scale(0.0);
  residual->Scale(0.0);

  constitutive::ViewAccessor<real64> fluidBulkModulus = constitutiveManager->GetParameterData<real64>("BulkModulus");
  constitutive::ViewAccessor<real64> fluidViscosity   = constitutiveManager->GetParameterData<real64>("fluidViscosity");
  constitutive::ViewAccessor<real64_array> density    = constitutiveManager->GetStateData<real64_array>("fluidDensity");

  auto constitutiveMap = elemManager->
                         ConstructViewAccessor<std::pair<Array2dT<localIndex>, Array2dT<localIndex>>>
                           (CellBlockSubRegion::viewKeyStruct::constitutiveMapString, string());

  auto
  elemGhostRank = elemManager->ConstructViewAccessor<integer_array>( ObjectManagerBase::
                                                                     viewKeyStruct::
                                                                     ghostRankString );

  auto cellLocalIndex = elemManager->ConstructViewAccessor<localIndex_array>(viewKeyStruct::cellLocalIndexString);

  auto pressure_n = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::fluidPressureString);
  auto dP         = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidPressureString);
  auto porosity_n = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::porosityString);
  auto dPorosity  = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaPorosityString);
  auto volume_n   = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::volumeString);
  auto dVolume    = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaVolumeString);
  auto dRho       = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidDensityString);
  auto gravDepth  = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::gravityDepthString);

  const real64_array & trans = faceManager->getReference<real64_array>(viewKeyStruct::transmissibilityString);

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

      elemDOF(0) = cellLocalIndex[er][esr][k];

      // under the assumption that pore volume change = volume change
      dPorosity[er][esr][k] = ( dVolume[er][esr][k] * ( 1.0 - porosity_n[er][esr][k]) )
                              /( volume_n[er][esr][k] + dVolume[er][esr][k] );

      // Residual contribution is mass conservation in the cell
      localElemResidual(0) = ( (dPorosity[er][esr][k] + porosity_n[er][esr][k])
                             * (dRho[er][esr][k] + density[matIndex1][matIndex2])
                             * (volume_n[er][esr][k] + dVolume[er][esr][k] ) )
                           - porosity_n[er][esr][k] * density[matIndex1][matIndex2] * volume_n[er][esr][k] ;

      real64 dVdP = 0.0; // is this always zero, even in a coupled problem?
      real64 dPorositydP = volume_n[er][esr][k] * ( 1.0 - porosity_n[er][esr][k] )
                         / ( ( volume_n[er][esr][k] + dVolume[er][esr][k] )
                           * ( volume_n[er][esr][k] + dVolume[er][esr][k] ) ) * dVdP ;

      // Derivative of residual wrt to pressure in the cell
      localElem_dRdP(0, 0) = ( dPorositydP * (dRho[er][esr][k] + density[matIndex1][matIndex2])
                                           * (volume_n[er][esr][k] + dVolume[er][esr][k] ) )
                           + ( m_dDens_dPres[er][esr][k] * (dPorosity[er][esr][k] + porosity_n[er][esr][k])
                                                   * (volume_n[er][esr][k] + dVolume[er][esr][k] ) )
                           + ( dVdP * (dPorosity[er][esr][k] + porosity_n[er][esr][k])
                                    * (dRho[er][esr][k] + density[matIndex1][matIndex2]) );

      // add contribution to global residual and dRdP
      residual->SumIntoGlobalValues(elemDOF, localElemResidual);
      jacobian->SumIntoGlobalValues(elemDOF, localElem_dRdP);
      ++numLocalDOF;
    }
  });


  constexpr localIndex numElems = 2;
  Epetra_IntSerialDenseVector eqnRowIndices(numElems);
  Epetra_IntSerialDenseVector dofColIndices;
  Epetra_SerialDenseVector localFlux(numElems);
  Epetra_SerialDenseMatrix localFluxJacobian;

  // temporary working arrays
  real64 densWeight[numElems] = { 0.5, 0.5 };
  real64 dens[numElems], dDens_dP[numElems];
  real64 visc[numElems], dVisc_dP[numElems];
  real64 mobi[numElems], dMobi_dP[numElems];
  real64_array dDensMean_dP, dFlux_dP;

  //***** Now loop over all faces/connectors to calculate the flux contributions *****
  for (auto const & conn : m_faceConnectors)
  {
    const localIndex stencilSize = conn.stencilCellIndices.size();

    // resize and clear local working arrays
    dDensMean_dP.resize(stencilSize); // doesn't need to be that large, but it's convenient
    dFlux_dP.resize(stencilSize);

    // clear working arrays
    dDensMean_dP = 0.0;

    // resize local matrices and vectors
    dofColIndices.Resize(integer_conversion<int>(stencilSize));
    localFluxJacobian.Shape(numElems, integer_conversion<int>(stencilSize));

    // get the column indices of participating DOFs
    for (localIndex ke = 0; ke < stencilSize; ++ke)
    {
      dofColIndices[ke] = cellLocalIndex[conn.stencilCellIndices[ke].region   ]
                                        [conn.stencilCellIndices[ke].subRegion]
                                        [conn.stencilCellIndices[ke].index    ];
    }

    // calculate quantities on primary connected cells
    real64 densMean = 0.0;
    for (localIndex ke = 0; ke < numElems; ++ke)
    {
      // get the maps to the element indices
      localIndex const er  = conn.connectedCellIndices[ke].region;
      localIndex const esr = conn.connectedCellIndices[ke].subRegion;
      localIndex const ei  = conn.connectedCellIndices[ke].index;

      eqnRowIndices[ke] = cellLocalIndex[er][esr][ei];

      localIndex const constModelIndex      = constitutiveMap[er][esr].get().first[ei][0];
      localIndex const constModelArrayIndex = constitutiveMap[er][esr].get().second[ei][0];

      // density
      dens[ke]     = density[constModelIndex][constModelArrayIndex];
      dDens_dP[ke] = m_dDens_dPres[er][esr][ei];

      // viscosity (TODO non-constant viscosity)
      visc[ke]     = fluidViscosity[constModelIndex];
      dVisc_dP[ke] = 0.0;

      // mobility
      mobi[ke]      = dens[ke] / visc[ke];
      dMobi_dP[ke]  = dDens_dP[ke] / visc[ke] - mobi[ke] / visc[ke] * dVisc_dP[ke];

      // average density
      densMean += densWeight[ke] * dens[ke];
      dDensMean_dP[ke] = densWeight[ke] * dDens_dP[ke];
    }

    //***** calculation of flux *****

    // compute potential difference MPFA-style
    real64 potDif = 0.0;
    for (localIndex ke = 0; ke < stencilSize; ++ke)
    {
      localIndex er  = conn.stencilCellIndices[ke].region;
      localIndex esr = conn.stencilCellIndices[ke].subRegion;
      localIndex ei  = conn.stencilCellIndices[ke].index;

      real64 const pres     = pressure_n[er][esr][ei] + dP[er][esr][ei];
      real64 const gravD    = gravDepth[er][esr][ei];
      real64 const gravTerm = m_gravityFlag ? densMean * gravD : 0.0;
      real64 const dGrav_dP = m_gravityFlag ? dDensMean_dP[ke] * gravD : 0.0;

      potDif += conn.stencilWeights[ke] * (pres + gravTerm);
      dFlux_dP[ke] = conn.stencilWeights[ke] * (1.0 + dGrav_dP);
    }

    // upwinding of fluid properties (make this an option?)
    localIndex const k_up = (potDif >= 0) ? 0 : 1;

    // compute the final flux and derivatives
    const real64 flux = mobi[k_up] * potDif;
    for (localIndex ke = 0; ke < stencilSize; ++ke)
      dFlux_dP[ke] *= mobi[k_up];
    dFlux_dP[k_up] += dMobi_dP[k_up] * potDif;

    //***** end flux terms *****

    // populate local flux vector and derivatives
    localFlux(0) =  flux;
    localFlux(1) = -flux;

    for (localIndex ke = 0; ke < stencilSize; ++ke)
    {
      localFluxJacobian(0, ke) =  dFlux_dP[ke];
      localFluxJacobian(1, ke) = -dFlux_dP[ke];
    }

    // Multiply by dt
    localFlux.Scale(dt);
    localFluxJacobian.Scale(dt);

    // Add to global residual/jacobian
    jacobian->SumIntoGlobalValues(eqnRowIndices, localFluxJacobian);
    residual->SumIntoGlobalValues(eqnRowIndices, localFlux);
  }

  jacobian->GlobalAssemble(true);
  residual->GlobalAssemble();




  // apply pressure boundary conditions.
  ApplyDirichletBC_implicit( elemManager,
                             time_n + dt,
                             *m_linearSystem );


  if( verboseLevel() >= 2 )
  {
    jacobian->Print(std::cout);
    residual->Print(std::cout);
  }




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
  cellLocalIndex = elementRegionManager->
                   ConstructViewAccessor<localIndex_array>( viewKeys.cellLocalIndex.Key() );

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
  elemGhostRank = elementRegionManager->ConstructViewAccessor<integer_array>( ObjectManagerBase::
                                                                              viewKeyStruct::
                                                                              ghostRankString );

  // loop over all elements to update incremental pressure
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k)->void
  {
    if( elemGhostRank[er][esr][k]<0 )
    {
      // extract solution and apply to dP
      int const lid = rowMap->LID(integer_conversion<int>(cellLocalIndex[er][esr][k]));
      dP[er][esr][k] += scalingFactor * local_solution[lid];
    }
  });


  // TODO Sync dP once element field syncing is reimplmented.
  //std::map<string, array<string> > fieldNames;
  //fieldNames["element"].push_back(viewKeyStruct::deltaFluidPressureString);
  //CommunicationTools::SynchronizeFields(fieldNames,
  //                            mesh,
  //                            domain->getReference< array<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k)->void
  {
    localIndex const matIndex1 = constitutiveMap[er][esr].get().first[k][0];
    localIndex const matIndex2 = constitutiveMap[er][esr].get().second[k][0];
    // update dRho though EOS call
    ConstitutiveBase * const EOS = constitutiveManager->GetGroup<ConstitutiveBase>(matIndex1);
    EOS->EquationOfStateDensityUpdate( dP[er][esr][k],
                                       matIndex2,
                                       dRho[er][esr][k],
                                       m_dDens_dPres[er][esr][k] );
  });
}


void SinglePhaseFlow_TPFA::PrecomputeData(DomainPartition *const domain)
{
  if (m_precomputeDone)
    return;

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

  auto volume       = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::volumeString);
  auto gravityDepth = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::gravityDepthString);
  auto permeability = elemManager->ConstructViewAccessor<r1_array>(viewKeyStruct::permeabilityString);

  integer_array const &
  faceGhostRank = faceManager->getReference<integer_array>(ObjectManagerBase::
                                                           viewKeyStruct::
                                                           ghostRankString);

  // Loop over all the elements and calculate element centers, and element volumes
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k )->void
  {
    arrayView1d<localIndex> nodeList = elemsToNodes[er][esr][k];
    array< R1Tensor > Xlocal;
    Xlocal.resize(nodeList.size());

    R1Tensor & center = elemCenter[er][esr][k];
    center = 0.0;

    for (localIndex a = 0; a < nodeList.size(); ++a)
    {
      Xlocal[a] = X[nodeList[a]];
      center += Xlocal[a];
    }
    center /= nodeList.size();

    volume[er][esr][k] = computationalGeometry::HexVolume(Xlocal);

    gravityDepth[er][esr][k] = Dot(center, getGravityVector());
  });


  r1_array & faceCenter    = faceManager->getReference<r1_array>(viewKeyStruct::faceCenterString);
  real64_array & faceArea  = faceManager->getReference<real64_array>(viewKeyStruct::faceAreaString);
  real64_array & faceTrans = faceManager->getReference<real64_array>(viewKeyStruct::transmissibilityString);
  array<array<localIndex>> const & faceToNodes = faceManager->nodeList();


  R1Tensor faceNormal, faceConormal, cellToFaceVec;
  R2Tensor permTensor, permTemp;
  faceTrans = 0.0;

  // loop over faces and calculate faceArea, faceNormal and faceCenter
  localIndex numFaceConnectors = 0;
  m_faceConnectors.resize(faceManager->size());
  for (localIndex kf = 0; kf < faceManager->size(); ++kf)
  {
    faceArea[kf] = computationalGeometry::Centroid_3DPolygon(faceToNodes[kf],
                                                             X,
                                                             faceCenter[kf],
                                                             faceNormal);

    constexpr localIndex numElems = 2;
    real64 transInv = 0.0;

    for (localIndex ke = 0; ke < numElems; ++ke)
    {
      if (elemRegionList[kf][ke] != -1)
      {
        const localIndex er  = elemRegionList[kf][ke];
        const localIndex esr = elemSubRegionList[kf][ke];
        const localIndex ei  = elemList[kf][ke];

        cellToFaceVec = faceCenter[kf];
        cellToFaceVec -= elemCenter[er][esr][ei];

        if (ke == 1)
          cellToFaceVec *= -1.0;

        const real64 c2fDistance = cellToFaceVec.Normalize();

        // assemble full permeability tensor from principal axis/components
        R1Tensor & permValues = permeability[er][esr][ei];
        permTensor = 0.0;

        // XXX: unroll loop manually?
        for (unsigned icoord = 0; icoord < 3; ++icoord)
        {
          // in future principal axis may be specified in input
          R1Tensor axis(0.0); axis(icoord) = 1.0;

          // XXX: is there a more elegant way to do this?
          permTemp.dyadic_aa(axis);
          permTemp *= permValues(icoord);
          permTensor += permTemp;
        }

        faceConormal.AijBj(permTensor, faceNormal);
        const real64 ht = Dot(cellToFaceVec, faceConormal) * faceArea[kf] / c2fDistance;

        transInv += 1.0 / ht; // XXX: safeguard against div by zero?
      }
    }

    faceTrans[kf] = 1.0 / transInv; // XXX: safeguard against div by zero?

    // ensure consistent normal orientation
    if (Dot(cellToFaceVec, faceNormal) < 0)
      faceTrans[kf] *= -1;

    if (faceGhostRank[kf] < 0 && elemRegionList[kf][0] != -1 && elemRegionList[kf][1] != -1)
    {
      CellConnection & conn = m_faceConnectors[numFaceConnectors];
      conn.resize(numElems);
      for (localIndex ke = 0; ke < numElems; ++ke)
      {
        conn.faceIndex = kf;
        conn.connectedCellIndices[ke] = { elemRegionList[kf][ke], elemSubRegionList[kf][ke], elemList[kf][ke] };

        // TPFA-specific things
        conn.stencilCellIndices[ke] = conn.connectedCellIndices[ke];
        conn.stencilWeights[ke] = faceTrans[kf] * (ke == 0 ? 1 : -1);
      }
      ++numFaceConnectors;
    }
  }
  m_faceConnectors.resize(numFaceConnectors);

  m_precomputeDone = true;
}

void SinglePhaseFlow_TPFA::AllocateAuxStorage(DomainPartition *const domain)
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  // temporary storage for m_dRho_dP on each element
  m_dDens_dPres.resize(elemManager->numRegions());
  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegion const * const elemRegion = elemManager->GetRegion(er);
    m_dDens_dPres[er].resize(elemRegion->numSubRegions());
    for( localIndex esr=0 ; esr<elemRegion->numSubRegions() ; ++esr )
    {
      CellBlockSubRegion const * const cellBlockSubRegion = elemRegion->GetSubRegion(esr);
      m_dDens_dPres[er][esr].resize(cellBlockSubRegion->size());
    }
  }
}

void SinglePhaseFlow_TPFA::SolveSystem( EpetraBlockSystem * const blockSystem,
                                        SystemSolverParameters const * const params )
{
  Epetra_FEVector * const
  solution = m_linearSystem->GetSolutionVector( EpetraBlockSystem::BlockIDs::fluidPressureBlock );

  Epetra_FEVector * const
  residual = blockSystem->GetResidualVector( EpetraBlockSystem::BlockIDs::fluidPressureBlock );
  residual->Scale(-1.0);

  solution->Scale(0.0);

  m_linearSolverWrapper->SolveSingleBlockSystem( blockSystem,
                                                 params,
                                                 EpetraBlockSystem::BlockIDs::fluidPressureBlock );

  if( verboseLevel() >= 2 )
  {
    solution->Print(std::cout);
  }

}

void SinglePhaseFlow_TPFA::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elementRegionManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor<real64_array>
  dP = elementRegionManager->
       ConstructViewAccessor<real64_array>( viewKeyStruct::deltaFluidPressureString );

  ElementRegionManager::ElementViewAccessor<real64_array>
  dRho = elementRegionManager->
         ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidDensityString);


  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k)->void
  {
    dP[er][esr][k] = 0.0;
    dRho[er][esr][k] = 0.0;
  });

}


REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseFlow_TPFA, std::string const &, ManagedGroup * const )
} /* namespace ANST */
