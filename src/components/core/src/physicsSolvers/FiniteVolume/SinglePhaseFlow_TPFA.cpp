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
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
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
  m_precomputeDone(false),
  m_faceConnectors(),
  m_gravityFlag(),
  m_dDens_dPres()
{
  // set the blockID for the block system interface
  getLinearSystemRepository()->
  SetBlockID( BlockIDs::fluidPressureBlock, this->getName() );
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

  docNode->AllocateChildNode( viewKeyStruct::gravityFlagString,
                              viewKeyStruct::gravityFlagString,
                              -1,
                              "integer",
                              "integer",
                              "Flag that enables/disables gravity",
                              "Flag that enables/disables gravity",
                              "1",
                              "",
                              1,
                              1,
                              0 );
}

void SinglePhaseFlow_TPFA::FillOtherDocumentationNodes( dataRepository::ManagedGroup * const rootGroup )
{
  SolverBase::FillOtherDocumentationNodes( rootGroup );
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

        docNode->AllocateChildNode( viewKeyStruct::deltaFluidPressureString,
                                    viewKeyStruct::deltaFluidPressureString,
                                    -1,
                                    "real64_array",
                                    "real64_array",
                                    "Change in fluid pressure",
                                    "Change in fluid pressure",
                                    "",
                                    elemManager->getName(),
                                    0,
                                    0,
                                    0 );

        docNode->AllocateChildNode( viewKeyStruct::fluidDensityString,
                                    viewKeyStruct::fluidDensityString,
                                    -1,
                                    "real64_array",
                                    "real64_array",
                                    "Fluid density",
                                    "Fluid density",
                                    "",
                                    elemManager->getName(),
                                    0,
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
                                    0,
                                    0,
                                    0 );

        docNode->AllocateChildNode( viewKeyStruct::fluidViscosityString,
                                    viewKeyStruct::fluidViscosityString,
                                    -1,
                                    "real64_array",
                                    "real64_array",
                                    "Fluid viscosity",
                                    "Fluid viscosity",
                                    "",
                                    elemManager->getName(),
                                    0,
                                    0,
                                    0 );

        docNode->AllocateChildNode( viewKeyStruct::deltaFluidViscosityString,
                                    viewKeyStruct::deltaFluidViscosityString,
                                    -1,
                                    "real64_array",
                                    "real64_array",
                                    "Change in fluid viscosity",
                                    "Change in fluid viscosity",
                                    "",
                                    elemManager->getName(),
                                    0,
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
                                    0,
                                    0,
                                    0 );

        docNode->AllocateChildNode( viewKeyStruct::referencePorosityString,
                                    viewKeyStruct::referencePorosityString,
                                    -1,
                                    "real64_array",
                                    "real64_array",
                                    "Reference porosity",
                                    "Reference porosity",
                                    "",
                                    elemManager->getName(),
                                    0,
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

        docNode->AllocateChildNode( viewKeyStruct::gravityDepthString,
                                    viewKeyStruct::gravityDepthString,
                                    -1,
                                    "real64_array",
                                    "real64_array",
                                    "Precomputed (gravity dot depth)",
                                    "Precomputed (gravity dot depth)",
                                    "",
                                    elemManager->getName(),
                                    0,
                                    0,
                                    0 );

        docNode->AllocateChildNode( viewKeyStruct::blockLocalDofNumberString,
                                    viewKeyStruct::blockLocalDofNumberString,
                                    -1,
                                    "localIndex_array",
                                    "localIndex_array",
                                    "DOF index",
                                    "DOF index",
                                    "0",
                                    "",
                                    0,
                                    0,
                                    0 );

      });
  }
}

void SinglePhaseFlow_TPFA::FinalInitialization( ManagedGroup * const problemManager )
{
  DomainPartition * domain = problemManager->GetGroup<DomainPartition>(keys::domain);

  // obtain gravity flag from the input
  m_gravityFlag = static_cast<bool>(this->getReference<integer>(viewKeyStruct::gravityFlagString));

  // Allocate additional storage for derivatives
  AllocateAuxStorage(domain);
}

real64 SinglePhaseFlow_TPFA::SolverStep( real64 const& time_n,
                                         real64 const& dt,
                                         const int cycleNumber,
                                         ManagedGroup * domain )
{
  // Call function to fill geometry parameters for use forming system
  // Can't call this in FinalInitialization() as field data has not been loaded there yet
  PrecomputeData(ManagedGroup::group_cast<DomainPartition *>(domain));

  // currently the only method is implcit time integration
  return this->NonlinearImplicitStep( time_n,
                                      dt,
                                      cycleNumber,
                                      domain->group_cast<DomainPartition*>(),
                                      getLinearSystemRepository() );
}


/**
 * This function currently applies Dirichlet boundary conditions on the elements/zones as they
 * hold the DOF. Futher work will need to be done to apply a Dirichlet bc to the connectors (faces)
 */
void SinglePhaseFlow_TPFA::ApplyDirichletBC_implicit( ManagedGroup * object,
                                                      real64 const time,
                                                      EpetraBlockSystem * const blockSystem )
{
  BoundaryConditionManager * bcManager = BoundaryConditionManager::get();
  ElementRegionManager * elemManager = object->group_cast<ElementRegionManager *>();

  ElementRegionManager::ElementViewAccessor<localIndex_array>
  blockLocalDofNumber = elemManager->
                  ConstructViewAccessor<localIndex_array>( viewKeyStruct::blockLocalDofNumberString );

  ElementRegionManager::ElementViewAccessor<real64_array>
  pres = elemManager->ConstructViewAccessor<real64_array>( viewKeyStruct::fluidPressureString );

  ElementRegionManager::ElementViewAccessor<real64_array>
  dPres = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidPressureString);


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
                                                          blockLocalDofNumber[er][esr].get(),
                                                          1,
                                                          blockSystem,
                                                          BlockIDs::fluidPressureBlock,
                                                          [&](localIndex const a)->real64
        {
        return pres[er][esr][a] + dPres[er][esr][a];
        });
      });
    }
  }
}


void SinglePhaseFlow_TPFA::
ImplicitStepSetup( real64 const& time_n,
                   real64 const& dt,
                   DomainPartition * const domain,
                   systemSolverInterface::EpetraBlockSystem * const blockSystem)
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  ConstitutiveManager * const
  constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);

  auto pres = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::fluidPressureString);
  auto dPres = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidPressureString);

  auto dens = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::fluidDensityString);
  auto dDens = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidDensityString);

  auto visc = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::fluidViscosityString);
  auto dVisc = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidViscosityString);

  auto poro = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::porosityString);
  auto dPoro = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaPorosityString);
  auto refPoro = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::referencePorosityString);

  auto
  constitutiveMap = elemManager->
                    ConstructViewAccessor< std::pair< Array2dT<localIndex>,Array2dT<localIndex> > >( CellBlockSubRegion::viewKeyStruct::constitutiveMapString,                                                                                                          string() );

  //***** loop over all elements and initialize the derivative arrays *****
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k)->void
  {
    dPres[er][esr][k] = 0.0;
    dDens[er][esr][k] = 0.0;
    dVisc[er][esr][k] = 0.0;
    dPoro[er][esr][k] = 0.0;

    // initialize dDens_dPres
    localIndex const matIndex1 = constitutiveMap[er][esr].get().first[k][0];
    localIndex const matIndex2 = constitutiveMap[er][esr].get().second[k][0];
    ConstitutiveBase * const EOS = constitutiveManager->GetGroup<ConstitutiveBase>(matIndex1);

    real64 dummy;
    real64 const pressure = pres[er][esr][k] + dPres[er][esr][k];

    EOS->FluidDensityUpdate(pressure, matIndex2, dens[er][esr][k], m_dDens_dPres[er][esr][k]);
    EOS->FluidViscosityUpdate(pressure, matIndex2, visc[er][esr][k], m_dVisc_dPres[er][esr][k]);
    EOS->SimplePorosityUpdate(pressure, refPoro[er][esr][k], matIndex2, poro[er][esr][k], m_dPoro_dPres[er][esr][k]);
  });


  // setup dof numbers and linear system
  SetupSystem( domain, blockSystem );
}

void SinglePhaseFlow_TPFA::ImplicitStepComplete( real64 const & time_n,
                                                     real64 const & dt,
                                                     DomainPartition * const domain)
{

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  auto pres = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::fluidPressureString);
  auto dPres = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidPressureString);

  auto dens = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::fluidDensityString);
  auto dDens = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidDensityString);

  auto visc = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::fluidViscosityString);
  auto dVisc = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidViscosityString);

  auto poro = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::porosityString);
  auto dPoro = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaPorosityString);

  //***** Loop over all elements and update the pressure and density *****
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k)->void
  {
    // update the fluid pressure and density.
    pres[er][esr][k] += dPres[er][esr][k];
    dens[er][esr][k] += dDens[er][esr][k];
    visc[er][esr][k] += dVisc[er][esr][k];
    poro[er][esr][k] += dPoro[er][esr][k];

    if (verboseLevel() >= 1)
    {
      std::cout << "pressure[" << er << "][" << esr << "][" << k << "] = " << pres[er][esr][k] << std::endl;
    }
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
  blockLocalDofNumber = elementRegionManager->
                        ConstructViewAccessor<localIndex_array>( viewKeys.blockLocalDofNumber.Key(), string() );

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
  m_linearSolverWrapper.m_epetraComm.GatherAll( &numLocalRowsToSend,
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
      blockLocalDofNumber[er][esr] = -1;
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
      blockLocalDofNumber[er][esr][k] = firstLocalRow+localCount+offset;
      ++localCount;
    }
    else
    {
      blockLocalDofNumber[er][esr][k] = -1;
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

  //TODO element sync doesn't work yet
//  std::map<string, array<string> > fieldNames;
//  fieldNames["element"].push_back(viewKeys.blockLocalDofNumber.Key());
//
//  CommunicationTools::
//  SynchronizeFields(fieldNames,
//                    mesh,
//                    domain->getReference< array<NeighborCommunicator> >( domain->viewKeys.neighbors ) );
//

  // construct row map, and set a pointer to the row map
  Epetra_Map * const
  rowMap = blockSystem->
           SetRowMap( BlockIDs::fluidPressureBlock,
                      std::make_unique<Epetra_Map>( static_cast<int>(numGlobalRows),
                                                    static_cast<int>(numLocalRows),
                                                    0,
                                                    m_linearSolverWrapper.m_epetraComm ) );

  // construct sparisty matrix, set a pointer to the sparsity pattern matrix
  Epetra_FECrsGraph * const
  sparsity = blockSystem->SetSparsity( BlockIDs::fluidPressureBlock,
                                       BlockIDs::fluidPressureBlock,
                                       std::make_unique<Epetra_FECrsGraph>(Copy,*rowMap,0) );



  // set the sparsity patter
  SetSparsityPattern( domain, sparsity );

  // assemble the global sparsity matrix
  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();

  // construct system matrix
  blockSystem->SetMatrix( BlockIDs::fluidPressureBlock,
                          BlockIDs::fluidPressureBlock,
                          std::make_unique<Epetra_FECrsMatrix>(Copy,*sparsity) );

  // construct solution vector
  blockSystem->SetSolutionVector( BlockIDs::fluidPressureBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

  // construct residual vector
  blockSystem->SetResidualVector( BlockIDs::fluidPressureBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

}

void SinglePhaseFlow_TPFA::SetSparsityPattern( DomainPartition const * const domain,
                                               Epetra_FECrsGraph * const sparsity )
{
  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elementRegionManager = meshLevel->getElemManager();
  ElementRegionManager::ElementViewAccessor<localIndex_array const>
  blockLocalDofNumber = elementRegionManager->
                  ConstructViewAccessor<localIndex_array>( viewKeys.blockLocalDofNumber.Key() );

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
      elementLocalDofIndexRow[ke] = blockLocalDofNumber[conn.connectedCellIndices[ke].region   ]
                                                       [conn.connectedCellIndices[ke].subRegion]
                                                       [conn.connectedCellIndices[ke].index    ];
    }

    const localIndex stencilSize = conn.stencilCellIndices.size();
    elementLocalDofIndexCol.resize(stencilSize);
    for (localIndex ke = 0; ke < stencilSize; ++ke)
    {
      // DOF index is equal to cell index since only one DOF per cell
      elementLocalDofIndexCol[ke] = blockLocalDofNumber[conn.stencilCellIndices[ke].region   ]
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
      elementLocalDofIndexRow[0] = blockLocalDofNumber[er][esr][k];

      sparsity->InsertGlobalIndices( 1,
                                     elementLocalDofIndexRow.data(),
                                     1,
                                     elementLocalDofIndexRow.data());
    }
  });

}



void SinglePhaseFlow_TPFA::AssembleSystem ( DomainPartition * const  domain,
                                            EpetraBlockSystem * const blockSystem,
                                            real64 const time_n,
                                            real64 const dt )
{
  //***** extract data required for assembly of system *****
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  ElementRegionManager * const elemManager         = mesh->getElemManager();
  FaceManager const * const    faceManager         = mesh->getFaceManager();

  integer_array const & faceGhostRank = faceManager->getReference<integer_array>(ObjectManagerBase::viewKeyStruct::ghostRankString);

  Array2dT<localIndex> const & faceToElemRegionList     = faceManager->elementRegionList();
  Array2dT<localIndex> const & faceToElemSubRegionList  = faceManager->elementSubRegionList();
  Array2dT<localIndex> const & faceToElemList           = faceManager->elementList();

  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix(BlockIDs::fluidPressureBlock,
                                                               BlockIDs::fluidPressureBlock);
  Epetra_FEVector * const residual = blockSystem->GetResidualVector(BlockIDs::fluidPressureBlock);

  jacobian->Scale(0.0);
  residual->Scale(0.0);

  auto
  elemGhostRank = elemManager->ConstructViewAccessor<integer_array>( ObjectManagerBase::
                                                                     viewKeyStruct::
                                                                     ghostRankString );

  auto blockLocalDofNumber = elemManager->
                             ConstructViewAccessor<localIndex_array>(viewKeyStruct::
                                                                     blockLocalDofNumberString);

  auto pres      = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::fluidPressureString);
  auto dPres     = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidPressureString);
  auto dens      = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::fluidDensityString);
  auto dDens     = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidDensityString);
  auto visc      = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::fluidViscosityString);
  auto dVisc     = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidViscosityString);
  auto poro      = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::porosityString);
  auto dPoro     = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaPorosityString);
  auto volume    = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::volumeString);
  auto gravDepth = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::gravityDepthString);

  auto trans     = faceManager->getReference<real64_array>(viewKeyStruct::transmissibilityString);

  //***** Loop over all elements and assemble the change in volume/density terms *****
  Epetra_IntSerialDenseVector elemDOF(1);
  Epetra_SerialDenseVector localElemResidual(1);
  Epetra_SerialDenseMatrix localElem_dRdP(1, 1);

  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k)->void
  {
    if( elemGhostRank[er][esr][k]<0 )
    {
      elemDOF(0) = blockLocalDofNumber[er][esr][k];

      real64 const dens_new = dens[er][esr][k] + dDens[er][esr][k];
      real64 const poro_new = poro[er][esr][k] + dPoro[er][esr][k];
      real64 const vol      = volume[er][esr][k];

      // Residual contribution is mass conservation in the cell
      localElemResidual(0) = poro_new         * dens_new         * vol
                           - poro[er][esr][k] * dens[er][esr][k] * vol;

      // Derivative of residual wrt to pressure in the cell
      localElem_dRdP(0, 0) = (m_dPoro_dPres[er][esr][k] * dens_new * vol)
                           + (m_dDens_dPres[er][esr][k] * poro_new * vol);

      // add contribution to global residual and dRdP
      residual->SumIntoGlobalValues(elemDOF, localElemResidual);
      jacobian->SumIntoGlobalValues(elemDOF, localElem_dRdP);
    }
  });


  constexpr localIndex numElems = 2;
  Epetra_IntSerialDenseVector eqnRowIndices(numElems);
  Epetra_IntSerialDenseVector dofColIndices;
  Epetra_SerialDenseVector localFlux(numElems);
  Epetra_SerialDenseMatrix localFluxJacobian;

  // temporary working arrays
  real64 densWeight[numElems] = { 0.5, 0.5 };
  real64 mobility[numElems], dMobility_dP[numElems];
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
      dofColIndices[ke] = blockLocalDofNumber[conn.stencilCellIndices[ke].region   ]
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

      eqnRowIndices[ke] = blockLocalDofNumber[er][esr][ei];

      // density
      real64 const density   = dens[er][esr][ei];
      real64 const dDens_dP  = m_dDens_dPres[er][esr][ei];

      // viscosity
      real64 const viscosity = visc[er][esr][ei];
      real64 const dVisc_dP  = m_dVisc_dPres[er][esr][ei];

      // mobility
      mobility[ke]  = density / viscosity;
      dMobility_dP[ke]  = dDens_dP / viscosity - mobility[ke] / viscosity * dVisc_dP;

      // average density
      densMean += densWeight[ke] * density;
      dDensMean_dP[ke] = densWeight[ke] * dDens_dP;
    }

    //***** calculation of flux *****

    // compute potential difference MPFA-style
    real64 potDif = 0.0;
    for (localIndex ke = 0; ke < stencilSize; ++ke)
    {
      localIndex er  = conn.stencilCellIndices[ke].region;
      localIndex esr = conn.stencilCellIndices[ke].subRegion;
      localIndex ei  = conn.stencilCellIndices[ke].index;

      real64 const gravD    = gravDepth[er][esr][ei];
      real64 const gravTerm = m_gravityFlag ? densMean * gravD : 0.0;
      real64 const dGrav_dP = m_gravityFlag ? dDensMean_dP[ke] * gravD : 0.0;

      potDif += conn.stencilWeights[ke] * (pres[er][esr][ei] + dPres[er][esr][ei] + gravTerm);
      dFlux_dP[ke] = conn.stencilWeights[ke] * (1.0 + dGrav_dP);
    }

    // upwinding of fluid properties (make this an option?)
    localIndex const k_up = (potDif >= 0) ? 0 : 1;

    // compute the final flux and derivatives
    real64 const flux = mobility[k_up] * potDif;
    for (localIndex ke = 0; ke < stencilSize; ++ke)
      dFlux_dP[ke] *= mobility[k_up];
    dFlux_dP[k_up] += dMobility_dP[k_up] * potDif;

    //***** end flux terms *****

    // populate local flux vector and derivatives
    localFlux(0) =  flux;
    localFlux(1) = -flux;

    for (localIndex ke = 0; ke < stencilSize; ++ke)
    {
      localFluxJacobian(0, integer_conversion<int>(ke)) =  dFlux_dP[ke];
      localFluxJacobian(1, integer_conversion<int>(ke)) = -dFlux_dP[ke];
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

  if( verboseLevel() >= 2 )
  {
    jacobian->Print(std::cout);
    residual->Print(std::cout);
  }

}

void SinglePhaseFlow_TPFA::ApplyBoundaryConditions( DomainPartition * const domain,
                                                    systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                                    real64 const time_n,
                                                    real64 const dt )
{

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  // apply pressure boundary conditions.
  ApplyDirichletBC_implicit( elemManager,
                             time_n + dt,
                             blockSystem );


  if( verboseLevel() >= 2 )
  {
    Epetra_FECrsMatrix * const dRdP = blockSystem->GetMatrix( BlockIDs::fluidPressureBlock,
                                                              BlockIDs::fluidPressureBlock );
    Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );

    dRdP->Print(std::cout);
    residual->Print(std::cout);
  }

}

real64
SinglePhaseFlow_TPFA::
CalculateResidualNorm(systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                      DomainPartition * const domain)
{

  Epetra_FEVector const * const residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );
  Epetra_Map      const * const rowMap   = blockSystem->GetRowMap( BlockIDs::fluidPressureBlock );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  auto elemGhostRank = elemManager->ConstructViewAccessor<integer_array>( ObjectManagerBase::
                                                                          viewKeyStruct::
                                                                          ghostRankString );

  auto blockLocalDofNumber = elemManager->
      ConstructViewAccessor<localIndex_array>(viewKeyStruct::blockLocalDofNumberString);

  auto refPoro = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::referencePorosityString);
  auto volume  = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::volumeString);

  // get a view into local residual vector
  int localSizeInt;
  double* localResidual = nullptr;
  residual->ExtractView(&localResidual, &localSizeInt);

  // compute the norm of local residual scaled by cell pore volume
  real64 localResidualNorm = sumOverElemsInMesh(mesh, [&] (localIndex const er,
                                                           localIndex const esr,
                                                           localIndex const k) -> real64
  {
    if (elemGhostRank[er][esr][k] < 0)
    {
      int const lid = rowMap->LID(integer_conversion<int>(blockLocalDofNumber[er][esr][k]));
      real64 const val = localResidual[lid] / (refPoro[er][esr][k] * volume[er][esr][k]);
      return val * val;
    }
    return 0.0;
  });

  // compute global residual norm
  realT globalResidualNorm;
  MPI_Allreduce(&localResidualNorm, &globalResidualNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return sqrt(globalResidualNorm);
}


void SinglePhaseFlow_TPFA::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                                real64 const scalingFactor,
                                                DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( BlockIDs::fluidPressureBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( BlockIDs::fluidPressureBlock );

  int dummy;
  double* local_solution = nullptr;
  solution->ExtractView(&local_solution,&dummy);

  ElementRegionManager * const elementRegionManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor<localIndex_array>
  blockLocalDofNumber = elementRegionManager->
                        ConstructViewAccessor<localIndex_array>( viewKeys.blockLocalDofNumber.Key() );

  auto pres  = elementRegionManager->ConstructViewAccessor<real64_array>(viewKeyStruct::fluidPressureString);
  auto dPres = elementRegionManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidPressureString);

  auto dens  = elementRegionManager->ConstructViewAccessor<real64_array>(viewKeyStruct::fluidDensityString);
  auto dDens = elementRegionManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidDensityString);

  auto visc  = elementRegionManager->ConstructViewAccessor<real64_array>(viewKeyStruct::fluidViscosityString);
  auto dVisc = elementRegionManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidViscosityString);

  auto poro  = elementRegionManager->ConstructViewAccessor<real64_array>(viewKeyStruct::porosityString);
  auto dPoro = elementRegionManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaPorosityString);
  auto refPoro = elementRegionManager->ConstructViewAccessor<real64_array>(viewKeyStruct::referencePorosityString);

  ConstitutiveManager * const constitutiveManager =
      domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);

  auto
  constitutiveMap = elementRegionManager->
                    ConstructViewAccessor< std::pair< Array2dT<localIndex>,
                                                      Array2dT<localIndex> > >( CellBlockSubRegion::
                                                                                viewKeyStruct::
                                                                                constitutiveMapString );

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
      int const lid = rowMap->LID(integer_conversion<int>(blockLocalDofNumber[er][esr][k]));
      dPres[er][esr][k] += scalingFactor * local_solution[lid];
    }
  });


  // TODO Sync dP once element field syncing is reimplemented.
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
    ConstitutiveBase * const EOS = constitutiveManager->GetGroup<ConstitutiveBase>(matIndex1);

    // update dDens and derivatives
    real64 const new_pres = pres[er][esr][k] + dPres[er][esr][k];
    real64 new_value;

    EOS->FluidDensityUpdate(new_pres, matIndex2, new_value, m_dDens_dPres[er][esr][k]);
    dDens[er][esr][k] = new_value - dens[er][esr][k];

    EOS->FluidViscosityUpdate(new_pres, matIndex2, new_value, m_dVisc_dPres[er][esr][k]);
    dVisc[er][esr][k] = new_value - visc[er][esr][k];

    EOS->SimplePorosityUpdate(new_pres, refPoro[er][esr][k], matIndex2, new_value, m_dVisc_dPres[er][esr][k]);
    dPoro[er][esr][k] = new_value - poro[er][esr][k];
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
    localIndex numActual = 0;

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

        real64 const c2fDistance = cellToFaceVec.Normalize();

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
        real64 const ht = Dot(cellToFaceVec, faceConormal) * faceArea[kf] / c2fDistance;

        transInv += 1.0 / ht; // XXX: safeguard against div by zero?
        ++numActual;
      }
    }

    faceTrans[kf] = 1.0 * numActual / transInv; // XXX: safeguard against div by zero?

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
  m_dPoro_dPres.resize(elemManager->numRegions());
  m_dVisc_dPres.resize(elemManager->numRegions());
  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegion const * const elemRegion = elemManager->GetRegion(er);
    m_dDens_dPres[er].resize(elemRegion->numSubRegions());
    m_dPoro_dPres[er].resize(elemRegion->numSubRegions());
    m_dVisc_dPres[er].resize(elemRegion->numSubRegions());
    for( localIndex esr=0 ; esr<elemRegion->numSubRegions() ; ++esr )
    {
      CellBlockSubRegion const * const cellBlockSubRegion = elemRegion->GetSubRegion(esr);
      m_dDens_dPres[er][esr].resize(cellBlockSubRegion->size());
      m_dPoro_dPres[er][esr].resize(cellBlockSubRegion->size());
      m_dVisc_dPres[er][esr].resize(cellBlockSubRegion->size());
    }
  }
}

void SinglePhaseFlow_TPFA::SolveSystem( EpetraBlockSystem * const blockSystem,
                                        SystemSolverParameters const * const params )
{
  Epetra_FEVector * const
  solution = blockSystem->GetSolutionVector( BlockIDs::fluidPressureBlock );

  Epetra_FEVector * const
  residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );
  residual->Scale(-1.0);

  solution->Scale(0.0);

  m_linearSolverWrapper.SolveSingleBlockSystem( blockSystem,
                                                params,
                                                BlockIDs::fluidPressureBlock );

  if( verboseLevel() >= 2 )
  {
    solution->Print(std::cout);
  }

}

void SinglePhaseFlow_TPFA::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elementRegionManager = mesh->getElemManager();

  auto dPres = elementRegionManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidPressureString);
  auto dDens = elementRegionManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidDensityString);
  auto dVisc = elementRegionManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidViscosityString);
  auto dPoro = elementRegionManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaPorosityString);


  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k)->void
  {
    dPres[er][esr][k] = 0.0;
    dDens[er][esr][k] = 0.0;
    dVisc[er][esr][k] = 0.0;
    dPoro[er][esr][k] = 0.0;
  });

}


REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseFlow_TPFA, std::string const &, ManagedGroup * const )
} /* namespace ANST */
