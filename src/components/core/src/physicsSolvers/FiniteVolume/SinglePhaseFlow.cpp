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

#include "SinglePhaseFlow.hpp"

#include <vector>
#include <cmath>

#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
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


SinglePhaseFlow::SinglePhaseFlow( const std::string& name,
                                            ManagedGroup * const parent ):
  SolverBase(name, parent),
  m_precomputeDone(false),
  m_gravityFlag(1),
  m_dDens_dPres(),
  m_dPoro_dPres(),
  m_dVisc_dPres()
{
  // set the blockID for the block system interface
  getLinearSystemRepository()->
  SetBlockID( BlockIDs::fluidPressureBlock, this->getName() );

  this->RegisterViewWrapper(viewKeyStruct::gravityFlagString, &m_gravityFlag, false);
  this->RegisterViewWrapper(viewKeyStruct::discretizationString, &m_discretizationName, false);
}


void SinglePhaseFlow::FillDocumentationNode(  )
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

  docNode->AllocateChildNode( viewKeyStruct::discretizationString,
                              viewKeyStruct::discretizationString,
                              -1,
                              "string",
                              "string",
                              "Name of the finite volume discretization to use",
                              "Name of the finite volume discretization to use",
                              "",
                              "",
                              1,
                              1,
                              0 );
}

void SinglePhaseFlow::FillOtherDocumentationNodes( dataRepository::ManagedGroup * const rootGroup )
{
  SolverBase::FillOtherDocumentationNodes( rootGroup );
  DomainPartition * domain  = rootGroup->GetGroup<DomainPartition>(keys::domain);

  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel * meshLevel = ManagedGroup::group_cast<MeshBody*>(mesh.second)->getMeshLevel(0);
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
                                    1 );

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
                                    1 );

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
                                    1 );

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
                                    1 );

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
                                    1 );

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
                                    1 );

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
                                    1 );

        docNode->AllocateChildNode( viewKeyStruct::blockLocalDofNumberString,
                                    viewKeyStruct::blockLocalDofNumberString,
                                    -1,
                                    "globalIndex_array",
                                    "globalIndex_array",
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

void SinglePhaseFlow::FinalInitialization( ManagedGroup * const problemManager )
{
  DomainPartition * domain = problemManager->GetGroup<DomainPartition>(keys::domain);

  // Precompute solver-specific constant data (e.g. gravity-depth)
  PrecomputeData(domain);

  // Allocate additional storage for derivatives
  AllocateAuxStorage(domain);
}

real64 SinglePhaseFlow::SolverStep( real64 const& time_n,
                                         real64 const& dt,
                                         const int cycleNumber,
                                         DomainPartition * domain )
{
  if (!m_precomputeDone)
  {
    NumericalMethodsManager const * numericalMethodManager = domain->
      getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);

    FiniteVolumeManager const * fvManager = numericalMethodManager->
      GetGroup<FiniteVolumeManager>(keys::finiteVolumeManager);

    FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation(m_discretizationName);
    FluxApproximationBase::CellStencil const & stencilCollection = fluxApprox->getStencil();

    // TODO HACK, should be a separate init stage
    const_cast<FluxApproximationBase *>(fluxApprox)->compute(domain->group_cast<DomainPartition *>());
    m_precomputeDone = true;
  }

  // currently the only method is implcit time integration
  return this->NonlinearImplicitStep( time_n,
                                      dt,
                                      cycleNumber,
                                      domain->group_cast<DomainPartition*>(),
                                      getLinearSystemRepository() );
}


void SinglePhaseFlow::
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


void SinglePhaseFlow::ImplicitStepComplete( real64 const & time_n,
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
  });

}

void SinglePhaseFlow::SetNumRowsAndTrilinosIndices( MeshLevel * const meshLevel,
                                                         localIndex & numLocalRows,
                                                         localIndex & numGlobalRows,
                                                         localIndex_array& localIndices,
                                                         localIndex offset )
{

  ElementRegionManager * const elementRegionManager = meshLevel->getElemManager();
  ElementRegionManager::ElementViewAccessor<globalIndex_array>
  blockLocalDofNumber = elementRegionManager->
                  ConstructViewAccessor<globalIndex_array>( viewKeys.blockLocalDofNumber.Key(),
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

void SinglePhaseFlow::SetupSystem ( DomainPartition * const domain,
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
                      std::make_unique<Epetra_Map>( static_cast<long long>(numGlobalRows),
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


void SinglePhaseFlow::SetSparsityPattern( DomainPartition const * const domain,
                                               Epetra_FECrsGraph * const sparsity )
{
  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elementRegionManager = meshLevel->getElemManager();
  ElementRegionManager::ElementViewAccessor<globalIndex_array const>
  blockLocalDofNumber = elementRegionManager->
                  ConstructViewAccessor<globalIndex_array>( viewKeys.blockLocalDofNumber.Key() );

  ElementRegionManager::ElementViewAccessor< integer_array const >
  elemGhostRank = elementRegionManager->
              ConstructViewAccessor<integer_array>( ObjectManagerBase::viewKeyStruct::ghostRankString );

  NumericalMethodsManager const * numericalMethodManager = domain->
    getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);

  FiniteVolumeManager const * fvManager = numericalMethodManager->
    GetGroup<FiniteVolumeManager>(keys::finiteVolumeManager);

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation(m_discretizationName);
  FluxApproximationBase::CellStencil const & stencilCollection = fluxApprox->getStencil();

  globalIndex_array elementLocalDofIndexRow;
  globalIndex_array elementLocalDofIndexCol;

  //**** loop over all faces. Fill in sparsity for all pairs of DOF/elem that are connected by face
  constexpr localIndex numElems = 2;
  stencilCollection.forAll([&] (auto stencil) -> void
  {
    elementLocalDofIndexRow.resize(numElems);
    stencil.forConnected([&] (auto const cell, localIndex const i) -> void
    {
      elementLocalDofIndexRow[i] = blockLocalDofNumber[cell.region][cell.subRegion][cell.index];
    });

    localIndex const stencilSize = stencil.size();
    elementLocalDofIndexCol.resize(stencilSize);
    stencil.forAll([&] (auto const cell, real64 w, localIndex const i) -> void
    {
      elementLocalDofIndexCol[i] = blockLocalDofNumber[cell.region][cell.subRegion][cell.index];
    });

    sparsity->InsertGlobalIndices(integer_conversion<int>(numElems),
                                  elementLocalDofIndexRow.data(),
                                  integer_conversion<int>(stencilSize),
                                  elementLocalDofIndexCol.data());
  });

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

void SinglePhaseFlow::AssembleSystem(DomainPartition * const  domain,
                                     EpetraBlockSystem * const blockSystem,
                                     real64 const time_n,
                                     real64 const dt)
{
  //***** extract data required for assembly of system *****
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  ElementRegionManager * const elemManager = mesh->getElemManager();

  NumericalMethodsManager const * numericalMethodManager = domain->
    getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);

  FiniteVolumeManager const * fvManager = numericalMethodManager->
    GetGroup<FiniteVolumeManager>(keys::finiteVolumeManager);

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation(m_discretizationName);
  FluxApproximationBase::CellStencil const & stencilCollection = fluxApprox->getStencil();

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
                             ConstructViewAccessor<globalIndex_array>(viewKeyStruct::
                                                                      blockLocalDofNumberString);

  auto pres      = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::fluidPressureString);
  auto dPres     = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidPressureString);
  auto dens      = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::fluidDensityString);
  auto dDens     = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidDensityString);
  auto visc      = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::fluidViscosityString);
  auto dVisc     = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidViscosityString);
  auto poro      = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::porosityString);
  auto dPoro     = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaPorosityString);
  auto gravDepth = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::gravityDepthString);

  auto volume    = elemManager->ConstructViewAccessor<real64_array>(CellBlock::viewKeyStruct::elementVolumeString);

  //***** Loop over all elements and assemble the change in volume/density terms *****
  globalIndex elemDOF;
  real64 localAccum;
  real64 localAccumJacobian;

  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k)->void
  {
    if( elemGhostRank[er][esr][k]<0 )
    {
      elemDOF = blockLocalDofNumber[er][esr][k];

      real64 const dens_new = dens[er][esr][k] + dDens[er][esr][k];
      real64 const poro_new = poro[er][esr][k] + dPoro[er][esr][k];
      real64 const vol      = volume[er][esr][k];

      // Residual contribution is mass conservation in the cell
      localAccum = poro_new         * dens_new         * vol
                 - poro[er][esr][k] * dens[er][esr][k] * vol;

      // Derivative of residual wrt to pressure in the cell
      localAccumJacobian = (m_dPoro_dPres[er][esr][k] * dens_new * vol)
                         + (m_dDens_dPres[er][esr][k] * poro_new * vol);

      // add contribution to global residual and dRdP
      residual->SumIntoGlobalValues(1, &elemDOF, &localAccum);
      jacobian->SumIntoGlobalValues(1, &elemDOF, 1, &elemDOF, &localAccumJacobian);
    }
  });


  constexpr localIndex numElems = 2;
  globalIndex eqnRowIndices[numElems];
  globalIndex_array dofColIndices;
  real64 localFlux[numElems];
  Array2dT<real64> localFluxJacobian;

  // temporary working arrays
  real64 densWeight[numElems] = { 0.5, 0.5 };
  real64 mobility[numElems], dMobility_dP[numElems];
  real64_array dDensMean_dP, dFlux_dP;

  stencilCollection.forAll([&] (auto stencil) -> void
  {
    const localIndex stencilSize = stencil.size();

    // resize and clear local working arrays
    dDensMean_dP.resize(stencilSize); // doesn't need to be that large, but it's convenient
    dFlux_dP.resize(stencilSize);

    // clear working arrays
    dDensMean_dP = 0.0;

    // resize local matrices and vectors
    dofColIndices.resize(stencilSize);
    localFluxJacobian.resize(numElems, stencilSize);

    // calculate quantities on primary connected cells
    real64 densMean = 0.0;
    stencil.forConnected([&] (auto cell, localIndex i) -> void
    {
      localIndex const er  = cell.region;
      localIndex const esr = cell.subRegion;
      localIndex const ei  = cell.index;

      eqnRowIndices[i] = blockLocalDofNumber[er][esr][ei];

      // density
      real64 const density   = dens[er][esr][ei];
      real64 const dDens_dP  = m_dDens_dPres[er][esr][ei];

      // viscosity
      real64 const viscosity = visc[er][esr][ei];
      real64 const dVisc_dP  = m_dVisc_dPres[er][esr][ei];

      // mobility
      mobility[i]  = density / viscosity;
      dMobility_dP[i]  = dDens_dP / viscosity - mobility[i] / viscosity * dVisc_dP;

      // average density
      densMean += densWeight[i] * density;
      dDensMean_dP[i] = densWeight[i] * dDens_dP;
    });

    //***** calculation of flux *****

    // compute potential difference MPFA-style
    real64 potDif = 0.0;
    stencil.forAll([&] (auto cell, auto w, localIndex i) -> void
    {
      localIndex const er  = cell.region;
      localIndex const esr = cell.subRegion;
      localIndex const ei  = cell.index;

      dofColIndices[i] = blockLocalDofNumber[er][esr][ei];

      real64 const gravD    = gravDepth[er][esr][ei];
      real64 const gravTerm = m_gravityFlag ? densMean * gravD : 0.0;
      real64 const dGrav_dP = m_gravityFlag ? dDensMean_dP[i] * gravD : 0.0;

      potDif += w * (pres[er][esr][ei] + dPres[er][esr][ei] + gravTerm);
      dFlux_dP[i] = w * (1.0 + dGrav_dP);
    });

    // upwinding of fluid properties (make this an option?)
    localIndex const k_up = (potDif >= 0) ? 0 : 1;

    // compute the final flux and derivatives
    real64 const flux = mobility[k_up] * potDif;
    for (localIndex ke = 0; ke < stencilSize; ++ke)
      dFlux_dP[ke] *= mobility[k_up];
    dFlux_dP[k_up] += dMobility_dP[k_up] * potDif;

    //***** end flux terms *****

    // populate local flux vector and derivatives
    localFlux[0] =  dt * flux;
    localFlux[1] = -localFlux[0];

    for (localIndex ke = 0; ke < stencilSize; ++ke)
    {
      localFluxJacobian[0][ke] =   dt * dFlux_dP[ke];
      localFluxJacobian[1][ke] = - dt * dFlux_dP[ke];
    }

    // Add to global residual/jacobian
    jacobian->SumIntoGlobalValues(2, eqnRowIndices, 2, dofColIndices.data(), localFluxJacobian.data());
    residual->SumIntoGlobalValues(2, eqnRowIndices, localFlux);
  });

  jacobian->GlobalAssemble(true);
  residual->GlobalAssemble();

  if( verboseLevel() >= 2 )
  {
    jacobian->Print(std::cout);
    residual->Print(std::cout);
  }

}



void SinglePhaseFlow::ApplyBoundaryConditions(DomainPartition * const domain,
                                              systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                              real64 const time_n,
                                              real64 const dt)
{

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  // apply pressure boundary conditions.
  ApplyDirichletBC_implicit(domain, time_n + dt, blockSystem);
  ApplyFaceDirichletBC_implicit(domain, time_n + dt, blockSystem);

  if (verboseLevel() >= 2)
  {
    Epetra_FECrsMatrix * const dRdP = blockSystem->GetMatrix( BlockIDs::fluidPressureBlock,
                                                              BlockIDs::fluidPressureBlock );
    Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );

    dRdP->Print(std::cout);
    residual->Print(std::cout);
  }

}

/**
 * This function currently applies Dirichlet boundary conditions on the elements/zones as they
 * hold the DOF. Futher work will need to be done to apply a Dirichlet bc to the connectors (faces)
 */
void SinglePhaseFlow::ApplyDirichletBC_implicit( DomainPartition * domain,
                                                 real64 const time,
                                                 EpetraBlockSystem * const blockSystem )
{
  BoundaryConditionManager * bcManager = BoundaryConditionManager::get();
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor<globalIndex_array>
    blockLocalDofNumber = elemManager->
    ConstructViewAccessor<globalIndex_array>( viewKeyStruct::blockLocalDofNumberString );

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
                                                            subRegion,
                                                            blockLocalDofNumber[er][esr].get(),
                                                            1,
                                                            blockSystem,
                                                            BlockIDs::fluidPressureBlock,
                                                            [&] (localIndex const a) -> real64
        {
          return pres[er][esr][a] + dPres[er][esr][a];
        });
      });
    }
  }
}

void SinglePhaseFlow::ApplyFaceDirichletBC_implicit(DomainPartition * domain, real64 const time,
                                                    systemSolverInterface::EpetraBlockSystem * const blockSystem)
{
  BoundaryConditionManager * bcManager = BoundaryConditionManager::get();
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  FaceManager * const faceManager = mesh->getFaceManager();

  NumericalMethodsManager * const numericalMethodManager = domain->
    getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);

  FiniteVolumeManager * const fvManager = numericalMethodManager->
    GetGroup<FiniteVolumeManager>(keys::finiteVolumeManager);

  FluxApproximationBase const * const fluxApprox = fvManager->getFluxApproximation(m_discretizationName);

  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix(BlockIDs::fluidPressureBlock,
                                                               BlockIDs::fluidPressureBlock);
  Epetra_FEVector * const residual = blockSystem->GetResidualVector(BlockIDs::fluidPressureBlock);

  auto
    elemGhostRank = elemManager->ConstructViewAccessor<integer_array>( ObjectManagerBase::
                                                                       viewKeyStruct::
                                                                       ghostRankString );

  auto blockLocalDofNumber = elemManager->
    ConstructViewAccessor<globalIndex_array>(viewKeyStruct::
                                             blockLocalDofNumberString);

  auto pres      = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::fluidPressureString);
  auto dPres     = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidPressureString);
  auto dens      = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::fluidDensityString);
  auto dDens     = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidDensityString);
  auto visc      = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::fluidViscosityString);
  auto dVisc     = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaFluidViscosityString);
  auto gravDepth = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::gravityDepthString);


}


real64
SinglePhaseFlow::
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
      ConstructViewAccessor<globalIndex_array>(viewKeyStruct::blockLocalDofNumberString);

  auto refPoro = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::referencePorosityString);
  auto volume  = elemManager->ConstructViewAccessor<real64_array>(CellBlock::viewKeyStruct::elementVolumeString);

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


void SinglePhaseFlow::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
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

  ElementRegionManager::ElementViewAccessor<globalIndex_array>
  blockLocalDofNumber = elementRegionManager->
                        ConstructViewAccessor<globalIndex_array>( viewKeys.blockLocalDofNumber.Key() );

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

void SinglePhaseFlow::PrecomputeData(DomainPartition *const domain)
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  R1Tensor const & gravityVector = getGravityVector();

  auto elemCenter = elemManager->ConstructViewAccessor< r1_array >( CellBlock::
                                                                    viewKeyStruct::
                                                                    elementCenterString );

  auto gravityDepth = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::gravityDepthString);

  // Loop over all the elements and calculate element centers, and element volumes
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k )->void
  {
    gravityDepth[er][esr][k] = Dot(elemCenter[er][esr][k], gravityVector);
  });
}

void SinglePhaseFlow::AllocateAuxStorage(DomainPartition *const domain)
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

void SinglePhaseFlow::SolveSystem( EpetraBlockSystem * const blockSystem,
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

void SinglePhaseFlow::ResetStateToBeginningOfStep( DomainPartition * const domain )
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


REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseFlow, std::string const &, ManagedGroup * const )
} /* namespace ANST */
