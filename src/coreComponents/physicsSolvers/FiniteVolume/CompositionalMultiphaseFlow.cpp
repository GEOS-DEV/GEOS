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
 * @file CompositionalMultiphaseFlow.cpp
 */

#include "CompositionalMultiphaseFlow.hpp"

#include "ArrayView.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/Fluid/CompositionalMultiphaseFluid.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/BoundaryConditions/BoundaryConditionManager.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "systemSolverInterface/LinearSolverWrapper.hpp"
#include "systemSolverInterface/EpetraBlockSystem.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace systemSolverInterface;
using namespace multidimensionalArray;

CompositionalMultiphaseFlow::CompositionalMultiphaseFlow(const string & name,
                                                         ManagedGroup * const parent)
  :
  FlowSolverBase(name, parent),
  m_numPhases(0),
  m_numComponents(0)
{
  // set the blockID for the block system interface
  //getLinearSystemRepository()->SetBlockID(BlockIDs::fluidPressureBlock, this->getName());
  getLinearSystemRepository()->SetBlockID(BlockIDs::compositionalBlock, this->getName());

  this->RegisterViewWrapper( viewKeys.temperature.Key(), &m_temperature, false );
}

void CompositionalMultiphaseFlow::FillDocumentationNode()
{
  FlowSolverBase::FillDocumentationNode();

  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName(CompositionalMultiphaseFlow::CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("A compositional multiphase flow solver");

  docNode->AllocateChildNode(viewKeys.temperature.Key(),
                             viewKeys.temperature.Key(),
                             -1,
                             "real64",
                             "real64",
                             "Temperature",
                             "Fluid pressure",
                             "REQUIRED",
                             "",
                             1,
                             1,
                             0);
}

void CompositionalMultiphaseFlow::FillOtherDocumentationNodes(dataRepository::ManagedGroup * const rootGroup)
{
  FlowSolverBase::FillOtherDocumentationNodes(rootGroup);

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  for (auto & mesh : domain->getMeshBodies()->GetSubGroups())
  {
    MeshLevel * meshLevel = ManagedGroup::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);
    ElementRegionManager * const elemManager = meshLevel->getElemManager();

    elemManager->forCellBlocks([&](CellBlockSubRegion * const cellBlock) -> void
                               {
                                 cxx_utilities::DocumentationNode * const docNode = cellBlock->getDocumentationNode();

                                 docNode->AllocateChildNode(viewKeys.pressure.Key(),
                                                            viewKeys.pressure.Key(),
                                                            -1,
                                                            "real64_array",
                                                            "real64_array",
                                                            "Fluid pressure",
                                                            "Fluid pressure",
                                                            "",
                                                            elemManager->getName(),
                                                            1,
                                                            0,
                                                            0);

                                 docNode->AllocateChildNode(viewKeys.deltaPressure.Key(),
                                                            viewKeys.deltaPressure.Key(),
                                                            -1,
                                                            "real64_array",
                                                            "real64_array",
                                                            "Change in fluid pressure",
                                                            "Change in fluid pressure",
                                                            "",
                                                            elemManager->getName(),
                                                            1,
                                                            0,
                                                            1);

                                 docNode->AllocateChildNode(viewKeys.globalCompDensity.Key(),
                                                            viewKeys.globalCompDensity.Key(),
                                                            -1,
                                                            "real64_array2d",
                                                            "real64_array2d",
                                                            "Global component density in mixture",
                                                            "Global component density in mixture",
                                                            "",
                                                            elemManager->getName(),
                                                            1,
                                                            0,
                                                            0);

                                 docNode->AllocateChildNode(viewKeys.deltaGlobalCompDensity.Key(),
                                                            viewKeys.deltaGlobalCompDensity.Key(),
                                                            -1,
                                                            "real64_array2d",
                                                            "real64_array2d",
                                                            "Change in global component density in mixture",
                                                            "Change in global component density in mixture",
                                                            "",
                                                            elemManager->getName(),
                                                            1,
                                                            0,
                                                            1);

                                 docNode->AllocateChildNode(viewKeys.globalCompMoleFraction.Key(),
                                                            viewKeys.globalCompMoleFraction.Key(),
                                                            -1,
                                                            "real64_array2d",
                                                            "real64_array2d",
                                                            "Global component mole fraction in mixture",
                                                            "Global component mole fraction in mixture",
                                                            "",
                                                            elemManager->getName(),
                                                            1,
                                                            0,
                                                            0);

                                 docNode->AllocateChildNode(viewKeys.dGlobalCompMoleFraction_dGlobalCompDensity.Key(),
                                                            viewKeys.dGlobalCompMoleFraction_dGlobalCompDensity.Key(),
                                                            -1,
                                                            "real64_array3d",
                                                            "real64_array3d",
                                                            "Derivatives of global component mole fraction",
                                                            "Derivatives of global component mole fraction",
                                                            "",
                                                            elemManager->getName(),
                                                            1,
                                                            0,
                                                            3);

                                 docNode->AllocateChildNode(viewKeys.phaseVolumeFraction.Key(),
                                                            viewKeys.phaseVolumeFraction.Key(),
                                                            -1,
                                                            "real64_array2d",
                                                            "real64_array2d",
                                                            "Fluid phase volume fraction",
                                                            "Fluid phase volume fraction",
                                                            "",
                                                            elemManager->getName(),
                                                            1,
                                                            0,
                                                            3);

                                 docNode->AllocateChildNode(viewKeys.phaseDensity.Key(),
                                                            viewKeys.phaseDensity.Key(),
                                                            -1,
                                                            "real64_array2d",
                                                            "real64_array2d",
                                                            "Fluid phase density",
                                                            "Fluid phase density",
                                                            "",
                                                            elemManager->getName(),
                                                            1,
                                                            0,
                                                            3);

                                 docNode->AllocateChildNode(viewKeys.phaseComponentMassFraction.Key(),
                                                            viewKeys.phaseComponentMassFraction.Key(),
                                                            -1,
                                                            "real64_array3d",
                                                            "real64_array3d",
                                                            "Fluid component-in-phase mass fraction",
                                                            "Fluid component-in-phase mass fraction",
                                                            "",
                                                            elemManager->getName(),
                                                            1,
                                                            0,
                                                            3);

                                 docNode->AllocateChildNode(viewKeys.porosity.Key(),
                                                            viewKeys.porosity.Key(),
                                                            -1,
                                                            "real64_array",
                                                            "real64_array",
                                                            "Porosity",
                                                            "Porosity",
                                                            "",
                                                            elemManager->getName(),
                                                            1,
                                                            0,
                                                            3);

                                 docNode->AllocateChildNode(viewKeys.blockLocalDofNumber.Key(),
                                                            viewKeys.blockLocalDofNumber.Key(),
                                                            -1,
                                                            "globalIndex_array2d",
                                                            "globalIndex_array2d",
                                                            "Pressure DOF index",
                                                            "Pressure DOF index",
                                                            "0",
                                                            "",
                                                            1,
                                                            0,
                                                            3);
                               });

    {
      FaceManager * const faceManager = meshLevel->getFaceManager();
      cxx_utilities::DocumentationNode * const docNode = faceManager->getDocumentationNode();

      docNode->AllocateChildNode(viewKeys.facePressure.Key(),
                                 viewKeys.facePressure.Key(),
                                 -1,
                                 "real64_array",
                                 "real64_array",
                                 "Fluid pressure",
                                 "Fluid pressure",
                                 "",
                                 faceManager->getName(),
                                 1,
                                 0,
                                 3);

      docNode->AllocateChildNode(viewKeys.globalCompMoleFraction.Key(),
                                 viewKeys.globalCompMoleFraction.Key(),
                                 -1,
                                 "real64_array2d",
                                 "real64_array2d",
                                 "Global component mole fraction in mixture",
                                 "Global component mole fraction in mixture",
                                 "",
                                 elemManager->getName(),
                                 1,
                                 0,
                                 3);

      docNode->AllocateChildNode(viewKeys.phaseDensity.Key(),
                                 viewKeys.phaseDensity.Key(),
                                 -1,
                                 "real64_array2d",
                                 "real64_array2d",
                                 "Fluid phase density",
                                 "Fluid phase density",
                                 "",
                                 faceManager->getName(),
                                 1,
                                 0,
                                 3);

      docNode->AllocateChildNode(viewKeys.phaseComponentMassFraction.Key(),
                                 viewKeys.phaseComponentMassFraction.Key(),
                                 -1,
                                 "real64_array3d",
                                 "real64_array3d",
                                 "Fluid component-in-phase mass fraction",
                                 "Fluid component-in-phase mass fraction",
                                 "",
                                 faceManager->getName(),
                                 1,
                                 0,
                                 3);

      docNode->AllocateChildNode(viewKeys.phaseViscosity.Key(),
                                 viewKeys.phaseViscosity.Key(),
                                 -1,
                                 "real64_array2d",
                                 "real64_array2d",
                                 "Fluid phase viscosity",
                                 "Fluid phase viscosity",
                                 "",
                                 faceManager->getName(),
                                 1,
                                 0,
                                 3);

      docNode->AllocateChildNode(viewKeys.phaseRelativePermeability.Key(),
                                 viewKeys.phaseRelativePermeability.Key(),
                                 -1,
                                 "real64_array2d",
                                 "real64_array2d",
                                 "Fluid phase relative permeability",
                                 "Fluid phase relative permeability",
                                 "",
                                 faceManager->getName(),
                                 1,
                                 0,
                                 3);
    }
  }
}

void CompositionalMultiphaseFlow::InitializePreSubGroups( ManagedGroup * const rootGroup )
{
  FlowSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  ConstitutiveManager * const cm = domain->getConstitutiveManager();
  ConstitutiveBase const * fluid = cm->GetConstitituveRelation( this->m_fluidName );
  CompositionalMultiphaseFluid const * mpFluid = fluid->group_cast<CompositionalMultiphaseFluid const *>();

  m_numPhases = mpFluid->numFluidPhases();
  m_numComponents = mpFluid->numFluidComponents();

  // compute number of DOF per cell
  m_numDofPerCell = m_numComponents + 1;

  ResizeFields(domain);
}

void CompositionalMultiphaseFlow::ResizeFields(DomainPartition * domain)
{
  for (auto & mesh : domain->getMeshBodies()->GetSubGroups())
  {
    MeshLevel * meshLevel = ManagedGroup::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);
    ElementRegionManager * const elemManager = meshLevel->getElemManager();

    elemManager->forCellBlocks([&](CellBlockSubRegion * const cellBlock) -> void
    {
      cellBlock->getReference<array2d<real64>>(viewKeys.globalCompDensity).resizeDimension<1>(m_numComponents);
      cellBlock->getReference<array2d<real64>>(viewKeys.deltaGlobalCompDensity).resizeDimension<1>(m_numComponents);
      cellBlock->getReference<array2d<real64>>(viewKeys.globalCompMoleFraction).resizeDimension<1>(m_numComponents);
      cellBlock->getReference<array2d<real64>>(viewKeys.phaseVolumeFraction).resizeDimension<1>(m_numPhases);
      cellBlock->getReference<array2d<real64>>(viewKeys.phaseDensity).resizeDimension<1>(m_numPhases);
      cellBlock->getReference<array3d<real64>>(viewKeys.phaseComponentMassFraction).resizeDimension<1,2>(m_numPhases, m_numComponents);
      cellBlock->getReference<array3d<real64>>(viewKeys.dGlobalCompMoleFraction_dGlobalCompDensity).resizeDimension<1,2>(m_numComponents, m_numComponents);
      cellBlock->getReference<array2d<real64>>(viewKeys.blockLocalDofNumber).resizeDimension<1>(m_numDofPerCell);
    });

    {
      FaceManager * const faceManager = meshLevel->getFaceManager();

      faceManager->getReference<array2d<real64>>(viewKeys.phaseDensity).resizeDimension<1>(m_numPhases);
      faceManager->getReference<array2d<real64>>(viewKeys.phaseViscosity).resizeDimension<1>(m_numPhases);
      faceManager->getReference<array2d<real64>>(viewKeys.phaseRelativePermeability).resizeDimension<1>(m_numPhases);
      faceManager->getReference<array3d<real64>>(viewKeys.phaseComponentMassFraction).resizeDimension<1,2>(m_numPhases, m_numComponents);;
    }
  }
}

void CompositionalMultiphaseFlow::FinalInitialization(ManagedGroup * const rootGroup)
{
  FlowSolverBase::FinalInitialization( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  //TODO this is a hack until the sets are fixed to include ghosts!!
  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back(viewKeys.pressure.Key());
  fieldNames["elems"].push_back(viewKeys.globalCompDensity.Key());
  CommunicationTools::SynchronizeFields(fieldNames,
                                        domain->getMeshBody(0)->getMeshLevel(0),
                                        domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );
}

real64 CompositionalMultiphaseFlow::SolverStep( real64 const & time_n,
                                                real64 const & dt,
                                                integer const cycleNumber,
                                                DomainPartition * domain )
{
  // currently the only method is implicit time integration
  return this->NonlinearImplicitStep( time_n,
                                      dt,
                                      cycleNumber,
                                      domain,
                                      getLinearSystemRepository() );
}

void CompositionalMultiphaseFlow::UpdateComponentFraction(DomainPartition * domain)
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  ConstitutiveManager * const cm = domain->getConstitutiveManager();
  ConstitutiveBase const * fluid = cm->GetConstitituveRelation( this->m_fluidName );
  CompositionalMultiphaseFluid const * mpFluid = fluid->group_cast<CompositionalMultiphaseFluid const *>();

  auto molarWeight = mpFluid->getReference<array1d<real64>>( CompositionalMultiphaseFluid::
                                                             viewKeyStruct::
                                                             componentMolarWeightString );

  auto compDens = elemManager->ConstructViewAccessor<array2d<real64>>(viewKeys.globalCompDensity.Key());
  auto dCompDens = elemManager->ConstructViewAccessor<array2d<real64>>(viewKeys.deltaGlobalCompDensity.Key());

  auto compMoleFrac =
    elemManager->ConstructViewAccessor<array2d<real64>>(viewKeys.globalCompMoleFraction.Key());
  auto dCompMoleFrac_dCompDens =
    elemManager->ConstructViewAccessor<array3d<real64>>(viewKeys.dGlobalCompMoleFraction_dGlobalCompDensity.Key());

  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei )->void
  {
    real64 molarDensSum = 0.0;
    for (localIndex ic = 0; ic < m_numComponents; ++ic)
    {
      molarDensSum += (compDens[er][esr][ei][ic] + dCompDens[er][esr][ei][ic]) / molarWeight[ic];
    }

    for (localIndex ic = 0; ic < m_numComponents; ++ic)
    {
      compMoleFrac[er][esr][ei][ic] = (compDens[er][esr][ei][ic] + dCompDens[er][esr][ei][ic]) / molarWeight[ic] / molarDensSum;
      for (localIndex jc = 0; jc < m_numComponents; ++jc)
      {
        dCompMoleFrac_dCompDens[er][esr][ei][ic][jc] = - compMoleFrac[er][esr][ei][ic] / molarWeight[jc] / molarDensSum;
      }
      dCompMoleFrac_dCompDens[er][esr][ei][ic][ic] += 1.0 / molarWeight[ic] / molarDensSum;
    }
  });
}

void CompositionalMultiphaseFlow::UpdateConstitutiveModels(DomainPartition * domain)
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  auto constitutiveRelations = elemManager->ConstructConstitutiveAccessor<ConstitutiveBase>( constitutiveManager );

  auto pres         = elemManager->ConstructViewAccessor<array1d<real64>>(viewKeys.pressure.Key());
  auto dPres        = elemManager->ConstructViewAccessor<array1d<real64>>(viewKeys.deltaPressure.Key());
  auto compMoleFrac = elemManager->ConstructViewAccessor<array2d<real64>>(viewKeys.globalCompMoleFraction.Key());

  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei )->void
  {
    constitutiveRelations[er][esr][m_fluidIndex]->StateUpdatePointMultiphaseFluid( pres[er][esr][ei] + dPres[er][esr][ei],
                                                                                   m_temperature,
                                                                                   compMoleFrac[er][esr][ei].data(),
                                                                                   ei, 0 ); // fluid

    constitutiveRelations[er][esr][m_solidIndex]->StateUpdatePointPressure(pres[er][esr][ei], ei, 0); // solid
  });
}

void CompositionalMultiphaseFlow::BackupFields(DomainPartition * domain)
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  auto phaseDensOld         = elemManager->ConstructViewAccessor<array2d<real64>>(viewKeys.phaseDensity.Key());
  auto phaseVolFracOld      = elemManager->ConstructViewAccessor<array2d<real64>>(viewKeys.phaseVolumeFraction.Key());
  auto phaseCompMassFracOld = elemManager->ConstructViewAccessor<array3d<real64>>(viewKeys.phaseComponentMassFraction.Key());
  auto poroOld              = elemManager->ConstructViewAccessor<array1d<real64>>(viewKeys.porosity.Key());
  auto poroRef              = elemManager->ConstructViewAccessor<array1d<real64>>(viewKeys.referencePorosity.Key());

  auto const pvmult =
    elemManager->ConstructMaterialViewAccessor< array2d<real64> >( ConstitutiveBase::
                                                                   viewKeyStruct::
                                                                   poreVolumeMultiplierString,
                                                                   constitutiveManager );


  auto const phaseDens =
    elemManager->ConstructMaterialViewAccessor< array3d<real64> >( CompositionalMultiphaseFluid::
                                                                   viewKeyStruct::
                                                                   phaseDensityString,
                                                                   constitutiveManager );

  auto const phaseVolFrac =
    elemManager->ConstructMaterialViewAccessor< array3d<real64> >( CompositionalMultiphaseFluid::
                                                                   viewKeyStruct::
                                                                   phaseVolumeFractionString,
                                                                   constitutiveManager );

  auto const phaseCompMassFrac =
    elemManager->ConstructMaterialViewAccessor< array4d<real64> >( CompositionalMultiphaseFluid::
                                                                   viewKeyStruct::
                                                                   phaseComponentMassFractionString,
                                                                   constitutiveManager );

// backup some fields used in time derivative approximation
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei )->void
  {
    for (localIndex ip = 0; ip < m_numPhases; ++ip)
    {
      phaseDensOld[er][esr][ei][ip] = phaseDens[er][esr][m_fluidIndex][ei][0][ip];
      phaseVolFracOld[er][esr][ei][ip] = phaseVolFrac[er][esr][m_fluidIndex][ei][0][ip];
      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        phaseCompMassFracOld[er][esr][ei][ip][ic] = phaseCompMassFrac[er][esr][m_fluidIndex][ei][0][ip][ic];
      }
    }

    poroOld[er][esr][ei] = poroRef[er][esr][ei] * pvmult[er][esr][m_solidIndex][ei][0];

  });
}

void
CompositionalMultiphaseFlow::ImplicitStepSetup( real64 const & time_n, real64 const & dt,
                                                DomainPartition * const domain,
                                                EpetraBlockSystem * const blockSystem )
{
  // set deltas to zero and recompute dependent quantities
  ResetStateToBeginningOfStep( domain );

  // backup fields used in time derivative approximation
  BackupFields( domain );

  // setup dof numbers and linear system
  SetupSystem( domain, blockSystem );
}

void CompositionalMultiphaseFlow::SetNumRowsAndTrilinosIndices( MeshLevel * const meshLevel,
                                                                localIndex & numLocalRows,
                                                                globalIndex & numGlobalRows,
                                                                localIndex offset )
{
  ElementRegionManager * const elementRegionManager = meshLevel->getElemManager();

  auto blockLocalDofNumber =
    elementRegionManager->ConstructViewAccessor<array2d<globalIndex>>( viewKeys.blockLocalDofNumber.Key(), string() );

  ElementRegionManager::ElementViewAccessor< integer_array >
    ghostRank = elementRegionManager->
    ConstructViewAccessor<integer_array>( ObjectManagerBase::viewKeyStruct::ghostRankString );

  int numMpiProcesses;
  MPI_Comm_size( MPI_COMM_GEOSX, &numMpiProcesses );

  int thisMpiProcess = 0;
  MPI_Comm_rank( MPI_COMM_GEOSX, &thisMpiProcess );

  localIndex numLocalRowsToSend = numLocalRows;
  array1d<localIndex> gather(numMpiProcesses);

  // communicate the number of local rows to each process
  m_linearSolverWrapper.m_epetraComm.GatherAll( &numLocalRowsToSend,
                                                gather.data(),
                                                1 );

  GEOS_ASSERT( numLocalRows == numLocalRowsToSend, "number of local rows inconsistent" );

  // find the first local row on this partition, and find the number of total global rows.
  localIndex firstLocalRow = 0;
  numGlobalRows = 0;

  for( integer p=0 ; p<numMpiProcesses ; ++p)
  {
    numGlobalRows += gather[p];
    if (p < thisMpiProcess)
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
  raja::ReduceSum< reducePolicy, localIndex  > localCount(0);
  forAllElemsInMesh<RAJA::seq_exec>( meshLevel, [=]( localIndex const er,
                                                     localIndex const esr,
                                                     localIndex const ei ) mutable ->void
  {
    if( ghostRank[er][esr][ei] < 0 )
    {
      for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
      {
        blockLocalDofNumber[er][esr][ei][idof] = firstLocalRow + localCount + offset;
        localCount += 1;
      }
    }
  });

  GEOS_ASSERT(localCount == numLocalRows, "Number of DOF assigned does not match numLocalRows" );
}

void CompositionalMultiphaseFlow::SetSparsityPattern( DomainPartition const * const domain,
                                                      Epetra_FECrsGraph * const sparsity )
{
  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elementRegionManager = meshLevel->getElemManager();

  auto blockLocalDofNumber =
    elementRegionManager->ConstructViewAccessor<array2d<globalIndex>>( viewKeys.blockLocalDofNumber.Key() );

  auto elemGhostRank =
    elementRegionManager->ConstructViewAccessor<integer_array>( ObjectManagerBase::viewKeyStruct::ghostRankString );

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
  stencilCollection.forAll<RAJA::seq_exec>([&] (StencilCollection<CellDescriptor, real64>::Accessor stencil) -> void
  {
    elementLocalDofIndexRow.resize(numElems * m_numDofPerCell);
    stencil.forConnected([&] (CellDescriptor const & cell, localIndex const i) -> void
    {
      for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
      {
        elementLocalDofIndexRow[i * m_numDofPerCell + idof] =
          blockLocalDofNumber[cell.region][cell.subRegion][cell.index][idof];
      }
    });

    localIndex const stencilSize = stencil.size();
    elementLocalDofIndexCol.resize(stencilSize * m_numDofPerCell);
    stencil.forAll([&] (CellDescriptor const & cell, real64 w, localIndex const i) -> void
    {
      for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
      {
        elementLocalDofIndexCol[i * m_numDofPerCell + idof] =
          blockLocalDofNumber[cell.region][cell.subRegion][cell.index][idof];
      }
    });

    sparsity->InsertGlobalIndices( integer_conversion<int>(numElems * m_numDofPerCell),
                                   elementLocalDofIndexRow.data(),
                                   integer_conversion<int>(stencilSize * m_numDofPerCell),
                                   elementLocalDofIndexCol.data() );
  });

  elementLocalDofIndexRow.resize(m_numDofPerCell);

  // loop over all elements and add all locals just in case the above connector loop missed some
  forAllElemsInMesh<RAJA::seq_exec>(meshLevel, [&] (localIndex const er,
                                                    localIndex const esr,
                                                    localIndex const ei) -> void
  {
    if (elemGhostRank[er][esr][ei] < 0)
    {
      for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
      {
        elementLocalDofIndexRow[idof] = blockLocalDofNumber[er][esr][ei][idof];
      }

      sparsity->InsertGlobalIndices( integer_conversion<int>(m_numDofPerCell),
                                     elementLocalDofIndexRow.data(),
                                     integer_conversion<int>(m_numDofPerCell),
                                     elementLocalDofIndexRow.data());
    }
  });

  // add additional connectivity resulting from boundary stencils
  fluxApprox->forBoundaryStencils([&] (FluxApproximationBase::BoundaryStencil const & boundaryStencilCollection) -> void
  {
    boundaryStencilCollection.forAll<RAJA::seq_exec>([=] (StencilCollection<PointDescriptor, real64>::Accessor stencil) mutable -> void
    {
      stencil.forConnected([&] (PointDescriptor const & point, localIndex const i) -> void
      {
        if (point.tag == PointDescriptor::Tag::CELL)
        {
          CellDescriptor const & c = point.cellIndex;
          for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
          {
            elementLocalDofIndexRow[idof] = blockLocalDofNumber[c.region][c.subRegion][c.index][idof];
          }
        }
      });

      localIndex const stencilSize = stencil.size();
      elementLocalDofIndexCol.resize(stencilSize * m_numDofPerCell);
      integer counter = 0;
      stencil.forAll([&] (PointDescriptor const & point, real64 w, localIndex i) -> void
      {
        if (point.tag == PointDescriptor::Tag::CELL)
        {
          CellDescriptor const & c = point.cellIndex;
          for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
          {
            elementLocalDofIndexCol[counter * m_numDofPerCell + idof] =
              blockLocalDofNumber[c.region][c.subRegion][c.index][idof];
          }
          ++counter;
        }
      });

      sparsity->InsertGlobalIndices( integer_conversion<int>(m_numDofPerCell),
                                     elementLocalDofIndexRow.data(),
                                     integer_conversion<int>(counter * m_numDofPerCell),
                                     elementLocalDofIndexCol.data() );
    });
  });
}

void CompositionalMultiphaseFlow::SetupSystem( DomainPartition * const domain,
                                               EpetraBlockSystem * const blockSystem )
{
  // assume that there is only a single MeshLevel for now
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elementRegionManager = mesh->getElemManager();

  // for this solver, the dof are on the cell center, and the block of rows corresponds to a cell
  localIndex numLocalRows  = 0;
  globalIndex numGlobalRows = 0;

  // get the number of local elements, and ghost elements...i.e. local rows and ghost rows
  elementRegionManager->forCellBlocks( [&]( CellBlockSubRegion * const subRegion )
  {
    numLocalRows += subRegion->size() - subRegion->GetNumberOfGhosts();
  });
  numLocalRows *= m_numDofPerCell;

  localIndex_array displacementIndices;
  SetNumRowsAndTrilinosIndices( mesh,
                                numLocalRows,
                                numGlobalRows,
                                0 );

  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back(viewKeys.blockLocalDofNumber.Key());
  CommunicationTools::
  SynchronizeFields(fieldNames,
                    mesh,
                    domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );


  // construct row map, and set a pointer to the row map
  Epetra_Map * const
    rowMap = blockSystem->
      SetRowMap( BlockIDs::compositionalBlock,
                 std::make_unique<Epetra_Map>( numGlobalRows,
                                               numLocalRows,
                                               0,
                                               m_linearSolverWrapper.m_epetraComm ) );

  // construct sparsity matrix, set a pointer to the sparsity pattern matrix
  Epetra_FECrsGraph * const
    sparsity = blockSystem->SetSparsity( BlockIDs::compositionalBlock,
                                         BlockIDs::compositionalBlock,
                                         std::make_unique<Epetra_FECrsGraph>(Copy,*rowMap,0) );



  // set the sparsity patter
  SetSparsityPattern( domain, sparsity );

  // assemble the global sparsity matrix
  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();

  // construct system matrix
  blockSystem->SetMatrix( BlockIDs::compositionalBlock,
                          BlockIDs::compositionalBlock,
                          std::make_unique<Epetra_FECrsMatrix>(Copy,*sparsity) );

  // construct solution vector
  blockSystem->SetSolutionVector( BlockIDs::compositionalBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

  // construct residual vector
  blockSystem->SetResidualVector( BlockIDs::compositionalBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );
}

void CompositionalMultiphaseFlow::AssembleSystem( DomainPartition * const domain,
                                                  EpetraBlockSystem * const blockSystem,
                                                  real64 const time_n, real64 const dt )
{

}

void CompositionalMultiphaseFlow::ApplyBoundaryConditions( DomainPartition * const domain,
                                                           EpetraBlockSystem * const blockSystem,
                                                           real64 const time_n, real64 const dt )
{
  // apply pressure boundary conditions.
  ApplyDirichletBC_implicit(domain, time_n, dt, blockSystem);
  ApplyFaceDirichletBC_implicit(domain, time_n, dt, blockSystem);

  if (verboseLevel() >= 2)
  {
    blockSystem->GetMatrix( BlockIDs::compositionalBlock, BlockIDs::compositionalBlock )->Print( std::cout );
    blockSystem->GetResidualVector( BlockIDs::compositionalBlock )->Print( std::cout );
  }
}

void
CompositionalMultiphaseFlow::ApplyDirichletBC_implicit( DomainPartition * object,
                                                        real64 const time_n, real64 const dt,
                                                        EpetraBlockSystem * const blockSystem )
{

}

void
CompositionalMultiphaseFlow::ApplyFaceDirichletBC_implicit( DomainPartition * domain,
                                                            real64 const time_n, real64 const dt,
                                                            EpetraBlockSystem * const blockSystem )
{

}

real64
CompositionalMultiphaseFlow::CalculateResidualNorm( EpetraBlockSystem const * const blockSystem,
                                                    DomainPartition * const domain )
{
  Epetra_FEVector const * const residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );
  Epetra_Map      const * const rowMap   = blockSystem->GetRowMap( BlockIDs::fluidPressureBlock );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  auto elemGhostRank = elemManager->ConstructViewAccessor<integer_array>( ObjectManagerBase::
                                                                          viewKeyStruct::
                                                                          ghostRankString );

  auto blockLocalDofNumber = elemManager->
    ConstructViewAccessor<array2d<globalIndex>>(viewKeys.blockLocalDofNumber.Key());

  auto refPoro = elemManager->ConstructViewAccessor<real64_array>(viewKeys.referencePorosity.Key());
  auto volume  = elemManager->ConstructViewAccessor<real64_array>(CellBlock::viewKeyStruct::elementVolumeString);

  // get a view into local residual vector
  int localSizeInt;
  double* localResidual = nullptr;
  residual->ExtractView(&localResidual, &localSizeInt);

  // compute the norm of local residual scaled by cell pore volume
  real64 localResidualNorm = sumOverElemsInMesh(mesh, [&] ( localIndex const er,
                                                            localIndex const esr,
                                                            localIndex const ei ) -> real64
  {
    if (elemGhostRank[er][esr][ei] < 0)
    {
      real64 cell_norm = 0.0;
      for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
      {
        int const lid = rowMap->LID(integer_conversion<int>(blockLocalDofNumber[er][esr][ei][idof]));
        real64 const val = localResidual[lid] / (refPoro[er][esr][ei] * volume[er][esr][ei]);
        cell_norm += val * val;
      }
      return cell_norm;
    }
    return 0.0;
  });

  // compute global residual norm
  realT globalResidualNorm;
  MPI_Allreduce(&localResidualNorm, &globalResidualNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_GEOSX);

  return sqrt(globalResidualNorm);
}

void CompositionalMultiphaseFlow::SolveSystem( EpetraBlockSystem * const blockSystem,
                                               SystemSolverParameters const * const params )
{
  Epetra_FEVector * const
    solution = blockSystem->GetSolutionVector( BlockIDs::compositionalBlock );

  Epetra_FEVector * const
    residual = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );

  residual->Scale(-1.0);
  solution->Scale(0.0);

  m_linearSolverWrapper.SolveSingleBlockSystem( blockSystem,
                                                params,
                                                BlockIDs::compositionalBlock );

  if( verboseLevel() >= 2 )
  {
    solution->Print(std::cout);
  }
}

void
CompositionalMultiphaseFlow::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                                  real64 const scalingFactor,
                                                  DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( BlockIDs::compositionalBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( BlockIDs::compositionalBlock );

  int dummy;
  double* local_solution = nullptr;
  solution->ExtractView(&local_solution,&dummy);

  auto blockLocalDofNumber =
    elemManager->ConstructViewAccessor<array2d<globalIndex>>( viewKeys.blockLocalDofNumber.Key() );

  auto dPres     = elemManager->ConstructViewAccessor<array1d<real64>>( viewKeys.deltaPressure.Key() );
  auto dCompDens = elemManager->ConstructViewAccessor<array1d<real64>>( viewKeys.deltaGlobalCompDensity.Key() );

  auto elemGhostRank = elemManager->ConstructViewAccessor<integer_array>( ObjectManagerBase::
                                                                          viewKeyStruct::
                                                                          ghostRankString );

  // loop over all elements to update incremental pressure
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei )->void
  {
    if( elemGhostRank[er][esr][ei] < 0 )
    {
      // extract solution and apply to dP
      {
        int const lid = rowMap->LID(integer_conversion<int>(blockLocalDofNumber[er][esr][ei][0]));
        dPres[er][esr][ei] += scalingFactor * local_solution[lid];
      }

      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        int const lid = rowMap->LID(integer_conversion<int>(blockLocalDofNumber[er][esr][ei][ic+1]));
        dCompDens[er][esr][ei] += scalingFactor * local_solution[lid];
      }
    }
  });

  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back(viewKeys.deltaPressure.Key());
  fieldNames["elems"].push_back(viewKeys.deltaGlobalCompDensity.Key());
  CommunicationTools::SynchronizeFields(fieldNames,
                                        mesh,
                                        domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  UpdateComponentFraction(domain);
  UpdateConstitutiveModels(domain);
}

void CompositionalMultiphaseFlow::ResetStateToBeginningOfStep(DomainPartition * const domain)
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  auto dPres     = elemManager->ConstructViewAccessor<array1d<real64>>(viewKeys.deltaPressure.Key());
  auto dCompDens = elemManager->ConstructViewAccessor<array2d<real64>>(viewKeys.deltaGlobalCompDensity.Key());

  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei )->void
  {
    dPres[er][esr][ei] = 0.0;
    for (localIndex ic = 0; ic < m_numComponents; ++ic)
      dCompDens[er][esr][ei][ic] = 0.0;
  });

  UpdateComponentFraction(domain);
  UpdateConstitutiveModels(domain);
}

void CompositionalMultiphaseFlow::ImplicitStepComplete( real64 const & time,
                                                        real64 const & dt,
                                                        DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  auto pres      = elemManager->ConstructViewAccessor<array1d<real64>>(viewKeys.pressure.Key());
  auto dPres     = elemManager->ConstructViewAccessor<array1d<real64>>(viewKeys.deltaPressure.Key());
  auto compDens  = elemManager->ConstructViewAccessor<array2d<real64>>(viewKeys.globalCompDensity.Key());
  auto dCompDens = elemManager->ConstructViewAccessor<array2d<real64>>(viewKeys.deltaGlobalCompDensity.Key());

  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei )->void
  {
    pres[er][esr][ei] += dPres[er][esr][ei];
    for (localIndex ic = 0; ic < m_numComponents; ++ic)
      compDens[er][esr][ei][ic] += dCompDens[er][esr][ei][ic];
  });
}


REGISTER_CATALOG_ENTRY(SolverBase, CompositionalMultiphaseFlow, string const &, ManagedGroup * const)
}// namespace geosx
