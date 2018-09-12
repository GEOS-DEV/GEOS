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
  getLinearSystemRepository()->SetBlockID(BlockIDs::fluidPressureBlock, this->getName());
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
                                                            "DOF index",
                                                            "DOF index",
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

  resizeFields( domain );
}

void CompositionalMultiphaseFlow::resizeFields( DomainPartition * domain )
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

void CompositionalMultiphaseFlow::updateComponentFraction(DomainPartition * domain)
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

void
CompositionalMultiphaseFlow::ImplicitStepSetup( real64 const & time_n, real64 const & dt,
                                                DomainPartition * const domain,
                                                EpetraBlockSystem * const blockSystem )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  auto pres = elemManager->ConstructViewAccessor<array1d<real64>>(viewKeys.pressure.Key());
  auto dPres = elemManager->ConstructViewAccessor<array1d<real64>>(viewKeys.deltaPressure.Key());

  auto compDens = elemManager->ConstructViewAccessor<array2d<real64>>(viewKeys.globalCompDensity.Key());
  auto dCompDens = elemManager->ConstructViewAccessor<array2d<real64>>(viewKeys.deltaGlobalCompDensity.Key());

  auto compMoleFrac =
    elemManager->ConstructViewAccessor<array2d<real64>>(viewKeys.globalCompMoleFraction.Key());
  auto dCompMoleFrac_dCompDens =
    elemManager->ConstructViewAccessor<array3d<real64>>(viewKeys.dGlobalCompMoleFraction_dGlobalCompDensity.Key());

  auto phaseDensOld = elemManager->ConstructViewAccessor<array2d<real64>>(viewKeys.phaseDensity.Key());
  auto phaseVolFracOld = elemManager->ConstructViewAccessor<array2d<real64>>(viewKeys.phaseVolumeFraction.Key());
  auto phaseCompMassFracOld = elemManager->ConstructViewAccessor<array3d<real64>>(viewKeys.phaseComponentMassFraction.Key());
  auto poroOld = elemManager->ConstructViewAccessor<array1d<real64>>(viewKeys.porosity.Key());
  auto poroRef = elemManager->ConstructViewAccessor<array1d<real64>>(viewKeys.referencePorosity.Key());

  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase>
    constitutiveRelations = elemManager->ConstructConstitutiveAccessor<ConstitutiveBase>( constitutiveManager );


  ElementRegionManager::MaterialViewAccessor< array2d<real64> >
    pvmult = elemManager->ConstructMaterialViewAccessor< array2d<real64> >( ConstitutiveBase::
                                                                            viewKeyStruct::
                                                                            poreVolumeMultiplierString,
                                                                            constitutiveManager );


  ElementRegionManager::MaterialViewAccessor< array3d<real64> > const
    phaseDens = elemManager->ConstructMaterialViewAccessor< array3d<real64> >( CompositionalMultiphaseFluid::
                                                                               viewKeyStruct::
                                                                               phaseDensityString,
                                                                               constitutiveManager);

  ElementRegionManager::MaterialViewAccessor< array3d<real64> > const
    phaseVolFrac = elemManager->ConstructMaterialViewAccessor< array3d<real64> >( CompositionalMultiphaseFluid::
                                                                                  viewKeyStruct::
                                                                                  phaseVolumeFractionString,
                                                                                  constitutiveManager);

  ElementRegionManager::MaterialViewAccessor< array4d<real64> > const
    phaseCompMassFrac = elemManager->ConstructMaterialViewAccessor< array4d<real64> >( CompositionalMultiphaseFluid::
                                                                                       viewKeyStruct::
                                                                                       phaseComponentMassFractionString,
                                                                                       constitutiveManager);

  //***** loop over all elements and initialize the derivative arrays *****
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei )->void
  {
    dPres[er][esr][ei] = 0.0;
    for (localIndex ic = 0; ic < m_numComponents; ++ic)
      dCompDens[er][esr][ei][ic] = 0.0;

    updateComponentFraction( domain );

    constitutiveRelations[er][esr][m_fluidIndex]->StateUpdatePointMultiphaseFluid( pres[er][esr][ei],
                                                                                   m_temperature,
                                                                                   compMoleFrac[er][esr][ei].data(),
                                                                                   ei, 0 ); // fluid

    constitutiveRelations[er][esr][m_solidIndex]->StateUpdatePointPressure(pres[er][esr][ei], ei, 0); // solid

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

  // setup dof numbers and linear system
  SetupSystem( domain, blockSystem );
}

void CompositionalMultiphaseFlow::SetupSystem( DomainPartition * const domain,
                                               EpetraBlockSystem * const blockSystem )
{

}

void CompositionalMultiphaseFlow::SetSparsityPattern( DomainPartition const * const domain,
                                                      Epetra_FECrsGraph * const sparsity )
{

}

void CompositionalMultiphaseFlow::SetNumRowsAndTrilinosIndices( MeshLevel * const meshLevel,
                                                                localIndex & numLocalRows,
                                                                globalIndex & numGlobalRows,
                                                                localIndex_array & localIndices,
                                                                localIndex offset )
{

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
    blockSystem->GetMatrix( BlockIDs::fluidPressureBlock, BlockIDs::fluidPressureBlock )->Print( std::cout );
    blockSystem->GetMatrix( BlockIDs::fluidPressureBlock, BlockIDs::compositionalBlock )->Print( std::cout );
    blockSystem->GetMatrix( BlockIDs::compositionalBlock, BlockIDs::fluidPressureBlock )->Print( std::cout );
    blockSystem->GetMatrix( BlockIDs::compositionalBlock, BlockIDs::compositionalBlock )->Print( std::cout );
    blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock )->Print( std::cout );
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
  return 0;
}

void CompositionalMultiphaseFlow::SolveSystem( EpetraBlockSystem * const blockSystem,
                                               SystemSolverParameters const * const params )
{

}

void
CompositionalMultiphaseFlow::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                                  real64 const scalingFactor,
                                                  DomainPartition * const domain )
{

}

void CompositionalMultiphaseFlow::ResetStateToBeginningOfStep(DomainPartition * const domain)
{

}

void CompositionalMultiphaseFlow::ImplicitStepComplete( real64 const & time,
                                                        real64 const & dt,
                                                        DomainPartition * const domain )
{

}


REGISTER_CATALOG_ENTRY(SolverBase, CompositionalMultiphaseFlow, string const &, ManagedGroup * const)
}// namespace geosx
