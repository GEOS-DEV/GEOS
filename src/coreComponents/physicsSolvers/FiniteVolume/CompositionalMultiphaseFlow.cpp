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
}

void CompositionalMultiphaseFlow::FillDocumentationNode()
{
  SolverBase::FillDocumentationNode();

  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName(CompositionalMultiphaseFlow::CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("A compositional multiphase flow solver");
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

                                 docNode->AllocateChildNode(viewKeys.globalComponentDensity.Key(),
                                                            viewKeys.globalComponentDensity.Key(),
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

                                 docNode->AllocateChildNode(viewKeys.deltaGlobalComponentDensity.Key(),
                                                            viewKeys.deltaGlobalComponentDensity.Key(),
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

                                 docNode->AllocateChildNode(viewKeys.phaseComponentDensity.Key(),
                                                            viewKeys.phaseComponentDensity.Key(),
                                                            -1,
                                                            "real64_array3d",
                                                            "real64_array3d",
                                                            "Fluid component-in-phase density",
                                                            "Fluid component-in-phase density",
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

      docNode->AllocateChildNode(viewKeys.phaseComponentDensity.Key(),
                                 viewKeys.phaseComponentDensity.Key(),
                                 -1,
                                 "real64_array3d",
                                 "real64_array3d",
                                 "Fluid component-in-phase density",
                                 "Fluid component-in-phase density",
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

  // TODO: this should be available without casting to a specific type
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
      cellBlock->getReference<array2d<real64>>(viewKeys.globalComponentDensity).resizeDimension<1>(m_numComponents);
      cellBlock->getReference<array2d<real64>>(viewKeys.deltaGlobalComponentDensity).resizeDimension<1>(m_numComponents);
      cellBlock->getReference<array2d<real64>>(viewKeys.phaseVolumeFraction).resizeDimension<1>(m_numPhases);
      cellBlock->getReference<array2d<real64>>(viewKeys.phaseDensity).resizeDimension<1>(m_numPhases);

      array3d<real64> & phaseCompDens = cellBlock->getReference<array3d<real64>>(viewKeys.phaseComponentDensity);
      phaseCompDens.resize(phaseCompDens.size(0), m_numPhases, m_numComponents);

      cellBlock->getReference<array2d<real64>>(viewKeys.blockLocalDofNumber).resizeDimension<1>(m_numDofPerCell);
    });

    {
      FaceManager * const faceManager = meshLevel->getFaceManager();

      faceManager->getReference<array2d<real64>>(viewKeys.phaseDensity).resizeDimension<1>(m_numPhases);
      faceManager->getReference<array2d<real64>>(viewKeys.phaseViscosity).resizeDimension<1>(m_numPhases);
      faceManager->getReference<array2d<real64>>(viewKeys.phaseRelativePermeability).resizeDimension<1>(m_numPhases);

      array3d<real64> & phaseCompDens = faceManager->getReference<array3d<real64>>(viewKeys.phaseComponentDensity);
      phaseCompDens.resize(phaseCompDens.size(0), m_numPhases, m_numComponents);
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
  fieldNames["elems"].push_back(viewKeys.globalComponentDensity.Key());
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

void
CompositionalMultiphaseFlow::ImplicitStepSetup( real64 const & time_n, real64 const & dt,
                                                DomainPartition * const domain,
                                                EpetraBlockSystem * const blockSystem )
{

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
