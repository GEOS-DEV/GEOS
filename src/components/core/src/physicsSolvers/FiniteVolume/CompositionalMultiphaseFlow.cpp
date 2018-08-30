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

  this->RegisterViewWrapper(viewKeyStruct::gravityFlagString, &m_gravityFlag, false);
  this->RegisterViewWrapper(viewKeyStruct::discretizationString, &m_discretizationName, false);
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

                                 docNode->AllocateChildNode(viewKeyStruct::pressureString,
                                                            viewKeyStruct::pressureString,
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

                                 docNode->AllocateChildNode(viewKeyStruct::deltaPressureString,
                                                            viewKeyStruct::deltaPressureString,
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

                                 docNode->AllocateChildNode(viewKeyStruct::globalComponentDensityString,
                                                            viewKeyStruct::globalComponentDensityString,
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

                                 docNode->AllocateChildNode(viewKeyStruct::deltaGlobalComponentDensityString,
                                                            viewKeyStruct::deltaGlobalComponentDensityString,
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

                                 docNode->AllocateChildNode(viewKeyStruct::phaseVolumeFractionString,
                                                            viewKeyStruct::phaseVolumeFractionString,
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

                                 docNode->AllocateChildNode(viewKeyStruct::phaseDensityString,
                                                            viewKeyStruct::phaseDensityString,
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

                                 docNode->AllocateChildNode(viewKeyStruct::phaseComponentDensityString,
                                                            viewKeyStruct::phaseComponentDensityString,
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

                                 docNode->AllocateChildNode(viewKeyStruct::porosityString,
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
                                                            3);

                                 docNode->AllocateChildNode(viewKeyStruct::blockLocalDofNumberString,
                                                            viewKeyStruct::blockLocalDofNumberString,
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

      docNode->AllocateChildNode(viewKeyStruct::facePressureString,
                                 viewKeyStruct::facePressureString,
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

      docNode->AllocateChildNode(viewKeyStruct::phaseDensityString,
                                 viewKeyStruct::phaseDensityString,
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

      docNode->AllocateChildNode(viewKeyStruct::phaseComponentDensityString,
                                 viewKeyStruct::phaseComponentDensityString,
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

      docNode->AllocateChildNode(viewKeyStruct::phaseViscosityString,
                                 viewKeyStruct::phaseViscosityString,
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

      docNode->AllocateChildNode(viewKeyStruct::phaseRelativePermeabilityString,
                                 viewKeyStruct::phaseRelativePermeabilityString,
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

void CompositionalMultiphaseFlow::FinalInitialization(ManagedGroup * const rootGroup)
{
  FlowSolverBase::FinalInitialization( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  ConstitutiveManager * const cm = domain->getConstitutiveManager();
  ConstitutiveBase const * fluid = cm->GetConstitituveRelation( this->m_fluidName );
  // TODO: extract number of phases/components from fluid
  m_numPhases = 2;
  m_numComponents = 3;

  // compute number of DOF per cell
  m_numDofPerCell = m_numComponents + 1;
}

real64 CompositionalMultiphaseFlow::SolverStep(real64 const & time_n,
                                               real64 const & dt,
                                               integer const cycleNumber,
                                               DomainPartition * domain)
{
  // currently the only method is implicit time integration
  return this->NonlinearImplicitStep( time_n,
                                      dt,
                                      cycleNumber,
                                      domain,
                                      getLinearSystemRepository() );
}

void
CompositionalMultiphaseFlow::ImplicitStepSetup(real64 const & time_n, real64 const & dt,
                                               DomainPartition * const domain,
                                               systemSolverInterface::EpetraBlockSystem * const blockSystem)
{

}

void CompositionalMultiphaseFlow::AssembleSystem(DomainPartition * const domain,
                                                 systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                                 real64 const time, real64 const dt)
{

}

void CompositionalMultiphaseFlow::ApplyBoundaryConditions(DomainPartition * const domain,
                                                          systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                                          real64 const time, real64 const dt)
{

}

real64
CompositionalMultiphaseFlow::CalculateResidualNorm(systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                                                   DomainPartition * const domain)
{
  return 0;
}

void CompositionalMultiphaseFlow::SolveSystem(systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                              SystemSolverParameters const * const params)
{

}

void
CompositionalMultiphaseFlow::ApplySystemSolution(systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                                                 real64 const scalingFactor, DomainPartition * const domain)
{

}

void CompositionalMultiphaseFlow::ResetStateToBeginningOfStep(DomainPartition * const domain)
{

}

void CompositionalMultiphaseFlow::ImplicitStepComplete(real64 const & time, real64 const & dt,
                                                       DomainPartition * const domain)
{

}

void CompositionalMultiphaseFlow::SetupSystem(DomainPartition * const domain,
                                              systemSolverInterface::EpetraBlockSystem * const blockSystem)
{

}

void CompositionalMultiphaseFlow::SetSparsityPattern(DomainPartition const * const domain,
                                                     Epetra_FECrsGraph * const sparsity)
{

}

void CompositionalMultiphaseFlow::SetNumRowsAndTrilinosIndices(MeshLevel * const meshLevel, localIndex & numLocalRows,
                                                               globalIndex & numGlobalRows,
                                                               localIndex_array & localIndices, localIndex offset)
{

}

void
CompositionalMultiphaseFlow::ApplyDirichletBC_implicit(DomainPartition * object, real64 const time, real64 const dt,
                                                       systemSolverInterface::EpetraBlockSystem * const blockSystem)
{

}

void
CompositionalMultiphaseFlow::ApplyFaceDirichletBC_implicit(DomainPartition * domain, real64 const time, real64 const dt,
                                                           systemSolverInterface::EpetraBlockSystem * const blockSystem)
{

}


REGISTER_CATALOG_ENTRY(SolverBase, CompositionalMultiphaseFlow, string const &, ManagedGroup * const)
}// namespace geosx
