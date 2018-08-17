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

CompositionalMultiphaseFlow::CompositionalMultiphaseFlow(const std::string & name,
                                                         dataRepository::ManagedGroup * const parent)
  : SolverBase(name, parent),
    m_precomputeDone(false),
    m_gravityFlag(1)
{
// set the blockID for the block system interface
  getLinearSystemRepository()->SetBlockID( BlockIDs::fluidPressureBlock, this->getName() );
  getLinearSystemRepository()->SetBlockID( BlockIDs::compositionalBlock, this->getName() );

  this->RegisterViewWrapper(viewKeyStruct::gravityFlagString, &m_gravityFlag, false);
  this->RegisterViewWrapper(viewKeyStruct::discretizationString, &m_discretizationName, false);
}

void CompositionalMultiphaseFlow::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  SolverBase::FillDocumentationNode();

  docNode->setName(this->CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("An example single phase flow solver");

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

void CompositionalMultiphaseFlow::FillOtherDocumentationNodes(dataRepository::ManagedGroup * const rootGroup)
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
                                  1,
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
                                  1,
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
                                  1,
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
                                  1,
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
                                  1,
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
                                  1,
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
                                  1,
                                  0,
                                  0 );

    });
}

void CompositionalMultiphaseFlow::FinalInitialization(dataRepository::ManagedGroup * const rootGroup)
{

}

real64 CompositionalMultiphaseFlow::SolverStep(real64 const & time_n, real64 const & dt, integer const cycleNumber,
                                               DomainPartition * domain)
{
  return dt;
}

void
CompositionalMultiphaseFlow::ImplicitStepSetup(real64 const & time_n, real64 const & dt, DomainPartition * const domain,
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

REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphaseFlow, std::string const &, ManagedGroup * const )
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

void CompositionalMultiphaseFlow::PrecomputeData(DomainPartition * const domain)
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

void CompositionalMultiphaseFlow::AllocateAuxStorage(DomainPartition * const domain)
{

}
// namespace geosx
