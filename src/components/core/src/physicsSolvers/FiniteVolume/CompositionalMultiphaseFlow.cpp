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

void CompositionalMultiphaseFlow::FillOtherDocumentationNodes(dataRepository::ManagedGroup * const group)
{

}

void CompositionalMultiphaseFlow::FinalInitialization(dataRepository::ManagedGroup * const problemManager)
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
} // namespace geosx
