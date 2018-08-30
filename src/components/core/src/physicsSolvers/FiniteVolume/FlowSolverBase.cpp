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
 * @file FlowSolverBase.cpp
 */

#include "FlowSolverBase.hpp"

#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace systemSolverInterface;
using namespace multidimensionalArray;

FlowSolverBase::FlowSolverBase( std::string const & name,
                                ManagedGroup * const parent )
  : SolverBase( name, parent ),
    m_gravityFlag(1)
{
  this->RegisterViewWrapper(viewKeyStruct::gravityFlagString,    &m_gravityFlag,        false);
  this->RegisterViewWrapper(viewKeyStruct::discretizationString, &m_discretizationName, false);

  this->RegisterViewWrapper(viewKeyStruct::fluidNameString,  &m_fluidName,  false);
  this->RegisterViewWrapper(viewKeyStruct::fluidIndexString, &m_fluidIndex, false);

  this->RegisterViewWrapper(viewKeyStruct::solidNameString,  &m_solidName,  false);
  this->RegisterViewWrapper(viewKeyStruct::solidIndexString, &m_solidIndex, false);
}

void FlowSolverBase::FillDocumentationNode()
{
  SolverBase::FillDocumentationNode();

  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->AllocateChildNode( viewKeyStruct::gravityFlagString,
                              viewKeyStruct::gravityFlagString,
                              -1,
                              "integer",
                              "integer",
                              "Flag that enables/disables gravity",
                              "Flag that enables/disables gravity",
                              "1",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::discretizationString,
                              viewKeyStruct::discretizationString,
                              -1,
                              "string",
                              "string",
                              "Name of the finite volume discretization to use",
                              "Name of the finite volume discretization to use",
                              "REQUIRED",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::fluidNameString,
                              viewKeyStruct::fluidNameString,
                              -1,
                              "string",
                              "string",
                              "Name of the fluid constitutive model to use",
                              "Name of the fluid constitutive model to use",
                              "REQUIRED",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::solidNameString,
                              viewKeyStruct::solidNameString,
                              -1,
                              "string",
                              "string",
                              "Name of the solid constitutive model to use",
                              "Name of the solid constitutive model to use",
                              "REQUIRED",
                              "",
                              1,
                              1,
                              0 );
}

void FlowSolverBase::FillOtherDocumentationNodes(dataRepository::ManagedGroup * const rootGroup)
{
  SolverBase::FillOtherDocumentationNodes(rootGroup);

  SolverBase::FillOtherDocumentationNodes( rootGroup );
  DomainPartition * domain  = rootGroup->GetGroup<DomainPartition>(keys::domain);

  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel * meshLevel = ManagedGroup::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);
    ElementRegionManager * const elemManager = meshLevel->getElemManager();

    elemManager->forCellBlocks([&](CellBlockSubRegion * const cellBlock) -> void
    {
      cxx_utilities::DocumentationNode * const docNode = cellBlock->getDocumentationNode();

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
                                  3 );
    });

    {
      FaceManager * const faceManager = meshLevel->getFaceManager();
      cxx_utilities::DocumentationNode * const docNode = faceManager->getDocumentationNode();

      docNode->AllocateChildNode( viewKeyStruct::gravityDepthString,
                                  viewKeyStruct::gravityDepthString,
                                  -1,
                                  "real64_array",
                                  "real64_array",
                                  "Precomputed (gravity dot depth)",
                                  "Precomputed (gravity dot depth)",
                                  "",
                                  faceManager->getName(),
                                  1,
                                  0,
                                  1 );
    }
  }
}

void FlowSolverBase::FinalInitialization(dataRepository::ManagedGroup * const rootGroup)
{
  SolverBase::FinalInitialization(rootGroup);

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  // Precompute solver-specific constant data (e.g. gravity-depth)
  PrecomputeData(domain);

  ConstitutiveManager * const cm = domain->getConstitutiveManager();
  m_fluidIndex = cm->GetConstitituveRelation( this->m_fluidName )->getIndexInParent();
  m_solidIndex = cm->GetConstitituveRelation( this->m_solidName )->getIndexInParent();
}

void FlowSolverBase::PrecomputeData(DomainPartition * const domain)
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  FaceManager * const faceManager = mesh->getFaceManager();

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


  r1_array & faceCenter = faceManager->getReference<r1_array>(FaceManager::viewKeyStruct::faceCenterString);
  real64_array & gravityDepthFace = faceManager->getReference<real64_array>(viewKeyStruct::gravityDepthString);

  for (localIndex kf = 0; kf < faceManager->size(); ++kf)
  {
    gravityDepthFace[kf] = Dot(faceCenter[kf], gravityVector);
  }
}

FlowSolverBase::~FlowSolverBase() = default;


} // namespace geosx
