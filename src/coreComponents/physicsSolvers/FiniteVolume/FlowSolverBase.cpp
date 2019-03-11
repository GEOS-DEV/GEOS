/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

FlowSolverBase::FlowSolverBase( std::string const & name,
                                ManagedGroup * const parent )
  : SolverBase( name, parent ),
    m_gravityFlag(1),
    m_fluidName(),
    m_solidName(),
    m_fluidIndex(),
    m_solidIndex(),
    m_poroElasticFlag(0),
    m_reservoirWellsSystemFlag(0),
    m_numDofPerCell(0),
    m_elemGhostRank(),
    m_volume(),
    m_gravDepth(),
    m_porosityRef()
{
  RegisterViewWrapper( viewKeyStruct::gravityFlagString, &m_gravityFlag, false )->
    setApplyDefaultValue(1)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Flag that enables/disables gravity");
  
  this->RegisterViewWrapper( viewKeyStruct::discretizationString, &m_discretizationName, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of discretization object to use for this solver.");

  this->RegisterViewWrapper( viewKeyStruct::fluidNameString,  &m_fluidName,  false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of fluid constitutive object to use for this solver.");

  this->RegisterViewWrapper( viewKeyStruct::solidNameString,  &m_solidName,  false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of solid constitutive object to use for this solver");

  this->RegisterViewWrapper( viewKeyStruct::fluidIndexString, &m_fluidIndex, false );
  this->RegisterViewWrapper( viewKeyStruct::solidIndexString, &m_solidIndex, false );
  
}

void FlowSolverBase::RegisterDataOnMesh( ManagedGroup * const MeshBodies )
{
  SolverBase::RegisterDataOnMesh( MeshBodies );

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * meshLevel = ManagedGroup::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);

    ElementRegionManager * const elemManager = meshLevel->getElemManager();
    elemManager->forElementSubRegions([&]( auto * const elementSubRegion) -> void
    {
      elementSubRegion->template RegisterViewWrapper< array1d<real64> >( viewKeyStruct::referencePorosityString )->setPlotLevel(PlotLevel::LEVEL_0);
      elementSubRegion->template RegisterViewWrapper< array1d<R1Tensor> >( viewKeyStruct::permeabilityString )->setPlotLevel(PlotLevel::LEVEL_0);
      elementSubRegion->template RegisterViewWrapper< array1d<real64> >( viewKeyStruct::gravityDepthString )->setApplyDefaultValue( 0.0 );
    });

    FaceManager * const faceManager = meshLevel->getFaceManager();
    faceManager->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::gravityDepthString )->setApplyDefaultValue( 0.0 );
  }
}

void FlowSolverBase::SetSparsityPattern( DomainPartition const * const domain,
                                         Epetra_FECrsGraph * const sparsity )
{
}

void FlowSolverBase::SetNumRowsAndTrilinosIndices( MeshLevel * const meshLevel,
                                                   localIndex & numLocalRows,
                                                   globalIndex & numGlobalRows,
                                                   localIndex offset )
{
}

void FlowSolverBase::InitializePreSubGroups(ManagedGroup * const rootGroup)
{
  SolverBase::InitializePreSubGroups(rootGroup);

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  ConstitutiveManager * const cm = domain->getConstitutiveManager();

  ConstitutiveBase const * fluid  = cm->GetConstitituveRelation<ConstitutiveBase>( m_fluidName );
  GEOS_ERROR_IF( fluid == nullptr, "Fluid model " + m_fluidName + " not found" );
  m_fluidIndex = fluid->getIndexInParent();

  ConstitutiveBase const * solid  = cm->GetConstitituveRelation<ConstitutiveBase>( m_solidName );
  GEOS_ERROR_IF( solid == nullptr, "Solid model " + m_solidName + " not found" );
  m_solidIndex = solid->getIndexInParent();
}

void FlowSolverBase::InitializePostInitialConditions_PreSubGroups(ManagedGroup * const rootGroup)
{
  SolverBase::InitializePostInitialConditions_PreSubGroups(rootGroup);

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  ResetViews( domain );

  // Precompute solver-specific constant data (e.g. gravity-depth)
  PrecomputeData(domain);
}

void FlowSolverBase::PrecomputeData(DomainPartition * const domain)
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  FaceManager * const faceManager = mesh->getFaceManager();

  R1Tensor const & gravityVector = getGravityVector();

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> const elemCenter =
    elemManager->ConstructViewAccessor<array1d<R1Tensor>, arrayView1d<R1Tensor>>( CellBlock::viewKeyStruct::elementCenterString );

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> gravityDepth =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>(viewKeyStruct::gravityDepthString);

  // Loop over all the elements and calculate element centers, and element volumes
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k )->void
  {
    gravityDepth[er][esr][k] = Dot(elemCenter[er][esr][k], gravityVector);
  });


  r1_array & faceCenter = faceManager->getReference<r1_array>(FaceManager::viewKeyStruct::faceCenterString);
  real64_array & gravityDepthFace = faceManager->getReference<real64_array>(viewKeysFlowSolverBase.gravityDepth);

  for (localIndex kf = 0; kf < faceManager->size(); ++kf)
  {
    gravityDepthFace[kf] = Dot(faceCenter[kf], gravityVector);
  }
}

FlowSolverBase::~FlowSolverBase() = default;

void FlowSolverBase::ResetViews( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();

  m_elemGhostRank =
    elemManager->ConstructViewAccessor<array1d<integer>, arrayView1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );
  m_volume =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( CellElementSubRegion::viewKeyStruct::elementVolumeString );
  m_gravDepth =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::gravityDepthString );
  m_porosityRef =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::referencePorosityString );
}


} // namespace geosx
