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

#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/Fluid/MultiFluidBase.hpp"
#include "constitutive/RelPerm/RelativePermeabilityBase.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/BoundaryConditions/BoundaryConditionManager.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "systemSolverInterface/LinearSolverWrapper.hpp"
#include "systemSolverInterface/EpetraBlockSystem.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace systemSolverInterface;

CompositionalMultiphaseFlow::CompositionalMultiphaseFlow( const string & name,
                                                          ManagedGroup * const parent )
  :
  FlowSolverBase( name, parent ),
  m_numPhases( 0 ),
  m_numComponents( 0 )
{
  // set the blockID for the block system interface
  getLinearSystemRepository()->SetBlockID(BlockIDs::compositionalBlock, this->getName());

  this->RegisterViewWrapper( viewKeyStruct::temperatureString, &m_temperature, false );
  this->RegisterViewWrapper( viewKeyStruct::useMassFlagString, &m_useMass, false );

  this->RegisterViewWrapper( viewKeyStruct::relPermNameString,  &m_relPermName,  false );
  this->RegisterViewWrapper( viewKeyStruct::relPermIndexString, &m_relPermIndex, false );
}

localIndex CompositionalMultiphaseFlow::numFluidComponents() const
{
  return m_numComponents;
}

localIndex CompositionalMultiphaseFlow::numFluidPhases() const
{
  return m_numPhases;
}

void CompositionalMultiphaseFlow::FillDocumentationNode()
{
  FlowSolverBase::FillDocumentationNode();

  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName(CompositionalMultiphaseFlow::CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("A compositional multiphase flow solver");

  docNode->AllocateChildNode( viewKeyStruct::temperatureString,
                              viewKeyStruct::temperatureString,
                              -1,
                              "real64",
                              "real64",
                              "Temperature",
                              "Temperature",
                              "REQUIRED",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::useMassFlagString,
                              viewKeyStruct::useMassFlagString,
                              -1,
                              "integer",
                              "integer",
                              "Use mass formulation instead of molar",
                              "Use mass formulation instead of molar",
                              "0",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::relPermNameString,
                              viewKeyStruct::relPermNameString,
                              -1,
                              "string",
                              "string",
                              "Name of the relative permeability constitutive model to use",
                              "Name of the relative permeability constitutive model to use",
                              "REQUIRED",
                              "",
                              1,
                              1,
                              0 );
}

void CompositionalMultiphaseFlow::FillOtherDocumentationNodes( ManagedGroup * const rootGroup )
{
  FlowSolverBase::FillOtherDocumentationNodes(rootGroup);

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  for (auto & mesh : domain->getMeshBodies()->GetSubGroups())
  {
    MeshLevel * meshLevel = ManagedGroup::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);
    ElementRegionManager * const elemManager = meshLevel->getElemManager();

    elemManager->forCellBlocks([&](CellBlockSubRegion * const cellBlock)
    {
      cxx_utilities::DocumentationNode * const docNode = cellBlock->getDocumentationNode();

      docNode->AllocateChildNode( viewKeyStruct::pressureString,
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
                                  0 );

      docNode->AllocateChildNode( viewKeyStruct::deltaPressureString,
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
                                  1 );

      // TODO this needs to be allocated on BC sets only
      docNode->AllocateChildNode( viewKeyStruct::bcPressureString,
                                  viewKeyStruct::bcPressureString,
                                  -1,
                                  "real64_array",
                                  "real64_array",
                                  "Boundary condition pressure",
                                  "Boundary condition pressure",
                                  "",
                                  elemManager->getName(),
                                  1,
                                  0,
                                  1 );

      docNode->AllocateChildNode( viewKeyStruct::globalCompDensityString,
                                  viewKeyStruct::globalCompDensityString,
                                  -1,
                                  "real64_array2d",
                                  "real64_array2d",
                                  "Global component density in mixture",
                                  "Global component density in mixture",
                                  "",
                                  elemManager->getName(),
                                  1,
                                  0,
                                  0 );

      docNode->AllocateChildNode( viewKeyStruct::deltaGlobalCompDensityString,
                                  viewKeyStruct::deltaGlobalCompDensityString,
                                  -1,
                                  "real64_array2d",
                                  "real64_array2d",
                                  "Change in global component density in mixture",
                                  "Change in global component density in mixture",
                                  "",
                                  elemManager->getName(),
                                  1,
                                  0,
                                  1 );

      docNode->AllocateChildNode( viewKeyStruct::globalCompFractionString,
                                  viewKeyStruct::globalCompFractionString,
                                  -1,
                                  "real64_array2d",
                                  "real64_array2d",
                                  "Global component mole fraction in mixture",
                                  "Global component mole fraction in mixture",
                                  "",
                                  elemManager->getName(),
                                  1,
                                  0,
                                  0 );

      docNode->AllocateChildNode( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString,
                                  viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString,
                                  -1,
                                  "real64_array3d",
                                  "real64_array3d",
                                  "Derivatives of global component mole fraction",
                                  "Derivatives of global component mole fraction",
                                  "",
                                  elemManager->getName(),
                                  1,
                                  0,
                                  3 );

      docNode->AllocateChildNode( viewKeyStruct::phaseVolumeFractionString,
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
                                  0 );

      docNode->AllocateChildNode( viewKeyStruct::dPhaseVolumeFraction_dPressureString,
                                  viewKeyStruct::dPhaseVolumeFraction_dPressureString,
                                  -1,
                                  "real64_array2d",
                                  "real64_array2d",
                                  "Fluid phase volume fraction derivative w.r.t. pressure",
                                  "Fluid phase volume fraction derivative w.r.t. pressure",
                                  "",
                                  elemManager->getName(),
                                  1,
                                  0,
                                  3 );

      docNode->AllocateChildNode( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString,
                                  viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString,
                                  -1,
                                  "real64_array3d",
                                  "real64_array3d",
                                  "Fluid phase volume fraction derivative w.r.t. component density",
                                  "Fluid phase volume fraction derivative w.r.t. component density",
                                  "",
                                  elemManager->getName(),
                                  1,
                                  0,
                                  3 );

      docNode->AllocateChildNode( viewKeyStruct::phaseVolumeFractionOldString,
                                  viewKeyStruct::phaseVolumeFractionOldString,
                                  -1,
                                  "real64_array2d",
                                  "real64_array2d",
                                  "Fluid phase volume fraction (old)",
                                  "Fluid phase volume fraction (old)",
                                  "",
                                  elemManager->getName(),
                                  1,
                                  0,
                                  3 );

      docNode->AllocateChildNode( viewKeyStruct::phaseDensityOldString,
                                  viewKeyStruct::phaseDensityOldString,
                                  -1,
                                  "real64_array2d",
                                  "real64_array2d",
                                  "Fluid phase density (old)",
                                  "Fluid phase density (old)",
                                  "",
                                  elemManager->getName(),
                                  1,
                                  0,
                                  3 );

      docNode->AllocateChildNode( viewKeyStruct::phaseComponentFractionOldString,
                                  viewKeyStruct::phaseComponentFractionOldString,
                                  -1,
                                  "real64_array3d",
                                  "real64_array3d",
                                  "Fluid component-in-phase fraction (old)",
                                  "Fluid component-in-phase fraction (old)",
                                  "",
                                  elemManager->getName(),
                                  1,
                                  0,
                                  3 );

      docNode->AllocateChildNode( viewKeyStruct::porosityOldString,
                                  viewKeyStruct::porosityOldString,
                                  -1,
                                  "real64_array",
                                  "real64_array",
                                  "Porosity (old)",
                                  "Porosity (old)",
                                  "",
                                  elemManager->getName(),
                                  1,
                                  0,
                                  3 );

      docNode->AllocateChildNode( viewKeyStruct::blockLocalDofNumberString,
                                  viewKeyStruct::blockLocalDofNumberString,
                                  -1,
                                  "globalIndex_array",
                                  "globalIndex_array",
                                  "Pressure DOF index",
                                  "Pressure DOF index",
                                  "0",
                                  "",
                                  1,
                                  0,
                                  3 );
    });

    {
      FaceManager * const faceManager = meshLevel->getFaceManager();
      cxx_utilities::DocumentationNode * const docNode = faceManager->getDocumentationNode();

      docNode->AllocateChildNode( viewKeyStruct::facePressureString,
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
                                  3 );

      docNode->AllocateChildNode( viewKeyStruct::globalCompFractionString,
                                  viewKeyStruct::globalCompFractionString,
                                  -1,
                                  "real64_array2d",
                                  "real64_array2d",
                                  "Global component mole fraction in mixture",
                                  "Global component mole fraction in mixture",
                                  "",
                                  elemManager->getName(),
                                  1,
                                  0,
                                  3 );

      docNode->AllocateChildNode( viewKeyStruct::phaseDensityOldString,
                                  viewKeyStruct::phaseDensityOldString,
                                  -1,
                                  "real64_array2d",
                                  "real64_array2d",
                                  "Fluid phase density",
                                  "Fluid phase density",
                                  "",
                                  faceManager->getName(),
                                  1,
                                  0,
                                  3 );

      docNode->AllocateChildNode( viewKeyStruct::phaseComponentFractionOldString,
                                  viewKeyStruct::phaseComponentFractionOldString,
                                  -1,
                                  "real64_array3d",
                                  "real64_array3d",
                                  "Fluid component-in-phase fraction",
                                  "Fluid component-in-phase fraction",
                                  "",
                                  faceManager->getName(),
                                  1,
                                  0,
                                  3 );

      docNode->AllocateChildNode( viewKeyStruct::phaseViscosityString,
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
                                  3 );

      docNode->AllocateChildNode( viewKeyStruct::phaseRelativePermeabilityString,
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
                                  3 );
    }
  }
}

void CompositionalMultiphaseFlow::InitializePreSubGroups( ManagedGroup * const rootGroup )
{
  FlowSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition     const * domain = rootGroup->GetGroup<DomainPartition>( keys::domain );
  ConstitutiveManager const * cm     = domain->getConstitutiveManager();

  MultiFluidBase const * fluid = cm->GetConstitituveRelation<MultiFluidBase>( m_fluidName );
  m_numPhases     = fluid->numFluidPhases();
  m_numComponents = fluid->numFluidComponents();
  m_numDofPerCell = m_numComponents + 1;

  RelativePermeabilityBase const * relPerm = cm->GetConstitituveRelation<RelativePermeabilityBase>( m_relPermName );
  GEOS_ERROR_IF( relPerm == nullptr, "Relative permeability model " + m_relPermName + " not found" );
  m_relPermIndex = relPerm->getIndexInParent();

  // Consistency check between the models
  GEOS_ERROR_IF( fluid->numFluidPhases() != relPerm->numFluidPhases(),
                 "Number of fluid phases differs between fluid model '" << m_fluidName
                 << "' and relperm model '" << m_relPermName << "'" );

  for (localIndex ip = 0; ip < m_numPhases; ++ip)
  {
    string const & phase_fl = fluid->phaseName( ip );
    string const & phase_rp = relPerm->phaseName( ip );
    GEOS_ERROR_IF( phase_fl != phase_rp, "Phase '" << phase_fl << "' in fluid model '" << m_fluidName
                   << "' does not match phase '" << phase_rp << "' in relperm model '" << m_relPermName << "'" );
  }
}

void CompositionalMultiphaseFlow::ResizeFields( DomainPartition * domain )
{
  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  for (auto & mesh : domain->getMeshBodies()->GetSubGroups())
  {
    MeshLevel * meshLevel = ManagedGroup::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);
    ElementRegionManager * const elemManager = meshLevel->getElemManager();

    elemManager->forCellBlocks( [&] ( CellBlockSubRegion * const subRegion )
    {
      subRegion->getReference<array2d<real64>>(viewKeyStruct::globalCompDensityString).resizeDimension<1>(NC);
      subRegion->getReference<array2d<real64>>(viewKeyStruct::deltaGlobalCompDensityString).resizeDimension<1>(NC);

      subRegion->getReference<array2d<real64>>(viewKeyStruct::globalCompFractionString).resizeDimension<1>(NC);
      subRegion->getReference<array3d<real64>>(viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString).resizeDimension<1,2>(NC, NC);

      subRegion->getReference<array2d<real64>>(viewKeyStruct::phaseVolumeFractionString).resizeDimension<1>(NP);
      subRegion->getReference<array2d<real64>>(viewKeyStruct::dPhaseVolumeFraction_dPressureString).resizeDimension<1>(NP);
      subRegion->getReference<array3d<real64>>(viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString).resizeDimension<1,2>(NP, NC);

      subRegion->getReference<array2d<real64>>(viewKeyStruct::phaseVolumeFractionOldString).resizeDimension<1>(NP);
      subRegion->getReference<array2d<real64>>(viewKeyStruct::phaseDensityOldString).resizeDimension<1>(NP);
      subRegion->getReference<array3d<real64>>(viewKeyStruct::phaseComponentFractionOldString).resizeDimension<1,2>(NP, NC);
    });

    {
      FaceManager * const faceManager = meshLevel->getFaceManager();

      faceManager->getReference<array2d<real64>>(viewKeyStruct::phaseDensityOldString).resizeDimension<1>(NP);
      faceManager->getReference<array2d<real64>>(viewKeyStruct::phaseViscosityString).resizeDimension<1>(NP);
      faceManager->getReference<array2d<real64>>(viewKeyStruct::phaseRelativePermeabilityString).resizeDimension<1>(NP);
      faceManager->getReference<array3d<real64>>(viewKeyStruct::phaseComponentFractionOldString).resizeDimension<1,2>(NP, NC);;
    }
  }
}

void CompositionalMultiphaseFlow::IntermediateInitializationPreSubGroups( ManagedGroup * const rootGroup )
{
  FlowSolverBase::IntermediateInitializationPreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  ResizeFields( domain );
}

MultiFluidBase * CompositionalMultiphaseFlow::GetFluidModel( ManagedGroup * dataGroup ) const
{
  ManagedGroup * const constitutiveModels = dataGroup->GetGroup( ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  MultiFluidBase * const fluid = constitutiveModels->GetGroup<MultiFluidBase>( m_fluidName );
  GEOS_ERROR_IF( fluid == nullptr, "Target group does not contain the fluid model" );

  return fluid;
}

MultiFluidBase const * CompositionalMultiphaseFlow::GetFluidModel( ManagedGroup const * dataGroup ) const
{
  ManagedGroup const * const constitutiveModels = dataGroup->GetGroup( ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  MultiFluidBase const * const fluid = constitutiveModels->GetGroup<MultiFluidBase>( m_fluidName );
  GEOS_ERROR_IF( fluid == nullptr, "Target group does not contain the fluid model" );

  return fluid;
}

ConstitutiveBase * CompositionalMultiphaseFlow::GetSolidModel( ManagedGroup * dataGroup ) const
{
  ManagedGroup * const contitutiveModels = dataGroup->GetGroup( ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( contitutiveModels == nullptr, "Target group does not contain constitutive models" );

  ConstitutiveBase * const solid = contitutiveModels->GetGroup<ConstitutiveBase>( m_solidName );
  GEOS_ERROR_IF( solid == nullptr, "Target group does not contain the solid model" );

  return solid;
}

ConstitutiveBase const * CompositionalMultiphaseFlow::GetSolidModel( ManagedGroup const * dataGroup ) const
{
  ManagedGroup const * const contitutiveModels = dataGroup->GetGroup( ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( contitutiveModels == nullptr, "Target group does not contain constitutive models" );

  ConstitutiveBase const * const solid = contitutiveModels->GetGroup<ConstitutiveBase>( m_solidName );
  GEOS_ERROR_IF( solid == nullptr, "Target group does not contain the solid model" );

  return solid;
}

RelativePermeabilityBase * CompositionalMultiphaseFlow::GetRelPermModel( ManagedGroup * dataGroup ) const
{
  ManagedGroup * const constitutiveModels = dataGroup->GetGroup( ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  RelativePermeabilityBase * const relPerm = constitutiveModels->GetGroup<RelativePermeabilityBase>( m_relPermName );
  GEOS_ERROR_IF( relPerm == nullptr, "Target group does not contain the relative permeability model" );

  return relPerm;
}

RelativePermeabilityBase const * CompositionalMultiphaseFlow::GetRelPermModel( ManagedGroup const * dataGroup ) const
{
  ManagedGroup const * const constitutiveModels = dataGroup->GetGroup( ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  RelativePermeabilityBase const * const relPerm = constitutiveModels->GetGroup<RelativePermeabilityBase>( m_relPermName );
  GEOS_ERROR_IF( relPerm == nullptr, "Target group does not contain the relative permeability model" );

  return relPerm;
}

namespace
{

template<typename LAMBDA>
typename std::enable_if<std::is_convertible<LAMBDA,
                                            std::function<void(localIndex,localIndex,ElementRegion *,CellBlockSubRegion *)>
                                           >::value, void>::type
applyToSubRegions( DomainPartition * const domain, LAMBDA && lambda )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  elemManager->forCellBlocksComplete( lambda );
}

template<typename LAMBDA>
typename std::enable_if<std::is_convertible<LAMBDA,std::function<void(CellBlockSubRegion *)>>::value, void>::type
applyToSubRegions( DomainPartition * const domain, LAMBDA && lambda )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  elemManager->forCellBlocks( lambda );
}

template<typename LAMBDA>
typename std::enable_if<std::is_convertible<LAMBDA,
                                            std::function<void(localIndex,localIndex,ElementRegion const *,CellBlockSubRegion const *)>
                                           >::value, void>::type
applyToSubRegions( DomainPartition const * const domain, LAMBDA && lambda )
{
  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  elemManager->forCellBlocksComplete( lambda );
}

}

void CompositionalMultiphaseFlow::UpdateComponentFraction( ManagedGroup * const dataGroup )
{
  arrayView2d<real64 const> const & compDens =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::globalCompDensityString );

  arrayView2d<real64 const> const & dCompDens =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

  arrayView2d<real64> const & compFrac =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::globalCompFractionString );

  arrayView3d<real64> const & dCompFrac_dCompDens =
    dataGroup->getReference<array3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  forall_in_range( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    real64 totalDensity = 0.0;

    for (localIndex ic = 0; ic < m_numComponents; ++ic)
    {
      totalDensity += compDens[a][ic] + dCompDens[a][ic];
    }

    real64 const totalDensityInv = 1.0 / totalDensity;

    for (localIndex ic = 0; ic < m_numComponents; ++ic)
    {
      compFrac[a][ic] = (compDens[a][ic] + dCompDens[a][ic]) * totalDensityInv;
      for (localIndex jc = 0; jc < m_numComponents; ++jc)
      {
        dCompFrac_dCompDens[a][ic][jc] = - compFrac[a][ic] * totalDensityInv;
      }
      dCompFrac_dCompDens[a][ic][ic] += totalDensityInv;
    }
  });
}

void CompositionalMultiphaseFlow::UpdateComponentFractionAll( DomainPartition * const domain )
{
  applyToSubRegions( domain, [&] ( CellBlockSubRegion * subRegion )
  {
    UpdateComponentFraction( subRegion );
  });
}

void CompositionalMultiphaseFlow::UpdatePhaseVolumeFraction( ManagedGroup * const dataGroup )
{
  MultiFluidBase * fluid = GetFluidModel( dataGroup );

  arrayView2d<real64> const & phaseVolFrac =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::phaseVolumeFractionString );

  arrayView2d<real64> const & dPhaseVolFrac_dPres =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dPressureString );

  arrayView3d<real64> const & dPhaseVolFrac_dComp =
    dataGroup->getReference<array3d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  arrayView3d<real64> const & dCompFrac_dCompDens =
    dataGroup->getReference<array3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  arrayView2d<real64 const> const & compDens =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::globalCompDensityString );

  arrayView2d<real64 const> const & dCompDens =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

  arrayView3d<real64 const> const & phaseFrac =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseFractionString );

  arrayView3d<real64 const> const & dPhaseFrac_dPres =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dPressureString );

  arrayView4d<real64 const> const & dPhaseFrac_dComp =
    fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dGlobalCompFractionString );

  arrayView3d<real64 const> const & phaseDens =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseDensityString );

  arrayView3d<real64 const> const & dPhaseDens_dPres =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString );

  arrayView4d<real64 const> const & dPhaseDens_dComp =
    fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dGlobalCompFractionString );

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  forall_in_range( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    stackArray1d<real64, maxNumComp> work( NC );

    // compute total density from component partial densities
    real64 totalDensity = 0.0;
    real64 const dTotalDens_dCompDens = 1.0;
    for (localIndex ic = 0; ic < NC; ++ic)
    {
      totalDensity += compDens[a][ic] + dCompDens[a][ic];
    }

    for (localIndex ip = 0; ip < NP; ++ip)
    {
      // Expression for volume fractions: S_p = (nu_p / rho_p) * rho_t
      real64 const phaseDensInv = 1.0 / phaseDens[a][0][ip];

      // compute saturation and derivatives except multiplying by the total density
      phaseVolFrac[a][ip] = phaseFrac[a][0][ip] * phaseDensInv;

      dPhaseVolFrac_dPres[a][ip] =
        (dPhaseFrac_dPres[a][0][ip] - phaseVolFrac[a][ip] * dPhaseDens_dPres[a][0][ip]) * phaseDensInv;

      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dPhaseVolFrac_dComp[a][ip][jc] =
          (dPhaseFrac_dComp[a][0][ip][jc] - phaseVolFrac[a][ip] * dPhaseDens_dComp[a][0][ip][jc]) * phaseDensInv;
      }

      // apply chain rule to convert derivatives from global component fractions to densities
      applyChainRuleInPlace( NC, dCompFrac_dCompDens[a], dPhaseVolFrac_dComp[a][ip], work );

      // now finalize the computation by multiplying by total density
      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dPhaseVolFrac_dComp[a][ip][jc] *= totalDensity;
        dPhaseVolFrac_dComp[a][ip][jc] += phaseVolFrac[a][ip] * dTotalDens_dCompDens;
      }

      phaseVolFrac[a][ip] *= totalDensity;
      dPhaseVolFrac_dPres[a][ip] *= totalDensity;
    }
  });
}

void CompositionalMultiphaseFlow::UpdatePhaseVolumeFractionAll( DomainPartition * const domain )
{
  applyToSubRegions( domain, [&] ( CellBlockSubRegion * subRegion )
  {
    UpdatePhaseVolumeFraction( subRegion );
  });
}

void CompositionalMultiphaseFlow::UpdateFluidModel( ManagedGroup * const dataGroup )
{
  MultiFluidBase * const fluid = GetFluidModel( dataGroup );

  arrayView1d<real64 const> const & pres = dataGroup->getReference<array1d<real64>>( viewKeyStruct::pressureString );
  arrayView1d<real64 const> const & dPres = dataGroup->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );
  arrayView2d<real64 const> const & compFrac = dataGroup->getReference<array2d<real64>>( viewKeyStruct::globalCompFractionString );

  // currently not thread-safe, force sequential policy
  forall_in_range<RAJA::seq_exec>( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    fluid->StateUpdatePointMultiFluid(pres[a] + dPres[a], m_temperature, compFrac[a], a, 0);
  });
}

void CompositionalMultiphaseFlow::UpdateFluidModelAll( DomainPartition * const domain )
{
  applyToSubRegions( domain, [&] ( CellBlockSubRegion * subRegion )
  {
    UpdateFluidModel( subRegion );
  });
}

void CompositionalMultiphaseFlow::UpdateSolidModel( ManagedGroup * dataGroup )
{
  ConstitutiveBase * const solid = GetSolidModel( dataGroup );

  arrayView1d<real64 const> const & pres  = dataGroup->getReference<array1d<real64>>( viewKeyStruct::pressureString );
  arrayView1d<real64 const> const & dPres = dataGroup->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

  forall_in_range( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    solid->StateUpdatePointPressure( pres[a] + dPres[a], a, 0 );
  });
}

void CompositionalMultiphaseFlow::UpdateSolidModelAll( DomainPartition * const domain )
{
  applyToSubRegions( domain, [&] ( CellBlockSubRegion * subRegion )
  {
    UpdateSolidModel( subRegion );
  });
}

void CompositionalMultiphaseFlow::UpdateRelPermModel( ManagedGroup * dataGroup )
{
  RelativePermeabilityBase * const relPerm = GetRelPermModel( dataGroup );

  arrayView2d<real64> const & phaseVolFrac =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::phaseVolumeFractionString );

  forall_in_range( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    relPerm->StateUpdatePointRelPerm( phaseVolFrac[a], a, 0 );
  });
}

void CompositionalMultiphaseFlow::UpdateRelPermModelAll( DomainPartition * const domain )
{
  applyToSubRegions( domain, [&] ( CellBlockSubRegion * subRegion )
  {
    UpdateRelPermModel( subRegion );
  });
}

void CompositionalMultiphaseFlow::UpdateState( ManagedGroup * const dataGroup )
{
  UpdateComponentFraction( dataGroup );
  UpdateFluidModel( dataGroup );
  UpdatePhaseVolumeFraction( dataGroup );
  UpdateSolidModel( dataGroup );
  UpdateRelPermModel( dataGroup );
}

void CompositionalMultiphaseFlow::UpdateStateAll( DomainPartition * const domain )
{
  UpdateComponentFractionAll( domain );
  UpdateFluidModelAll( domain );
  UpdatePhaseVolumeFractionAll( domain );
  UpdateSolidModelAll( domain );
  UpdateRelPermModelAll( domain );
}

void CompositionalMultiphaseFlow::InitializeFluidState( DomainPartition * const domain )
{
  // 1. Assume global component fractions have been prescribed.
  // Initialize constitutive state to get fluid density.
  UpdateFluidModelAll( domain );

  // 2. Back-calculate global component densities from fractions and total fluid density
  // in order to initialize the primary solution variables
  applyToSubRegions( domain, [&] ( localIndex er, localIndex esr,
                                   ElementRegion * const region,
                                   CellBlockSubRegion * const subRegion )
  {
    arrayView2d<real64 const> const & compFrac  = m_compFrac[er][esr];
    arrayView2d<real64 const> const & totalDens = m_totalDens[er][esr][m_fluidIndex];

    arrayView2d<real64> compDens = m_globalCompDensity[er][esr];

    for_elems_in_subRegion( subRegion, GEOSX_LAMBDA ( localIndex ei )
    {
      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        compDens[ei][ic] = totalDens[ei][0] * compFrac[ei][ic];
      }
    });
  });

  // 3. Calculate phase saturations
  UpdatePhaseVolumeFractionAll(domain);

  // 4. Initialize solid state
  UpdateSolidModelAll( domain );

  // 5. Initialize rel perm state
  UpdateRelPermModelAll( domain );
}

void CompositionalMultiphaseFlow::FinalInitializationPreSubGroups( ManagedGroup * const rootGroup )
{
  FlowSolverBase::FinalInitializationPreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  //TODO this is a hack until the sets are fixed to include ghosts!!
  std::map<string, string_array> fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::pressureString );
  fieldNames["elems"].push_back( viewKeyStruct::globalCompDensityString );
  CommunicationTools::SynchronizeFields(fieldNames,
                                        domain->getMeshBody(0)->getMeshLevel(0),
                                        domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  ConstitutiveManager * constitutiveManager = domain->getConstitutiveManager();

  // set mass fraction flag on main model
  // TODO find a way to set this before constitutive model is duplicated and attached to subregions?
  {
    MultiFluidBase * fluid = constitutiveManager->GetConstitituveRelation<MultiFluidBase>( m_fluidName );
    fluid->setMassFlag( static_cast<bool>(m_useMass) );
  }

  // set mass fraction flag on subregion models
  applyToSubRegions( domain, [&] ( CellBlockSubRegion * subRegion )
  {
    MultiFluidBase * fluid = GetFluidModel( subRegion );
    fluid->setMassFlag( static_cast<bool>(m_useMass) );
  });

  // Initialize primary variables from applied initial conditions
  ResetViews( domain );
  InitializeFluidState( domain );
}

real64 CompositionalMultiphaseFlow::SolverStep( real64 const & time_n,
                                                real64 const & dt,
                                                integer const cycleNumber,
                                                DomainPartition * const domain )
{
  real64 dt_return = dt;

  ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );

  // currently the only method is implicit time integration
  dt_return= this->NonlinearImplicitStep( time_n,
                                          dt,
                                          cycleNumber,
                                          domain,
                                          getLinearSystemRepository() );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

void CompositionalMultiphaseFlow::BackupFields( DomainPartition * const domain )
{
  // backup some fields used in time derivative approximation
  applyToSubRegions( domain, [&] ( localIndex er, localIndex esr,
                                   ElementRegion * const region,
                                   CellBlockSubRegion * const subRegion )
  {
    arrayView1d<integer const> const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<real64  const> const & poroRef       = m_porosityRef[er][esr];

    arrayView2d<real64  const> const & phaseVolFrac = m_phaseVolFrac[er][esr];

    arrayView3d<real64  const> const & phaseDens     = m_phaseDens[er][esr][m_fluidIndex];
    arrayView4d<real64  const> const & phaseCompFrac = m_phaseCompFrac[er][esr][m_fluidIndex];
    arrayView2d<real64  const> const & pvMult        = m_pvMult[er][esr][m_solidIndex];

    arrayView2d<real64> const & phaseDensOld     = m_phaseDensOld[er][esr];
    arrayView2d<real64> const & phaseVolFracOld  = m_phaseVolFracOld[er][esr];
    arrayView3d<real64> const & phaseCompFracOld = m_phaseCompFracOld[er][esr];
    arrayView1d<real64> const & poroOld          = m_porosityOld[er][esr];

    for_elems_in_subRegion( subRegion, GEOSX_LAMBDA ( localIndex const ei )
    {
      if (elemGhostRank[ei] >= 0)
        return;

      for (localIndex ip = 0; ip < m_numPhases; ++ip)
      {
        phaseDensOld[ei][ip] = phaseDens[ei][0][ip];
        phaseVolFracOld[ei][ip] = phaseVolFrac[ei][ip];

        for (localIndex ic = 0; ic < m_numComponents; ++ic)
        {
          phaseCompFracOld[ei][ip][ic] = phaseCompFrac[ei][0][ip][ic];
        }
      }

      poroOld[ei] = poroRef[ei] * pvMult[ei][0];
    });
  });
}

void
CompositionalMultiphaseFlow::ImplicitStepSetup( real64 const & time_n, real64 const & dt,
                                                DomainPartition * const domain,
                                                EpetraBlockSystem * const blockSystem )
{
  // bind the stored views to the current domain
  ResetViews( domain );

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
  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & dofNumber     = m_dofNumber;
  ElementRegionManager::ElementViewAccessor<arrayView1d<integer>>     const & elemGhostRank = m_elemGhostRank;

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

  GEOS_ERROR_IF( numLocalRows != numLocalRowsToSend, "number of local rows inconsistent" );

  // find the first local row on this partition, and find the number of total global rows.
  localIndex firstLocalRow = 0;
  numGlobalRows = 0;

  for( integer p = 0 ; p < numMpiProcesses ; ++p)
  {
    numGlobalRows += gather[p];
    if (p < thisMpiProcess)
      firstLocalRow += gather[p];
  }

  // create trilinos dof indexing, setting initial values to -1 to indicate unset values.
  for( localIndex er=0 ; er < elemGhostRank.size() ; ++er )
  {
    for( localIndex esr=0 ; esr < elemGhostRank[er].size() ; ++esr )
    {
      dofNumber[er][esr] = -1;
    }
  }

  // loop over all elements and set the dof number if the element is not a ghost
  ReduceSum< reducePolicy, localIndex  > localCount(0);
  forAllElemsInMesh<RAJA::seq_exec>( meshLevel, GEOSX_LAMBDA ( localIndex const er,
                                                               localIndex const esr,
                                                               localIndex const ei )
  {
    if( elemGhostRank[er][esr][ei] < 0 )
    {
      dofNumber[er][esr][ei] = firstLocalRow + localCount + offset;
      localCount += 1;
    }
  });

  GEOS_ERROR_IF( localCount != numLocalRows, "Number of DOF assigned does not match numLocalRows" );
}

void CompositionalMultiphaseFlow::SetSparsityPattern( DomainPartition const * const domain,
                                                      Epetra_FECrsGraph * const sparsity )
{
  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elementRegionManager = meshLevel->getElemManager();

  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & dofNumber =
    elementRegionManager->ConstructViewAccessor<array1d<globalIndex>, arrayView1d<globalIndex>>( viewKeyStruct::blockLocalDofNumberString );

  ElementRegionManager::ElementViewAccessor<arrayView1d<integer>> const & elemGhostRank =
    elementRegionManager->ConstructViewAccessor<array1d<integer>, arrayView1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

  NumericalMethodsManager const * numericalMethodManager = domain->
    getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager const * fvManager = numericalMethodManager->
    GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation( m_discretizationName );
  FluxApproximationBase::CellStencil const & stencilCollection = fluxApprox->getStencil();

  localIndex constexpr numElems   = StencilCollection<CellDescriptor, real64>::NUM_POINT_IN_FLUX;
  localIndex constexpr maxStencil = StencilCollection<CellDescriptor, real64>::MAX_STENCIL_SIZE;
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex const NDOF = m_numDofPerCell;

  //**** loop over all faces. Fill in sparsity for all pairs of DOF/elem that are connected by face
  stencilCollection.forAll( GEOSX_LAMBDA ( StencilCollection<CellDescriptor, real64>::Accessor stencil )
  {
    localIndex const stencilSize = stencil.size();
    stackArray1d<globalIndex, numElems   * maxNumDof> elementLocalDofIndexRow( numElems * NDOF );
    stackArray1d<globalIndex, maxStencil * maxNumDof> elementLocalDofIndexCol( stencilSize * NDOF );

    stencil.forConnected( [&] ( CellDescriptor const & cell, localIndex const i )
    {
      globalIndex const offset = NDOF * dofNumber[cell.region][cell.subRegion][cell.index];
      for (localIndex idof = 0; idof < NDOF; ++idof)
      {
        elementLocalDofIndexRow[i * NDOF + idof] = offset + idof;
      }
    } );

    stencil.forAll( [&] ( CellDescriptor const & cell, real64 w, localIndex const i )
    {
      globalIndex const offset = NDOF * dofNumber[cell.region][cell.subRegion][cell.index];
      for (localIndex idof = 0; idof < NDOF; ++idof)
      {
        elementLocalDofIndexCol[i * NDOF + idof] = offset + idof;
      }
    } );

    sparsity->InsertGlobalIndices( integer_conversion<int>(numElems * m_numDofPerCell),
                                   elementLocalDofIndexRow.data(),
                                   integer_conversion<int>(stencilSize * m_numDofPerCell),
                                   elementLocalDofIndexCol.data() );
  } );

  // loop over all elements and add all locals just in case the above connector loop missed some
  forAllElemsInMesh( meshLevel, [&] ( localIndex const er,
                                      localIndex const esr,
                                      localIndex const ei )
  {
    if (elemGhostRank[er][esr][ei] < 0)
    {
      stackArray1d<globalIndex, maxNumDof> elementLocalDofIndexRow( NDOF );

      globalIndex const offset = NDOF * dofNumber[er][esr][ei];
      for (localIndex idof = 0; idof < NDOF; ++idof)
      {
        elementLocalDofIndexRow[idof] = offset + idof;
      }

      sparsity->InsertGlobalIndices( integer_conversion<int>( NDOF ),
                                     elementLocalDofIndexRow.data(),
                                     integer_conversion<int>( NDOF ),
                                     elementLocalDofIndexRow.data() );
    }
  });

  // add additional connectivity resulting from boundary stencils
  fluxApprox->forBoundaryStencils( [&] ( FluxApproximationBase::BoundaryStencil const & boundaryStencilCollection )
  {
    boundaryStencilCollection.forAll( GEOSX_LAMBDA ( StencilCollection<PointDescriptor, real64>::Accessor stencil )
    {
      localIndex const stencilSize = stencil.size();
      stackArray1d<globalIndex, numElems   * maxNumDof> elementLocalDofIndexRow( numElems * NDOF );
      stackArray1d<globalIndex, maxStencil * maxNumDof> elementLocalDofIndexCol( stencilSize * NDOF );

      stencil.forConnected( [&] ( PointDescriptor const & point, localIndex const i )
      {
        if (point.tag == PointDescriptor::Tag::CELL)
        {
          CellDescriptor const & cell = point.cellIndex;
          globalIndex const offset = NDOF * dofNumber[cell.region][cell.subRegion][cell.index];
          for (localIndex idof = 0; idof < NDOF; ++idof)
          {
            elementLocalDofIndexRow[idof] = offset + idof;
          }
        }
      });

      integer counter = 0;
      stencil.forAll( [&] ( PointDescriptor const & point, real64 w, localIndex i )
      {
        if (point.tag == PointDescriptor::Tag::CELL)
        {
          CellDescriptor const & cell = point.cellIndex;
          globalIndex const offset = NDOF * dofNumber[cell.region][cell.subRegion][cell.index];
          for (localIndex idof = 0; idof < NDOF; ++idof)
          {
            elementLocalDofIndexCol[counter * NDOF + idof] = offset + idof;
          }
          ++counter;
        }
      });

      sparsity->InsertGlobalIndices( integer_conversion<int>( NDOF ),
                                     elementLocalDofIndexRow.data(),
                                     integer_conversion<int>( counter * NDOF ),
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

  localIndex_array displacementIndices;
  SetNumRowsAndTrilinosIndices( mesh,
                                numLocalRows,
                                numGlobalRows,
                                0 );

  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::blockLocalDofNumberString );
  CommunicationTools::
  SynchronizeFields(fieldNames,
                    mesh,
                    domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );


  // construct row map, and set a pointer to the row map
  Epetra_Map * const rowMap =
    blockSystem->SetRowMap( BlockIDs::compositionalBlock,
                            std::make_unique<Epetra_Map>( numGlobalRows * m_numDofPerCell,
                                                          numLocalRows * m_numDofPerCell,
                                                          0,
                                                          m_linearSolverWrapper.m_epetraComm ) );

  // construct sparsity matrix, set a pointer to the sparsity pattern matrix
  Epetra_FECrsGraph * const sparsity =
    blockSystem->SetSparsity( BlockIDs::compositionalBlock,
                              BlockIDs::compositionalBlock,
                              std::make_unique<Epetra_FECrsGraph>(Copy, *rowMap, 0) );



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
  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::compositionalBlock,
                                                                BlockIDs::compositionalBlock );
  Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );

  jacobian->Scale(0.0);
  residual->Scale(0.0);

  AssembleAccumulationTerms( domain, jacobian, residual, time_n, dt );
  AssembleFluxTerms( domain, jacobian, residual, time_n, dt );
  AssembleVolumeBalanceTerms( domain, jacobian, residual, time_n, dt );

  jacobian->GlobalAssemble( true );
  residual->GlobalAssemble();

  if( verboseLevel() >= 3 )
  {
    GEOS_LOG_RANK("After CompositionalMultiphaseFlow::AssembleSystem");
    GEOS_LOG_RANK("\nJacobian:\n" << *jacobian);
    GEOS_LOG_RANK("\nResidual:\n" << *residual);
  }

}

void CompositionalMultiphaseFlow::AssembleAccumulationTerms( DomainPartition const * const domain,
                                                             Epetra_FECrsMatrix * const jacobian,
                                                             Epetra_FEVector * const residual,
                                                             real64 const time_n,
                                                             real64 const dt )
{
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex const NC   = m_numComponents;
  localIndex const NP   = m_numPhases;
  localIndex const NDOF = m_numDofPerCell;

  applyToSubRegions( domain, [&] ( localIndex er, localIndex esr,
                                   ElementRegion const * const region,
                                   CellBlockSubRegion const * const subRegion )
  {
    arrayView1d<integer     const> const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<globalIndex const> const & dofNumber     = m_dofNumber[er][esr];

    arrayView1d<real64 const> const & volume      = m_volume[er][esr];
    arrayView1d<real64 const> const & porosityRef = m_porosityRef[er][esr];

    arrayView2d<real64 const> const & phaseVolFrac            = m_phaseVolFrac[er][esr];
    arrayView2d<real64 const> const & dPhaseVolFrac_dPres     = m_dPhaseVolFrac_dPres[er][esr];
    arrayView3d<real64 const> const & dPhaseVolFrac_dCompDens = m_dPhaseVolFrac_dCompDens[er][esr];
    arrayView3d<real64 const> const & dCompFrac_dCompDens     = m_dCompFrac_dCompDens[er][esr];

    arrayView1d<real64 const> const & porosityOld      = m_porosityOld[er][esr];
    arrayView2d<real64 const> const & phaseVolFracOld  = m_phaseVolFracOld[er][esr];
    arrayView2d<real64 const> const & phaseDensOld     = m_phaseDensOld[er][esr];
    arrayView3d<real64 const> const & phaseCompFracOld = m_phaseCompFracOld[er][esr];

    arrayView2d<real64 const> const & pvMult               = m_pvMult[er][esr][m_solidIndex];
    arrayView2d<real64 const> const & dPvMult_dPres        = m_dPvMult_dPres[er][esr][m_solidIndex];
    arrayView3d<real64 const> const & phaseDens            = m_phaseDens[er][esr][m_fluidIndex];
    arrayView3d<real64 const> const & dPhaseDens_dPres     = m_dPhaseDens_dPres[er][esr][m_fluidIndex];
    arrayView4d<real64 const> const & dPhaseDens_dComp     = m_dPhaseDens_dComp[er][esr][m_fluidIndex];
    arrayView4d<real64 const> const & phaseCompFrac        = m_phaseCompFrac[er][esr][m_fluidIndex];
    arrayView4d<real64 const> const & dPhaseCompFrac_dPres = m_dPhaseCompFrac_dPres[er][esr][m_fluidIndex];
    arrayView5d<real64 const> const & dPhaseCompFrac_dComp = m_dPhaseCompFrac_dComp[er][esr][m_fluidIndex];

    for_elems_in_subRegion( subRegion, GEOSX_LAMBDA ( localIndex const ei )
    {
      if (elemGhostRank[ei] >= 0)
        return;

      stackArray1d<long long, maxNumDof>           localAccumDOF( NDOF );
      stackArray1d<double, maxNumComp>             localAccum( NC );
      stackArray2d<double, maxNumComp * maxNumDof> localAccumJacobian( NC, NDOF );

      // temporary work arrays
      stackArray1d<double, maxNumComp> dPhaseAmount_dC( NC );
      stackArray1d<double, maxNumComp> dPhaseCompFrac_dC( NC );

      // reset the local values
      localAccum = 0.0;
      localAccumJacobian = 0.0;

      // set DOF indices for this block
      globalIndex const offset = NDOF * dofNumber[ei];
      for (localIndex idof = 0; idof < NDOF; ++idof)
      {
        localAccumDOF[idof] = integer_conversion<long long>(offset + idof);
      }

      // compute fluid-independent (pore volume) part
      real64 const volNew   = volume[ei];
      real64 const volOld   = volume[ei];
      real64 const dVol_dP  = 0.0; // used in poroelastic solver

      real64 const poroNew  = porosityRef[ei] * pvMult[ei][0];
      real64 const poroOld  = porosityOld[ei];
      real64 const dPoro_dP = porosityRef[ei] * dPvMult_dPres[ei][0];

      real64 const poreVolNew = volNew * poroNew;
      real64 const poreVolOld = volOld * poroOld;
      real64 const dPoreVol_dP = dVol_dP * poroNew + volNew * dPoro_dP;

      // sum contributions to component accumulation from each phase
      for (localIndex ip = 0; ip < NP; ++ip)
      {
        real64 const phaseAmountNew = poreVolNew * phaseVolFrac[ei][ip] * phaseDens[ei][0][ip];
        real64 const phaseAmountOld = poreVolOld * phaseVolFracOld[ei][ip] * phaseDensOld[ei][ip];

        real64 const dPhaseAmount_dP = dPoreVol_dP * phaseVolFrac[ei][ip] * phaseDens[ei][0][ip]
                                      + poreVolNew * (dPhaseVolFrac_dPres[ei][ip] * phaseDens[ei][0][ip]
                                                           + phaseVolFrac[ei][ip] * dPhaseDens_dPres[ei][0][ip]);

        // assemble density dependence
        applyChainRule( NC, dCompFrac_dCompDens[ei], dPhaseDens_dComp[ei][0][ip], dPhaseAmount_dC );
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc]     * phaseVolFrac[ei][ip]
                              + phaseDens[ei][0][ip] * dPhaseVolFrac_dCompDens[ei][ip][jc];
          dPhaseAmount_dC[jc] *= poreVolNew;
        }

        // ic - index of component whose conservation equation is assembled
        // (i.e. row number in local matrix)
        for (localIndex ic = 0; ic < NC; ++ic)
        {
          real64 const phaseCompAmountNew = phaseAmountNew * phaseCompFrac[ei][0][ip][ic];
          real64 const phaseCompAmountOld = phaseAmountOld * phaseCompFracOld[ei][ip][ic];

          real64 const dPhaseCompAmount_dP = dPhaseAmount_dP * phaseCompFrac[ei][0][ip][ic]
                                           + phaseAmountNew  * dPhaseCompFrac_dPres[ei][0][ip][ic];

          localAccum[ic] += phaseCompAmountNew - phaseCompAmountOld;
          localAccumJacobian[ic][0] += dPhaseCompAmount_dP;

          // jc - index of component w.r.t. whose compositional var the derivative is being taken
          // (i.e. col number in local matrix)

          // assemble phase composition dependence
          applyChainRule( NC, dCompFrac_dCompDens[ei], dPhaseCompFrac_dComp[ei][0][ip][ic], dPhaseCompFrac_dC );
          for (localIndex jc = 0; jc < NC; ++jc)
          {
            real64 const dPhaseCompAmount_dC =           dPhaseCompFrac_dC[jc] * phaseAmountNew
                                             + phaseCompFrac[ei][0][ip][ic] * dPhaseAmount_dC[jc];
            localAccumJacobian[ic][jc+1] += dPhaseCompAmount_dC;
          }
        }
      }

      // TODO: apply equation/variable change transformation(s)

      // add contribution to global residual and dRdP
      residual->SumIntoGlobalValues( integer_conversion<int>( NC ),
                                     localAccumDOF.data(),
                                     localAccum.data() );

      jacobian->SumIntoGlobalValues( integer_conversion<int>( NC ),
                                     localAccumDOF.data(),
                                     integer_conversion<int>( NDOF ),
                                     localAccumDOF.data(),
                                     localAccumJacobian.data(),
                                     Epetra_FECrsMatrix::ROW_MAJOR );
    });
  });
}

void CompositionalMultiphaseFlow::AssembleFluxTerms( DomainPartition const * const domain,
                                                     Epetra_FECrsMatrix * const jacobian,
                                                     Epetra_FEVector * const residual,
                                                     real64 const time_n,
                                                     real64 const dt )
{
  NumericalMethodsManager const * const numericalMethodManager =
    domain->getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager const * const fvManager =
    numericalMethodManager->GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  FluxApproximationBase const * const fluxApprox = fvManager->getFluxApproximation( m_discretizationName );
  FluxApproximationBase::CellStencil const & stencilCollection = fluxApprox->getStencil();

  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & blockLocalDofNumber = m_dofNumber;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & pres = m_pressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & dPres = m_deltaPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & gravDepth = m_gravDepth;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & dPhaseVolFrac_dPres = m_dPhaseVolFrac_dPres;
  ElementRegionManager::ElementViewAccessor<arrayView3d<real64>> const & dPhaseVolFrac_dComp = m_dPhaseVolFrac_dCompDens;
  ElementRegionManager::ElementViewAccessor<arrayView3d<real64>> const & dCompFrac_dCompDens = m_dCompFrac_dCompDens;

  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & phaseDens        = m_phaseDens;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dPhaseDens_dPres = m_dPhaseDens_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & dPhaseDens_dComp = m_dPhaseDens_dComp;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & phaseVisc        = m_phaseVisc;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dPhaseVisc_dPres = m_dPhaseVisc_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & dPhaseVisc_dComp = m_dPhaseVisc_dComp;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & phaseCompFrac    = m_phaseCompFrac;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & dPhaseCompFrac_dPres = m_dPhaseCompFrac_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView5d<real64>> const & dPhaseCompFrac_dComp = m_dPhaseCompFrac_dComp;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & phaseRelPerm                = m_phaseRelPerm;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & dPhaseRelPerm_dPhaseVolFrac = m_dPhaseRelPerm_dPhaseVolFrac;

  localIndex constexpr numElems   = StencilCollection<CellDescriptor, real64>::NUM_POINT_IN_FLUX;
  localIndex constexpr maxStencil = StencilCollection<CellDescriptor, real64>::MAX_STENCIL_SIZE;
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex const NC   = m_numComponents;
  localIndex const NP   = m_numPhases;
  localIndex const NDOF = m_numDofPerCell;

  real64 constexpr densWeight[numElems] = { 0.5, 0.5 };

  stencilCollection.forAll<RAJA::seq_exec>(GEOSX_LAMBDA (StencilCollection<CellDescriptor, real64>::Accessor stencil)
  {
    localIndex const stencilSize = stencil.size();

    // create local work arrays
    stackArray1d<long long, numElems * maxNumComp>  eqnRowIndices( numElems * NC );
    stackArray1d<long long, maxStencil * maxNumDof> dofColIndices( stencilSize * NDOF );

    stackArray1d<double, numElems * maxNumComp>                          localFlux( numElems * NC );
    stackArray2d<double, numElems * maxNumComp * maxStencil * maxNumDof> localFluxJacobian( numElems * NC, stencilSize * NDOF );

    stackArray1d<real64, maxNumComp> dPhaseCompFrac_dCompDens( NC );

    stackArray1d<real64, maxStencil>              dPhaseFlux_dP( stencilSize );
    stackArray2d<real64, maxStencil * maxNumComp> dPhaseFlux_dC( stencilSize, NC );

    stackArray1d<real64, maxNumComp>                           compFlux( NC );
    stackArray2d<real64, maxStencil * maxNumComp>              dCompFlux_dP( stencilSize, NC );
    stackArray3d<real64, maxStencil * maxNumComp * maxNumComp> dCompFlux_dC( stencilSize, NC, NC );

    stackArray1d<real64, maxNumComp> dRelPerm_dC( NC );
    stackArray1d<real64, maxNumComp> dDens_dC( NC );
    stackArray1d<real64, maxNumComp> dVisc_dC( NC );

    stackArray1d<real64, numElems>              mobility( numElems );
    stackArray1d<real64, numElems>              dMobility_dP( numElems );
    stackArray2d<real64, numElems * maxNumComp> dMobility_dC( numElems, NC );

    stackArray1d<real64, numElems>              dDensMean_dP( numElems );
    stackArray2d<real64, numElems * maxNumComp> dDensMean_dC( numElems, NC );

    stackArray1d<real64, maxStencil>            dPresGrad_dP( stencilSize );

    stackArray1d<real64, numElems>              dGravHead_dP( numElems );
    stackArray2d<real64, numElems * maxNumComp> dGravHead_dC( numElems, NC );

    // reset the local values
    compFlux = 0.0;
    dCompFlux_dP = 0.0;
    dCompFlux_dC = 0.0;

    localFlux = 0.0;
    localFluxJacobian = 0.0;

    // set equation indices for both connected cells
    stencil.forConnected( [&] ( auto const & cell,
                                localIndex i )
    {
      localIndex const er  = cell.region;
      localIndex const esr = cell.subRegion;
      localIndex const ei  = cell.index;

      globalIndex const offset = NDOF * blockLocalDofNumber[er][esr][ei];
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        eqnRowIndices[i * NC + ic] = offset + ic;
      }
    });

    // loop over phases, compute and upwind phase flux and sum contributions to each component's flux
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      // clear working arrays
      real64 densMean = 0.0;
      dDensMean_dP = 0.0;
      dDensMean_dC = 0.0;

      real64 presGrad = 0.0;
      dPresGrad_dP = 0.0;

      real64 gravHead = 0.0;
      dGravHead_dP = 0.0;
      dGravHead_dC = 0.0;

      real64 phaseFlux = 0.0;
      dPhaseFlux_dP = 0.0;
      dPhaseFlux_dC = 0.0;

      // calculate quantities on primary connected cells
      stencil.forConnected( [&] ( auto const & cell,
                                 localIndex i )
      {
        localIndex const er  = cell.region;
        localIndex const esr = cell.subRegion;
        localIndex const ei  = cell.index;

        // density
        real64 const density  = phaseDens[er][esr][m_fluidIndex][ei][0][ip];
        real64 const dDens_dP = dPhaseDens_dPres[er][esr][m_fluidIndex][ei][0][ip];
        applyChainRule( NC, dCompFrac_dCompDens[er][esr][ei], dPhaseDens_dComp[er][esr][m_fluidIndex][ei][0][ip], dDens_dC );

        // viscosity
        real64 const viscosity = phaseVisc[er][esr][m_fluidIndex][ei][0][ip];
        real64 const dVisc_dP  = dPhaseVisc_dPres[er][esr][m_fluidIndex][ei][0][ip];
        applyChainRule( NC, dCompFrac_dCompDens[er][esr][ei], dPhaseVisc_dComp[er][esr][m_fluidIndex][ei][0][ip], dVisc_dC );

        //relative permeability
        real64 const relPerm = phaseRelPerm[er][esr][m_relPermIndex][ei][0][ip];
        real64 dRelPerm_dP = 0.0;
        dRelPerm_dC = 0.0;
        for (localIndex jp = 0; jp < NP; ++jp)
        {
          real64 const dRelPerm_dS = dPhaseRelPerm_dPhaseVolFrac[er][esr][m_relPermIndex][ei][0][ip][jp];
          dRelPerm_dP += dRelPerm_dS * dPhaseVolFrac_dPres[er][esr][ei][jp];

          for (localIndex jc = 0; jc < NC; ++jc)
          {
            dRelPerm_dC[jc] += dRelPerm_dS * dPhaseVolFrac_dComp[er][esr][ei][jp][jc];
          }
        }

        // mobility and pressure derivative
        mobility[i] = relPerm * density / viscosity;
        dMobility_dP[i] = dRelPerm_dP * density / viscosity
                        + mobility[i] * (dDens_dP / density - dVisc_dP / viscosity);

        // average density and pressure derivative
        densMean += densWeight[i] * density;
        dDensMean_dP[i] = densWeight[i] * dDens_dP;

        // compositional derivatives
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dDensMean_dC[i][jc] = densWeight[i] * dDens_dC[jc];

          dMobility_dC[i][jc] = dRelPerm_dC[jc] * density / viscosity
                              + mobility[i] * (dDens_dC[jc] / density - dVisc_dC[jc] / viscosity);
        }
      });

      //***** calculation of flux *****

      // compute potential difference MPFA-style
      stencil.forAll( [&] ( CellDescriptor cell,
                            real64 w,
                            localIndex i )
      {
        localIndex const er  = cell.region;
        localIndex const esr = cell.subRegion;
        localIndex const ei  = cell.index;

        globalIndex const offset = NDOF * blockLocalDofNumber[er][esr][ei];
        for (localIndex jdof = 0; jdof < NDOF; ++jdof)
        {
          dofColIndices[i * NDOF + jdof] = offset + jdof;
        }

        presGrad += w * (pres[er][esr][ei] + dPres[er][esr][ei]);
        dPresGrad_dP[i] += w;

        if (m_gravityFlag)
        {
          real64 const gravD = w * gravDepth[er][esr][ei];
          gravHead += densMean * gravD;

          // need to add contributions from both cells the mean density depends on
          stencil.forConnected( [&] ( CellDescriptor,
                                      localIndex j )
          {
            dGravHead_dP[j] += dDensMean_dP[j] * gravD;
            for (localIndex jc = 0; jc < NC; ++jc)
            {
              dGravHead_dC[j][jc] += dDensMean_dC[j][jc] * gravD;
            }
          });
        }
      });

      // *** upwinding ***

      // use PPU currently; advanced stuff like IHU would go here
      // TODO isolate into a kernel?

      // compute phase potential gradient
      real64 potGrad = presGrad + gravHead;

      // choose upstream cell
      localIndex const k_up = (potGrad >= 0) ? 0 : 1;

      // skip the phase flux if phase not present or immobile upstream
      if (std::fabs(mobility[k_up]) < 1e-20) // TODO better constant
      {
        break;
      }

      // pressure gradient depends on all points in the stencil
      for (localIndex ke = 0; ke < stencilSize; ++ke)
      {
        dPhaseFlux_dP[ke] += dPresGrad_dP[ke];
      }

      // gravitational head depends only on the two cells connected (same as mean density)
      for (localIndex ke = 0; ke < numElems; ++ke)
      {
        dPhaseFlux_dP[ke] += dGravHead_dP[ke];
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dPhaseFlux_dC[ke][jc] += dGravHead_dC[ke][jc];
        }
      }

      // compute the phase flux and derivatives using upstream cell mobility
      phaseFlux = mobility[k_up] * potGrad;
      for (localIndex ke = 0; ke < stencilSize; ++ke)
      {
        dPhaseFlux_dP[ke] *= mobility[k_up];
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dPhaseFlux_dC[ke][jc] *= mobility[k_up];
        }
      }

      // add contribution from upstream cell mobility derivatives
      dPhaseFlux_dP[k_up] += dMobility_dP[k_up] * potGrad;
      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dPhaseFlux_dC[k_up][jc] += dMobility_dC[k_up][jc] * potGrad;
      }

      // get global identifiers of the upstream cell
      CellDescriptor cell_up = stencil.connectedIndex( k_up );
      localIndex er_up  = cell_up.region;
      localIndex esr_up = cell_up.subRegion;
      localIndex ei_up  = cell_up.index;

      // slice some constitutive arrays to avoid too much indexing in component loop
      arraySlice1d<real64> phaseCompFracSub = phaseCompFrac[er_up][esr_up][m_fluidIndex][ei_up][0][ip];
      arraySlice1d<real64> dPhaseCompFrac_dPresSub = dPhaseCompFrac_dPres[er_up][esr_up][m_fluidIndex][ei_up][0][ip];
      arraySlice2d<real64> dPhaseCompFrac_dCompSub = dPhaseCompFrac_dComp[er_up][esr_up][m_fluidIndex][ei_up][0][ip];

      // compute component fluxes and derivatives using upstream cell composition
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        real64 const ycp = phaseCompFracSub[ic];
        compFlux[ic] += phaseFlux * ycp;

        // derivatives stemming from phase flux
        for (localIndex ke = 0; ke < stencilSize; ++ke)
        {
          dCompFlux_dP[ke][ic] += dPhaseFlux_dP[ke] * ycp;
          for (localIndex jc = 0; jc < NC; ++jc)
          {
            dCompFlux_dC[ke][ic][jc] += dPhaseFlux_dC[ke][jc] * ycp;
          }
        }

        // additional derivatives stemming from upstream cell phase composition
        dCompFlux_dP[k_up][ic] += phaseFlux * dPhaseCompFrac_dPresSub[ic];

        // convert derivatives of component fraction w.r.t. component fractions to derivatives w.r.t. component densities
        applyChainRule( NC, dCompFrac_dCompDens[er_up][esr_up][ei_up], dPhaseCompFrac_dCompSub[ic], dPhaseCompFrac_dCompDens );
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dCompFlux_dC[k_up][ic][jc] += phaseFlux * dPhaseCompFrac_dCompDens[jc];
        }
      }
    }

    // *** end of upwinding

    // populate local flux vector and derivatives
    for (localIndex ic = 0; ic < NC; ++ic)
    {
      localFlux[ic]      =  dt * compFlux[ic];
      localFlux[NC + ic] = -dt * compFlux[ic];

      for (localIndex ke = 0; ke < stencilSize; ++ke)
      {
        localIndex const localDofIndexPres = ke * NDOF;
        localFluxJacobian[ic][localDofIndexPres] = dt * dCompFlux_dP[ke][ic];
        localFluxJacobian[NC + ic][localDofIndexPres] = -dt * dCompFlux_dP[ke][ic];

        for (localIndex jc = 0; jc < NC; ++jc)
        {
          localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
          localFluxJacobian[ic][localDofIndexComp] = dt * dCompFlux_dC[ke][ic][jc];
          localFluxJacobian[NC + ic][localDofIndexComp] = -dt * dCompFlux_dC[ke][ic][jc];
        }
      }
    }

    // TODO: apply equation/variable change transformation(s)

    // Add to global residual/jacobian
    residual->SumIntoGlobalValues( integer_conversion<int>( numElems * NC ),
                                   eqnRowIndices.data(),
                                   localFlux.data() );

    jacobian->SumIntoGlobalValues( integer_conversion<int>( numElems * NC ),
                                   eqnRowIndices.data(),
                                   integer_conversion<int>( stencilSize * NDOF ),
                                   dofColIndices.data(),
                                   localFluxJacobian.data(),
                                   Epetra_FECrsMatrix::ROW_MAJOR );

  });
}

void CompositionalMultiphaseFlow::AssembleVolumeBalanceTerms( DomainPartition const * const domain,
                                                              Epetra_FECrsMatrix * const jacobian,
                                                              Epetra_FEVector * const residual,
                                                              real64 const time_n,
                                                              real64 const dt )
{
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex const NC   = m_numComponents;
  localIndex const NP   = m_numPhases;
  localIndex const NDOF = m_numDofPerCell;

  applyToSubRegions( domain, [&] ( localIndex er, localIndex esr,
                                   ElementRegion const * const region,
                                   CellBlockSubRegion const * const subRegion )
  {
    arrayView1d<integer const>     const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<globalIndex const> const & dofNumber     = m_dofNumber[er][esr];

    arrayView1d<real64 const> const & volume      = m_volume[er][esr];
    arrayView1d<real64 const> const & porosityRef = m_porosityRef[er][esr];

    arrayView2d<real64 const> const & phaseVolFrac            = m_phaseVolFrac[er][esr];
    arrayView2d<real64 const> const & dPhaseVolFrac_dPres     = m_dPhaseVolFrac_dPres[er][esr];
    arrayView3d<real64 const> const & dPhaseVolFrac_dCompDens = m_dPhaseVolFrac_dCompDens[er][esr];

    arrayView2d<real64 const> const & pvMult        = m_pvMult[er][esr][m_solidIndex];
    arrayView2d<real64 const> const & dPvMult_dPres = m_dPvMult_dPres[er][esr][m_solidIndex];

    for_elems_in_subRegion( subRegion, GEOSX_LAMBDA ( localIndex const ei )
    {
      if (elemGhostRank[ei] >= 0)
        return;

      stackArray1d<long long, maxNumDof> localVolBalanceDOF( NDOF );
      stackArray1d<double, maxNumDof>    localVolBalanceJacobian( NDOF );

      // compute pore volume
      real64 const vol      = volume[ei];
      real64 const dVol_dP  = 0.0; // used in poroelastic solver

      real64 const poro     = porosityRef[ei] * pvMult[ei][0];
      real64 const dPoro_dP = porosityRef[ei] * dPvMult_dPres[ei][0];

      real64 const poreVol     = vol * poro;
      real64 const dPoreVol_dP = dVol_dP * poro + vol * dPoro_dP;

      // get equation/dof indices
      globalIndex const offset = NDOF * dofNumber[ei];
      globalIndex const localVolBalanceEqnIndex = offset + NC;
      for (localIndex jdof = 0; jdof < NDOF; ++jdof)
      {
        localVolBalanceDOF[jdof] = offset + jdof;
      }

      real64 localVolBalance = 1.0;
      localVolBalanceJacobian = 0.0;

      // sum contributions to component accumulation from each phase
      for (localIndex ip = 0; ip < NP; ++ip)
      {
        localVolBalance -= phaseVolFrac[ei][ip];
        localVolBalanceJacobian[0] -= dPhaseVolFrac_dPres[ei][ip];

        for (localIndex jc = 0; jc < NC; ++jc)
        {
          localVolBalanceJacobian[jc+1] -= dPhaseVolFrac_dCompDens[ei][ip][jc];
        }
      }

      // scale saturation-based volume balance by pore volume (for better scaling w.r.t. other equations)
      for (localIndex idof = 0; idof < NDOF; ++idof)
      {
        localVolBalanceJacobian[idof] *= poreVol;
      }
      localVolBalanceJacobian[0] += dPoreVol_dP * localVolBalance;
      localVolBalance *= poreVol;

      // TODO: apply equation/variable change transformation(s)

      // add contribution to global residual and dRdP
      residual->SumIntoGlobalValues( 1,
                                     &localVolBalanceEqnIndex,
                                     &localVolBalance );

      jacobian->SumIntoGlobalValues( 1,
                                     &localVolBalanceEqnIndex,
                                     integer_conversion<int>( NDOF ),
                                     localVolBalanceDOF.data(),
                                     localVolBalanceJacobian.data(),
                                     Epetra_FECrsMatrix::ROW_MAJOR );
    });
  });
}

void CompositionalMultiphaseFlow::ApplyBoundaryConditions( DomainPartition * const domain,
                                                           EpetraBlockSystem * const blockSystem,
                                                           real64 const time_n, real64 const dt )
{
  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::compositionalBlock,
                                                                BlockIDs::compositionalBlock );
  Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );

  // apply pressure boundary conditions.
  ApplyDirichletBC_implicit(domain, time_n, dt, blockSystem);
  ApplyFaceDirichletBC_implicit(domain, time_n, dt, blockSystem);

  jacobian->GlobalAssemble();
  residual->GlobalAssemble();

  if (verboseLevel() >= 2)
  {
    GEOS_LOG_RANK( "After CompositionalMultiphaseFlow::ApplyBoundaryCondition" );
    GEOS_LOG_RANK( "\nJacobian\n" << *jacobian );
    GEOS_LOG_RANK( "\nResidual\n" << *residual );
  }
}

void
CompositionalMultiphaseFlow::ApplyDirichletBC_implicit( DomainPartition * const domain,
                                                        real64 const time_n, real64 const dt,
                                                        EpetraBlockSystem * const blockSystem )
{
  BoundaryConditionManager * bcManager = BoundaryConditionManager::get();

  unordered_map<string, array1d<bool>> bcStatusMap; // map to check consistent application of BC

  // 1. apply pressure Dirichlet BCs
  bcManager->ApplyBoundaryCondition( time_n + dt,
                                     domain,
                                     "ElementRegions",
                                     viewKeyStruct::pressureString,
                                     [&]( BoundaryConditionBase const * const bc,
                                          string const & setName,
                                          set<localIndex> const & targetSet,
                                          ManagedGroup * subRegion,
                                          string const & )
  {
    // 1.0. Check whether pressure has already been applied to this set
    GEOS_ERROR_IF( bcStatusMap.count( setName ) > 0, "Conflicting pressure boundary conditions on set " << setName );
    bcStatusMap[setName].resize( m_numComponents );
    bcStatusMap[setName] = false;

    // 1.1. Apply BC to set the field values
    bc->ApplyBoundaryConditionToField<BcEqual>( targetSet,
                                                time_n + dt,
                                                subRegion,
                                                viewKeyStruct::bcPressureString );
  });

  // 2. Apply composition BC (global component fraction) and store them for constitutive call
  bcManager->ApplyBoundaryCondition( time_n + dt,
                                     domain,
                                     "ElementRegions",
                                     viewKeyStruct::globalCompFractionString,
                                     [&] ( BoundaryConditionBase const * const bc,
                                           string const & setName,
                                           set<localIndex> const & targetSet,
                                           ManagedGroup * subRegion,
                                           string const & )
  {
    // 2.0. Check pressure and record composition bc application
    localIndex const comp = bc->GetComponent();
    GEOS_ERROR_IF( bcStatusMap.count( setName ) == 0, "Pressure boundary condition not prescribed on set '" << setName << "'" );
    GEOS_ERROR_IF( bcStatusMap[setName][comp], "Conflicting composition[" << comp << "] boundary conditions on set '" << setName << "'" );
    bcStatusMap[setName][comp] = true;

    // 2.1. Apply BC to set the field values
    bc->ApplyBoundaryConditionToField<BcEqual>( targetSet,
                                                time_n + dt,
                                                subRegion,
                                                viewKeyStruct::globalCompFractionString );
  });

  // 2.3 Check consistency between composition BC applied to sets
  bool bcConsistent = true;
  for (auto const & bcEntry : bcStatusMap)
  {
    for (localIndex ic = 0; ic < m_numComponents; ++ic)
    {
      bcConsistent &= bcEntry.second[ic];
      GEOS_WARNING_IF( !bcConsistent, "Composition boundary condition not applied to component "
                                      << ic << " on set '" << bcEntry.first << "'" );
    }
  }
  GEOS_ERROR_IF( !bcConsistent, "Inconsistent composition boundary conditions" );

  // 3. Call constitutive update, back-calculate target global component densities and apply to the system
  bcManager->ApplyBoundaryCondition( time_n + dt,
                                     domain,
                                     "ElementRegions",
                                     viewKeyStruct::pressureString,
                                     [&] ( BoundaryConditionBase const * const bc,
                                           string const & setName,
                                           set<localIndex> const & targetSet,
                                           ManagedGroup * subRegion,
                                           string const & )
  {
    MultiFluidBase * const fluid = GetFluidModel(subRegion );

    arrayView1d<globalIndex const> const & dofNumber = subRegion->getReference<array1d<globalIndex>>( viewKeyStruct::blockLocalDofNumberString );

    arrayView1d<real64 const> const & pres      = subRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );
    arrayView1d<real64 const> const & dPres     = subRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );
    arrayView1d<real64 const> const & bcPres    = subRegion->getReference<array1d<real64>>( viewKeyStruct::bcPressureString );
    arrayView2d<real64 const> const & compFrac  = subRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompFractionString );
    arrayView2d<real64 const> const & compDens  = subRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompDensityString );
    arrayView2d<real64 const> const & dCompDens = subRegion->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

    arrayView2d<real64> const & totalDens = fluid->getReference<array2d<real64>>( MultiFluidBase::viewKeyStruct::totalDensityString );

    Epetra_FEVector * const rhs = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );

    array1d<real64> rhsContribution( targetSet.size() * m_numDofPerCell );
    array1d<globalIndex> dof( targetSet.size() * m_numDofPerCell );

    integer counter = 0;

    for (localIndex a : targetSet)
    {
      fluid->StateUpdatePointMultiFluid(bcPres[a], m_temperature, compFrac[a], a, 0);

      globalIndex const offset = m_numDofPerCell * dofNumber[a];
      dof[counter] = offset;

      // 4.1. Apply pressure to the matrix
      BcEqual::ApplyBcValue( dof[counter],
                             blockSystem,
                             BlockIDs::compositionalBlock,
                             rhsContribution[counter],
                             bcPres[a],
                             pres[a] + dPres[a] );

      ++counter;

      // 4.2. For each component, apply target global density value
      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        dof[counter] = offset + ic + 1;
        real64 const targetCompDens = totalDens[a][0] * compFrac[a][ic];

        BcEqual::ApplyBcValue( dof[counter],
                               blockSystem,
                               BlockIDs::compositionalBlock,
                               rhsContribution[counter],
                               targetCompDens,
                               compDens[a][ic] + dCompDens[a][ic] );

        ++counter;
      }

      // 4.3. Apply accumulated rhs values
      BcEqual::ReplaceGlobalValues( rhs,
                                    counter,
                                    dof.data(),
                                    rhsContribution.data() );
    }
  });

}

void
CompositionalMultiphaseFlow::ApplyFaceDirichletBC_implicit( DomainPartition * const domain,
                                                            real64 const time_n, real64 const dt,
                                                            EpetraBlockSystem * const blockSystem )
{

}

real64
CompositionalMultiphaseFlow::CalculateResidualNorm( EpetraBlockSystem const * const blockSystem,
                                                    DomainPartition * const domain )
{
  Epetra_FEVector const * const residual = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );
  Epetra_Map      const * const rowMap   = blockSystem->GetRowMap( BlockIDs::compositionalBlock );

  // get a view into local residual vector
  int localSizeInt;
  double* localResidual = nullptr;
  residual->ExtractView(&localResidual, &localSizeInt);

  real64 localResidualNorm = 0.0;
  applyToSubRegions( domain, [&] ( localIndex er, localIndex esr,
                                   ElementRegion const * const region,
                                   CellBlockSubRegion const * const subRegion )
  {
    arrayView1d<integer const>     const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<globalIndex const> const & dofNumber     = m_dofNumber[er][esr];
    arrayView1d<real64 const>      const & refPoro       = m_porosityRef[er][esr];
    arrayView1d<real64 const>      const & volume        = m_volume[er][esr];

    localResidualNorm += sum_in_range( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] < 0)
      {
        real64 cell_norm = 0.0;
        globalIndex const offset = m_numDofPerCell * dofNumber[ei];
        for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
        {
          int const lid = rowMap->LID(integer_conversion<int>(offset + idof));
          real64 const val = localResidual[lid] / (refPoro[ei] * volume[ei]);
          cell_norm += val * val;
        }
        return cell_norm;
      }
      return 0.0;
    } );
  } );

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
    GEOS_LOG_RANK("\nSolution:\n" << *solution);
  }
}

bool
CompositionalMultiphaseFlow::CheckSystemSolution( EpetraBlockSystem const * const blockSystem,
                                                  real64 const scalingFactor,
                                                  DomainPartition * const domain )
{
  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( BlockIDs::compositionalBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( BlockIDs::compositionalBlock );

  int dummy;
  double* local_solution = nullptr;
  solution->ExtractView(&local_solution,&dummy);

  RAJA::ReduceMin<reducePolicy, integer> result(1);

  applyToSubRegions( domain, [&] ( localIndex er, localIndex esr,
                                   ElementRegion * const region,
                                   CellBlockSubRegion * const subRegion )
  {
    arrayView1d<integer     const> const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<globalIndex const> const & dofNumber     = m_dofNumber[er][esr];

    arrayView1d<real64 const> const & pres      = m_pressure[er][esr];
    arrayView1d<real64 const> const & dPres     = m_deltaPressure[er][esr];
    arrayView2d<real64 const> const & compDens  = m_globalCompDensity[er][esr];
    arrayView2d<real64 const> const & dCompDens = m_deltaGlobalCompDensity[er][esr];

    for_elems_in_subRegion( subRegion, GEOSX_LAMBDA ( localIndex const ei )
    {
      if (elemGhostRank[ei] >= 0)
        return;

      globalIndex const offset = m_numDofPerCell * dofNumber[ei];
      // extract solution and apply to dP
      {
        int const lid = rowMap->LID(integer_conversion<int>(offset));
        real64 const newPres = pres[ei] + dPres[ei] + scalingFactor * local_solution[lid];
        if (newPres < 0.0)
        {
          result.min(0);
        }
      }

      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        int const lid = rowMap->LID(integer_conversion<int>(offset + ic + 1));
        real64 const newDens = compDens[ei][ic] + dCompDens[ei][ic] + scalingFactor * local_solution[lid];
        if (newDens < 0.0)
        {
          result.min(0);
        }
      }
    });
  });

  return result.get() > 0;
}

void
CompositionalMultiphaseFlow::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                                  real64 const scalingFactor,
                                                  DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( BlockIDs::compositionalBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( BlockIDs::compositionalBlock );

  int dummy;
  double * local_solution = nullptr;
  solution->ExtractView( &local_solution, &dummy );

  applyToSubRegions( domain, [&] ( localIndex er, localIndex esr,
                                   ElementRegion * const region,
                                   CellBlockSubRegion * const subRegion )
  {
    arrayView1d<integer     const> const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<globalIndex const> const & dofNumber     = m_dofNumber[er][esr];

    arrayView1d<real64> const & dPres     = m_deltaPressure[er][esr];
    arrayView2d<real64> const & dCompDens = m_deltaGlobalCompDensity[er][esr];

    for_elems_in_subRegion( subRegion, GEOSX_LAMBDA ( localIndex const ei )
    {
      if (elemGhostRank[ei] >= 0)
        return;

      globalIndex const offset = m_numDofPerCell * dofNumber[ei];
      // extract solution and apply to dP
      {
        int const lid = rowMap->LID( integer_conversion<int>( offset ) );
        dPres[ei] += scalingFactor * local_solution[lid];
      }

      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        int const lid = rowMap->LID( integer_conversion<int>( offset + ic + 1 ) );
        dCompDens[ei][ic] += scalingFactor * local_solution[lid];
      }
    });
  });

  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::deltaPressureString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaGlobalCompDensityString );
  CommunicationTools::SynchronizeFields( fieldNames,
                                         mesh,
                                         domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  UpdateStateAll(domain);
}

void CompositionalMultiphaseFlow::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  applyToSubRegions( domain, [&] ( localIndex er, localIndex esr,
                                   ElementRegion * elementRegion,
                                   CellBlockSubRegion * const subRegion )
  {
    arrayView1d<real64> const & dPres     = m_deltaPressure[er][esr];
    arrayView2d<real64> const & dCompDens = m_deltaGlobalCompDensity[er][esr];

    for_elems_in_subRegion( subRegion, GEOSX_LAMBDA ( localIndex const ei )
    {
      dPres[ei] = 0.0;
      for (localIndex ic = 0; ic < m_numComponents; ++ic)
        dCompDens[ei][ic] = 0.0;
    });
  });

  UpdateStateAll(domain);
}

void CompositionalMultiphaseFlow::ImplicitStepComplete( real64 const & time,
                                                        real64 const & dt,
                                                        DomainPartition * const domain )
{
  applyToSubRegions( domain, [&] ( localIndex er, localIndex esr,
                                   ElementRegion * elementRegion,
                                   CellBlockSubRegion * const subRegion )
  {
    arrayView1d<real64 const> const & dPres     = m_deltaPressure[er][esr];
    arrayView2d<real64 const> const & dCompDens = m_deltaGlobalCompDensity[er][esr];

    arrayView1d<real64> const & pres     = m_pressure[er][esr];
    arrayView2d<real64> const & compDens = m_globalCompDensity[er][esr];

    for_elems_in_subRegion( subRegion, GEOSX_LAMBDA ( localIndex const ei )
    {
      pres[ei] += dPres[ei];
      for (localIndex ic = 0; ic < m_numComponents; ++ic)
        compDens[ei][ic] += dCompDens[ei][ic];
    });
  });
}

void CompositionalMultiphaseFlow::ResetViews(DomainPartition * const domain)
{
  FlowSolverBase::ResetViews( domain );

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  m_dofNumber =
    elemManager->ConstructViewAccessor<array1d<globalIndex>, arrayView1d<globalIndex>>( viewKeyStruct::blockLocalDofNumberString );
  m_pressure =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::pressureString );
  m_deltaPressure =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::deltaPressureString );
  m_globalCompDensity =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::globalCompDensityString );
  m_deltaGlobalCompDensity =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

  m_compFrac =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::globalCompFractionString );
  m_dCompFrac_dCompDens =
    elemManager->ConstructViewAccessor<array3d<real64>, arrayView3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );
  m_phaseVolFrac =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::phaseVolumeFractionString );
  m_dPhaseVolFrac_dPres =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dPressureString );
  m_dPhaseVolFrac_dCompDens =
    elemManager->ConstructViewAccessor<array3d<real64>, arrayView3d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  m_porosityOld =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::porosityOldString );
  m_phaseVolFracOld =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::phaseVolumeFractionOldString );
  m_phaseDensOld =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::phaseDensityOldString );
  m_phaseCompFracOld =
    elemManager->ConstructViewAccessor<array3d<real64>, arrayView3d<real64>>( viewKeyStruct::phaseComponentFractionOldString );

  m_pvMult =
    elemManager->ConstructMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString,
                                                                                      constitutiveManager );
  m_dPvMult_dPres =
    elemManager->ConstructMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( ConstitutiveBase::viewKeyStruct::dPVMult_dPresString,
                                                                                      constitutiveManager );
  m_phaseFrac =
    elemManager->ConstructMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::phaseFractionString,
                                                                                      constitutiveManager );
  m_dPhaseFrac_dPres =
    elemManager->ConstructMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dPressureString,
                                                                                      constitutiveManager );
  m_dPhaseFrac_dComp =
    elemManager->ConstructMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dGlobalCompFractionString,
                                                                                      constitutiveManager );
  m_phaseDens =
    elemManager->ConstructMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::phaseDensityString,
                                                                                      constitutiveManager );
  m_dPhaseDens_dPres =
    elemManager->ConstructMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString,
                                                                                      constitutiveManager );
  m_dPhaseDens_dComp =
    elemManager->ConstructMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dGlobalCompFractionString,
                                                                                      constitutiveManager );
  m_phaseVisc =
    elemManager->ConstructMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::phaseViscosityString,
                                                                                      constitutiveManager );
  m_dPhaseVisc_dPres =
    elemManager->ConstructMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dPressureString,
                                                                                      constitutiveManager );
  m_dPhaseVisc_dComp =
    elemManager->ConstructMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dGlobalCompFractionString,
                                                                                      constitutiveManager );
  m_phaseCompFrac =
    elemManager->ConstructMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::phaseCompFractionString,
                                                                                      constitutiveManager );
  m_dPhaseCompFrac_dPres =
    elemManager->ConstructMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dPressureString,
                                                                                      constitutiveManager );
  m_dPhaseCompFrac_dComp =
    elemManager->ConstructMaterialViewAccessor<array5d<real64>, arrayView5d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dGlobalCompFractionString,
                                                                                      constitutiveManager );
  m_totalDens =
    elemManager->ConstructMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( MultiFluidBase::viewKeyStruct::totalDensityString,
                                                                                      constitutiveManager );
  m_phaseRelPerm =
    elemManager->ConstructMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( RelativePermeabilityBase::viewKeyStruct::phaseRelPermString,
                                                                                      constitutiveManager );
  m_dPhaseRelPerm_dPhaseVolFrac =
    elemManager->ConstructMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( RelativePermeabilityBase::viewKeyStruct::dPhaseRelPerm_dPhaseVolFractionString,
                                                                                      constitutiveManager );
}


REGISTER_CATALOG_ENTRY(SolverBase, CompositionalMultiphaseFlow, string const &, ManagedGroup * const)
}// namespace geosx
