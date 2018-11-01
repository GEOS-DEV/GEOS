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
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/Fluid/MultiFluidBase.hpp"
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

CompositionalMultiphaseFlow::CompositionalMultiphaseFlow( const string & name,
                                                          ManagedGroup * const parent )
  :
  FlowSolverBase( name, parent ),
  m_numPhases( 0 ),
  m_numComponents( 0 ),
  m_numDofPerCell( 0 )
{
  // set the blockID for the block system interface
  getLinearSystemRepository()->SetBlockID(BlockIDs::compositionalBlock, this->getName());

  this->RegisterViewWrapper( viewKeysCompMultiphaseFlow.temperature.Key(), &m_temperature, false );
  this->RegisterViewWrapper( viewKeysCompMultiphaseFlow.useMassFlag.Key(), &m_useMass, false );
}

void CompositionalMultiphaseFlow::FillDocumentationNode()
{
  FlowSolverBase::FillDocumentationNode();

  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName(CompositionalMultiphaseFlow::CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("A compositional multiphase flow solver");

  docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.temperature.Key(),
                              viewKeysCompMultiphaseFlow.temperature.Key(),
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

  docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.useMassFlag.Key(),
                              viewKeysCompMultiphaseFlow.useMassFlag.Key(),
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
}

void CompositionalMultiphaseFlow::FillOtherDocumentationNodes( ManagedGroup * const rootGroup )
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

      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.pressure.Key(),
                                  viewKeysCompMultiphaseFlow.pressure.Key(),
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

      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.deltaPressure.Key(),
                                  viewKeysCompMultiphaseFlow.deltaPressure.Key(),
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
      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.bcPressure.Key(),
                                  viewKeysCompMultiphaseFlow.bcPressure.Key(),
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

      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.globalCompDensity.Key(),
                                  viewKeysCompMultiphaseFlow.globalCompDensity.Key(),
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

      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.deltaGlobalCompDensity.Key(),
                                  viewKeysCompMultiphaseFlow.deltaGlobalCompDensity.Key(),
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

      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.globalCompFraction.Key(),
                                  viewKeysCompMultiphaseFlow.globalCompFraction.Key(),
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

      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.dGlobalCompFraction_dGlobalCompDensity.Key(),
                                  viewKeysCompMultiphaseFlow.dGlobalCompFraction_dGlobalCompDensity.Key(),
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

      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.phaseVolumeFraction.Key(),
                                  viewKeysCompMultiphaseFlow.phaseVolumeFraction.Key(),
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

      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.dPhaseVolumeFraction_dPressure.Key(),
                                  viewKeysCompMultiphaseFlow.dPhaseVolumeFraction_dPressure.Key(),
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

      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.dPhaseVolumeFraction_dGlobalCompDensity.Key(),
                                  viewKeysCompMultiphaseFlow.dPhaseVolumeFraction_dGlobalCompDensity.Key(),
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

      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.phaseVolumeFractionOld.Key(),
                                  viewKeysCompMultiphaseFlow.phaseVolumeFractionOld.Key(),
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

      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.phaseDensityOld.Key(),
                                  viewKeysCompMultiphaseFlow.phaseDensityOld.Key(),
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

      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.phaseComponentFractionOld.Key(),
                                  viewKeysCompMultiphaseFlow.phaseComponentFractionOld.Key(),
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

      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.porosityOld.Key(),
                                  viewKeysCompMultiphaseFlow.porosityOld.Key(),
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

      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.blockLocalDofNumber.Key(),
                                  viewKeysCompMultiphaseFlow.blockLocalDofNumber.Key(),
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

      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.facePressure.Key(),
                                  viewKeysCompMultiphaseFlow.facePressure.Key(),
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

      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.globalCompFraction.Key(),
                                  viewKeysCompMultiphaseFlow.globalCompFraction.Key(),
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

      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.phaseDensityOld.Key(),
                                  viewKeysCompMultiphaseFlow.phaseDensityOld.Key(),
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

      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.phaseComponentFractionOld.Key(),
                                  viewKeysCompMultiphaseFlow.phaseComponentFractionOld.Key(),
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

      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.phaseViscosity.Key(),
                                  viewKeysCompMultiphaseFlow.phaseViscosity.Key(),
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

      docNode->AllocateChildNode( viewKeysCompMultiphaseFlow.phaseRelativePermeability.Key(),
                                  viewKeysCompMultiphaseFlow.phaseRelativePermeability.Key(),
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

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>( keys::domain );

  ConstitutiveManager * const cm = domain->getConstitutiveManager();
  ConstitutiveBase const * fluid = cm->GetConstitituveRelation( this->m_fluidName );
  MultiFluidBase const * mpFluid = fluid->group_cast<MultiFluidBase const *>();

  m_numPhases = mpFluid->numFluidPhases();
  m_numComponents = mpFluid->numFluidComponents();

  // compute number of DOF per cell
  m_numDofPerCell = m_numComponents + 1;
}

void CompositionalMultiphaseFlow::ResizeFields( DomainPartition * domain )
{
  for (auto & mesh : domain->getMeshBodies()->GetSubGroups())
  {
    MeshLevel * meshLevel = ManagedGroup::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);
    ElementRegionManager * const elemManager = meshLevel->getElemManager();

    elemManager->forCellBlocks([&](CellBlockSubRegion * const cellBlock) -> void
    {
      cellBlock->getReference<array2d<real64>>(viewKeysCompMultiphaseFlow.globalCompDensity).resizeDimension<1>(m_numComponents);
      cellBlock->getReference<array2d<real64>>(viewKeysCompMultiphaseFlow.deltaGlobalCompDensity).resizeDimension<1>(m_numComponents);

      cellBlock->getReference<array2d<real64>>(viewKeysCompMultiphaseFlow.globalCompFraction).resizeDimension<1>(m_numComponents);
      cellBlock->getReference<array3d<real64>>(viewKeysCompMultiphaseFlow.dGlobalCompFraction_dGlobalCompDensity).resizeDimension<1,2>(m_numComponents, m_numComponents);

      cellBlock->getReference<array2d<real64>>(viewKeysCompMultiphaseFlow.phaseVolumeFraction).resizeDimension<1>(m_numPhases);
      cellBlock->getReference<array2d<real64>>(viewKeysCompMultiphaseFlow.dPhaseVolumeFraction_dPressure).resizeDimension<1>(m_numPhases);
      cellBlock->getReference<array3d<real64>>(viewKeysCompMultiphaseFlow.dPhaseVolumeFraction_dGlobalCompDensity).resizeDimension<1,2>(m_numPhases, m_numComponents);

      cellBlock->getReference<array2d<real64>>(viewKeysCompMultiphaseFlow.phaseVolumeFractionOld).resizeDimension<1>(m_numPhases);
      cellBlock->getReference<array2d<real64>>(viewKeysCompMultiphaseFlow.phaseDensityOld).resizeDimension<1>(m_numPhases);
      cellBlock->getReference<array3d<real64>>(viewKeysCompMultiphaseFlow.phaseComponentFractionOld).resizeDimension<1,2>(m_numPhases, m_numComponents);
    });

    {
      FaceManager * const faceManager = meshLevel->getFaceManager();

      faceManager->getReference<array2d<real64>>(viewKeysCompMultiphaseFlow.phaseDensityOld.Key()).resizeDimension<1>(m_numPhases);
      faceManager->getReference<array2d<real64>>(viewKeysCompMultiphaseFlow.phaseViscosity.Key()).resizeDimension<1>(m_numPhases);
      faceManager->getReference<array2d<real64>>(viewKeysCompMultiphaseFlow.phaseRelativePermeability.Key()).resizeDimension<1>(m_numPhases);
      faceManager->getReference<array3d<real64>>(viewKeysCompMultiphaseFlow.phaseComponentFractionOld.Key()).resizeDimension<1,2>(m_numPhases, m_numComponents);;
    }
  }
}

void CompositionalMultiphaseFlow::IntermediateInitializationPreSubGroups( ManagedGroup * const rootGroup )
{
  FlowSolverBase::IntermediateInitializationPreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  ResizeFields( domain );
}

void CompositionalMultiphaseFlow::UpdateComponentFraction( DomainPartition * domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  auto const compDens = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.globalCompDensity.Key() );
  auto const dCompDens = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.deltaGlobalCompDensity.Key() );

  auto compFrac =
    elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.globalCompFraction.Key() );

  auto dCompFrac_dCompDens =
    elemManager->ConstructViewAccessor<array3d<real64>>( viewKeysCompMultiphaseFlow.dGlobalCompFraction_dGlobalCompDensity.Key() );

  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei ) -> void
  {
    real64 totalDensity = 0.0;
    for (localIndex ic = 0; ic < m_numComponents; ++ic)
    {
      totalDensity += compDens[er][esr][ei][ic] + dCompDens[er][esr][ei][ic];
      dCompFrac_dCompDens[er][esr][ei][ic][ic] = 1.0;
    }

    real64 const totalDensityInv = 1.0 / totalDensity;

    for (localIndex ic = 0; ic < m_numComponents; ++ic)
    {
      compFrac[er][esr][ei][ic] = (compDens[er][esr][ei][ic] + dCompDens[er][esr][ei][ic]) * totalDensityInv;
      for (localIndex jc = 0; jc < m_numComponents; ++jc)
      {
        dCompFrac_dCompDens[er][esr][ei][ic][jc] = - compFrac[er][esr][ei][ic];
        dCompFrac_dCompDens[er][esr][ei][ic][jc] *= totalDensityInv;
      }
    }
  });
}

void CompositionalMultiphaseFlow::UpdatePhaseVolumeFraction( DomainPartition * domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  auto const compDens = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.globalCompDensity.Key() );
  auto const dCompDens = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.deltaGlobalCompDensity.Key() );

  auto phaseVolFrac = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.phaseVolumeFraction.Key() );
  auto dPhaseVolFrac_dPres = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.dPhaseVolumeFraction_dPressure.Key() );
  auto dPhaseVolFrac_dComp = elemManager->ConstructViewAccessor<array3d<real64>>( viewKeysCompMultiphaseFlow.dPhaseVolumeFraction_dGlobalCompDensity.Key() );

  auto phaseFrac = elemManager->ConstructMaterialViewAccessor<array3d<real64>>( MultiFluidBase::
                                                                                viewKeyStruct::
                                                                                phaseFractionString,
                                                                                constitutiveManager );

  auto dPhaseFrac_dPres = elemManager->ConstructMaterialViewAccessor<array3d<real64>>( MultiFluidBase::
                                                                                       viewKeyStruct::
                                                                                       dPhaseFraction_dPressureString,
                                                                                       constitutiveManager );

  auto dPhaseFrac_dComp = elemManager->ConstructMaterialViewAccessor<array4d<real64>>( MultiFluidBase::
                                                                                       viewKeyStruct::
                                                                                       dPhaseFraction_dGlobalCompFractionString,
                                                                                       constitutiveManager );

  auto phaseDens = elemManager->ConstructMaterialViewAccessor<array3d<real64>>( MultiFluidBase::
                                                                                viewKeyStruct::
                                                                                phaseDensityString,
                                                                                constitutiveManager );

  auto dPhaseDens_dPres = elemManager->ConstructMaterialViewAccessor<array3d<real64>>( MultiFluidBase::
                                                                                       viewKeyStruct::
                                                                                       dPhaseDensity_dPressureString,
                                                                                       constitutiveManager );

  auto dPhaseDens_dComp = elemManager->ConstructMaterialViewAccessor<array4d<real64>>( MultiFluidBase::
                                                                                       viewKeyStruct::
                                                                                       dPhaseDensity_dGlobalCompFractionString,
                                                                                       constitutiveManager );

  auto dCompFrac_dCompDens =
    elemManager->ConstructViewAccessor<array3d<real64>>( viewKeysCompMultiphaseFlow.dGlobalCompFraction_dGlobalCompDensity.Key() );

  localIndex const NC = m_numComponents;

  array1d<real64> work( NC );

//  forAllElemsInMesh( mesh, [&]( localIndex const er,
//                                localIndex const esr,
//                                localIndex const ei ) -> void
//  {

  for (localIndex er = 0; er < elemManager->numRegions(); ++er)
  {
    ElementRegion const * const elemRegion = elemManager->GetRegion(er);
    for (localIndex esr = 0; esr < elemRegion->numSubRegions(); ++esr)
    {
      CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion(esr);

      arrayView2d<real64 const> compDensSub = compDens[er][esr].get();
      arrayView2d<real64 const> dCompDensSub = dCompDens[er][esr].get();

      arrayView3d<real64 const> phaseFracSub = phaseFrac[er][esr][m_fluidIndex].get();
      arrayView3d<real64 const> dPhaseFrac_dPresSub = dPhaseFrac_dPres[er][esr][m_fluidIndex].get();
      arrayView4d<real64 const> dPhaseFrac_dCompSub = dPhaseFrac_dComp[er][esr][m_fluidIndex].get();

      arrayView3d<real64 const> phaseDensSub = phaseDens[er][esr][m_fluidIndex].get();
      arrayView3d<real64 const> dPhaseDens_dPresSub = dPhaseDens_dPres[er][esr][m_fluidIndex].get();
      arrayView4d<real64 const> dPhaseDens_dCompSub = dPhaseDens_dComp[er][esr][m_fluidIndex].get();

      arrayView2d<real64> phaseVolFracSub = phaseVolFrac[er][esr].get();
      arrayView2d<real64> dPhaseVolFrac_dPresSub = dPhaseVolFrac_dPres[er][esr].get();
      arrayView3d<real64> dPhaseVolFrac_dCompSub = dPhaseVolFrac_dComp[er][esr].get();

      arrayView3d<real64 const> dCompFrac_dCompDensSub = dCompFrac_dCompDens[er][esr].get();

      for (localIndex ei = 0; ei < subRegion->size(); ++ei)
      {

        // compute total density from component partial densities
        real64 totalDensity = 0.0;
        real64 const dTotalDens_dCompDens = 1.0;
        for (localIndex ic = 0; ic < m_numComponents; ++ic)
        {
          totalDensity += compDensSub[ei][ic] + dCompDensSub[ei][ic];
        }

        for (localIndex ip = 0; ip < m_numPhases; ++ip)
        {
          // Expression for volume fractions: S_p = (nu_p / rho_p) * rho_t
          real64 const phaseDensInv = 1.0 / phaseDensSub[ei][0][ip];

          // compute saturation and derivatives except multiplying by the total density
          phaseVolFracSub[ei][ip] = phaseFracSub[ei][0][ip] * phaseDensInv;

          dPhaseVolFrac_dPresSub[ei][ip] =
            (dPhaseFrac_dPresSub[ei][0][ip] - phaseVolFracSub[ei][ip] * dPhaseDens_dPresSub[ei][0][ip]) * phaseDensInv;

          for (localIndex jc = 0; jc < m_numComponents; ++jc)
          {
            dPhaseVolFrac_dCompSub[ei][ip][jc] =
              (dPhaseFrac_dCompSub[ei][0][ip][jc] - phaseVolFracSub[ei][ip] * dPhaseDens_dCompSub[ei][0][ip][jc]) * phaseDensInv;
          }

          // apply chain rule to convert derivatives from global component fractions to densities
          applyChainRuleInPlace( NC, dCompFrac_dCompDensSub[ei], dPhaseVolFrac_dCompSub[ei][ip], work );

          // now finalize the computation by multiplying by total density
          for (localIndex jc = 0; jc < m_numComponents; ++jc)
          {
            dPhaseVolFrac_dCompSub[ei][ip][jc] *= totalDensity;
            dPhaseVolFrac_dCompSub[ei][ip][jc] += phaseVolFracSub[ei][ip] * dTotalDens_dCompDens;
          }

          phaseVolFracSub[ei][ip] *= totalDensity;
          dPhaseVolFrac_dPresSub[ei][ip] *= totalDensity;
        }
      }
    }
  }//);
}

void CompositionalMultiphaseFlow::UpdateFluidModel( ManagedGroup * dataGroup,
                                                    set<localIndex> const & targetSet,
                                                    string const & pressureFieldName,
                                                    string const & deltaPressureFieldName,
                                                    string const & compFracFieldName )
{
  // TODO for now can only be a subregion; this will change in the future
  CellBlockSubRegion * subRegion = dataGroup->group_cast<CellBlockSubRegion *>();

  GEOS_ERROR_IF( subRegion == nullptr, "Target data group is not a subregion" );

  ManagedGroup * const constitutiveGroup = subRegion->GetConstitutiveModels();
  MultiFluidBase * const fluid = constitutiveGroup->GetGroup<MultiFluidBase>( m_fluidIndex ); // TODO could be incorrect

  GEOS_ERROR_IF( fluid == nullptr, "Fluid model does not exist in subregion " << subRegion->getName() );

  // TODO get rid of this branching if/when formulation changed from P_{n}+dP to P_{n+1}
  if ( !deltaPressureFieldName.empty() )
  {
    auto const & pres = dataGroup->getReference<array1d<real64>>(pressureFieldName);
    auto const & dPres = dataGroup->getReference<array1d<real64>>(deltaPressureFieldName);
    auto const & compFrac = dataGroup->getReference<array2d<real64>>(compFracFieldName);

    // TODO make this a RAJA loop?
    for (localIndex a : targetSet)
    {
      fluid->StateUpdatePointMultiphaseFluid( pres[a] + dPres[a], m_temperature,
                                              static_cast<real64 const *>(compFrac[a]),
                                              a, 0 );
    }
  }
  else
  {
    auto const & pres = dataGroup->getReference<array1d<real64>>(pressureFieldName);
    auto const & compFrac = dataGroup->getReference<array2d<real64>>(compFracFieldName);

    // TODO make this a RAJA loop?
    for (localIndex a : targetSet)
    {
      fluid->StateUpdatePointMultiphaseFluid( pres[a], m_temperature,
                                              static_cast<real64 const *>(compFrac[a]),
                                              a, 0 );
    }
  }
}

void CompositionalMultiphaseFlow::UpdateFluidModel( ManagedGroup * dataGroup, string const & setName )
{
  ManagedGroup * sets = dataGroup->GetGroup( keys::sets );

  GEOS_ERROR_IF( sets == nullptr, "Target group does not contain sets" );

  set<localIndex> const & targetSet = sets->getReference<set<localIndex>>( setName );
  UpdateFluidModel( dataGroup, targetSet );
}

void CompositionalMultiphaseFlow::UpdateFluidModels( DomainPartition * domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  elemManager->forCellBlocksComplete( [&] ( localIndex er,
                                            localIndex esr,
                                            ElementRegion * region,
                                            CellBlockSubRegion * subregion ) -> void
  {
    UpdateFluidModel( subregion, "all" );
  });
}

void CompositionalMultiphaseFlow::UpdateSolidModel( ManagedGroup * dataGroup, set<localIndex> const & targetSet )
{
  // TODO for now can only be a subregion; this will change in the future
  CellBlockSubRegion * subRegion = dataGroup->group_cast<CellBlockSubRegion *>();

  GEOS_ERROR_IF( subRegion == nullptr, "Target data group is not a subregion" );

  ManagedGroup * const constitutiveGroup = subRegion->GetConstitutiveModels();
  ConstitutiveBase * const solid = constitutiveGroup->GetGroup<ConstitutiveBase>( m_solidIndex ); // TODO could be incorrect

  GEOS_ERROR_IF( solid == nullptr, "Solid model does not exist in subregion " << subRegion->getName() );

  auto const & pres     = dataGroup->getReference<array1d<real64>>( viewKeysCompMultiphaseFlow.pressure.Key() );
  auto const & dPres    = dataGroup->getReference<array1d<real64>>( viewKeysCompMultiphaseFlow.deltaPressure.Key() );

  // TODO make this a RAJA loop?
  for ( localIndex a : targetSet)
  {
    solid->StateUpdatePointPressure( pres[a] + dPres[a], a, 0 );
  }
}

void CompositionalMultiphaseFlow::UpdateSolidModel( ManagedGroup * dataGroup, string const & setName )
{
  ManagedGroup * sets = dataGroup->GetGroup( keys::sets );

  GEOS_ERROR_IF( sets == nullptr, "Target group does not contain sets" );

  set<localIndex> const & targetSet = sets->getReference<set<localIndex>>( setName );
  UpdateSolidModel( dataGroup, targetSet );
}

void CompositionalMultiphaseFlow::UpdateSolidModels( DomainPartition * domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  elemManager->forCellBlocksComplete( [&] ( localIndex er,
                                            localIndex esr,
                                            ElementRegion * region,
                                            CellBlockSubRegion * subregion ) -> void
  {
    UpdateSolidModel( subregion, "all" );
  });
}

void CompositionalMultiphaseFlow::UpdateConstitutiveModels( DomainPartition * domain )
{
  UpdateFluidModels( domain );
  UpdateSolidModels( domain );
}

void CompositionalMultiphaseFlow::InitializeFluidState( DomainPartition * domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  auto compDens     = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.globalCompDensity.Key() );
  auto compFrac = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.globalCompFraction.Key() );

  auto totalDens    = elemManager->ConstructMaterialViewAccessor<array2d<real64>>( MultiFluidBase::
                                                                                   viewKeyStruct::
                                                                                   totalDensityString,
                                                                                   constitutiveManager );

  // 1. Assume global component fractions have been prescribed.
  // Update constitutive state to get fluid density.
  UpdateConstitutiveModels( domain );

  // 2. Back-calculate global component densities from fractions and total fluid density
  // in order to initialize the primary solution variables
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei ) -> void
  {
    for (localIndex ic = 0; ic < m_numComponents; ++ic)
    {
      compDens[er][esr][ei][ic] = totalDens[er][esr][m_fluidIndex][ei][0] * compFrac[er][esr][ei][ic];
    }
  });

  // 3. Calculate phase saturations
  UpdatePhaseVolumeFraction( domain );
}

void CompositionalMultiphaseFlow::FinalInitializationPreSubGroups( ManagedGroup * const rootGroup )
{
  FlowSolverBase::FinalInitializationPreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  //TODO this is a hack until the sets are fixed to include ghosts!!
  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back(viewKeysCompMultiphaseFlow.pressure.Key());
  fieldNames["elems"].push_back(viewKeysCompMultiphaseFlow.globalCompDensity.Key());
  CommunicationTools::SynchronizeFields(fieldNames,
                                        domain->getMeshBody(0)->getMeshLevel(0),
                                        domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  // set mass fraction flag on main model
  {
    ConstitutiveBase * const fluid = constitutiveManager->GetConstitituveRelation( this->m_fluidName );
    MultiFluidBase * const mpFluid = fluid->group_cast<MultiFluidBase *>();
    mpFluid->setMassFlag( m_useMass ); // this formulation uses mass fractions
  }

  // set mass fraction flag on subregion models
  auto constitutiveRelations = elemManager->ConstructConstitutiveAccessor<ConstitutiveBase>( constitutiveManager );
  elemManager->forCellBlocksComplete( [&] ( localIndex er,
                                            localIndex esr,
                                            ElementRegion * elementRegion,
                                            CellBlockSubRegion * cellBlock ) -> void
  {
    ConstitutiveBase * fluid = constitutiveRelations[er][esr][m_fluidIndex];
    MultiFluidBase * mpFluid = fluid->group_cast<MultiFluidBase *>();
    mpFluid->setMassFlag( m_useMass ); // this formulation uses mass fractions
  });

  // Initialize primary variables from applied initial conditions
  InitializeFluidState( domain );
}

real64 CompositionalMultiphaseFlow::SolverStep( real64 const & time_n,
                                                real64 const & dt,
                                                integer const cycleNumber,
                                                DomainPartition * domain )
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

void CompositionalMultiphaseFlow::BackupFields( DomainPartition * domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  auto phaseDensOld         = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.phaseDensityOld.Key() );
  auto phaseVolFrac         = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.phaseVolumeFraction.Key() );
  auto phaseVolFracOld      = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.phaseVolumeFractionOld.Key() );
  auto phaseCompMassFracOld = elemManager->ConstructViewAccessor<array3d<real64>>( viewKeysCompMultiphaseFlow.phaseComponentFractionOld.Key() );
  auto poroOld              = elemManager->ConstructViewAccessor<array1d<real64>>( viewKeysCompMultiphaseFlow.porosityOld.Key() );
  auto poroRef              = elemManager->ConstructViewAccessor<array1d<real64>>( viewKeysCompMultiphaseFlow.referencePorosity.Key() );

  auto const pvmult =
    elemManager->ConstructMaterialViewAccessor< array2d<real64> >( ConstitutiveBase::
                                                                   viewKeyStruct::
                                                                   poreVolumeMultiplierString,
                                                                   constitutiveManager );


  auto const phaseDens =
    elemManager->ConstructMaterialViewAccessor< array3d<real64> >( MultiFluidBase::
                                                                   viewKeyStruct::
                                                                   phaseDensityString,
                                                                   constitutiveManager );

  auto const phaseCompMassFrac =
    elemManager->ConstructMaterialViewAccessor< array4d<real64> >( MultiFluidBase::
                                                                   viewKeyStruct::
                                                                   phaseCompFractionString,
                                                                   constitutiveManager );

// backup some fields used in time derivative approximation
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei )->void
  {
    for (localIndex ip = 0; ip < m_numPhases; ++ip)
    {
      phaseDensOld[er][esr][ei][ip] = phaseDens[er][esr][m_fluidIndex][ei][0][ip];
      phaseVolFracOld[er][esr][ei][ip] = phaseVolFrac[er][esr][ei][ip];
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
    elementRegionManager->ConstructViewAccessor<array1d<globalIndex>>( viewKeysCompMultiphaseFlow.blockLocalDofNumber.Key(), string() );

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

  GEOS_ERROR_IF( numLocalRows != numLocalRowsToSend, "number of local rows inconsistent" );

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
      blockLocalDofNumber[er][esr][ei] = firstLocalRow + localCount + offset;
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

  auto blockLocalDofNumber =
    elementRegionManager->ConstructViewAccessor<array1d<globalIndex>>( viewKeysCompMultiphaseFlow.blockLocalDofNumber.Key() );

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
      globalIndex const offset = m_numDofPerCell * blockLocalDofNumber[cell.region][cell.subRegion][cell.index];
      for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
      {
        elementLocalDofIndexRow[i * m_numDofPerCell + idof] = offset + idof;
      }
    });

    localIndex const stencilSize = stencil.size();
    elementLocalDofIndexCol.resize(stencilSize * m_numDofPerCell);
    stencil.forAll([&] (CellDescriptor const & cell, real64 w, localIndex const i) -> void
    {
      globalIndex const offset = m_numDofPerCell * blockLocalDofNumber[cell.region][cell.subRegion][cell.index];
      for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
      {
        elementLocalDofIndexCol[i * m_numDofPerCell + idof] = offset + idof;
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
      globalIndex const offset = m_numDofPerCell * blockLocalDofNumber[er][esr][ei];
      for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
      {
        elementLocalDofIndexRow[idof] = offset + idof;
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
          CellDescriptor const & cell = point.cellIndex;
          globalIndex const offset = m_numDofPerCell * blockLocalDofNumber[cell.region][cell.subRegion][cell.index];
          for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
          {
            elementLocalDofIndexRow[idof] = offset + idof;
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
          CellDescriptor const & cell = point.cellIndex;
          globalIndex const offset = m_numDofPerCell * blockLocalDofNumber[cell.region][cell.subRegion][cell.index];
          for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
          {
            elementLocalDofIndexCol[counter * m_numDofPerCell + idof] = offset + idof;
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

  localIndex_array displacementIndices;
  SetNumRowsAndTrilinosIndices( mesh,
                                numLocalRows,
                                numGlobalRows,
                                0 );

  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back(viewKeysCompMultiphaseFlow.blockLocalDofNumber.Key());
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

  AssembleAccumulationTerms( domain, blockSystem, time_n, dt );
  AssembleFluxTerms( domain, blockSystem, time_n, dt );
  AssembleVolumeBalanceTerms( domain, blockSystem, time_n, dt );

  jacobian->GlobalAssemble(false);
  residual->GlobalAssemble();

  if( verboseLevel() >= 3 )
  {
    jacobian->Print(std::cout);
    residual->Print(std::cout);
  }

}

void CompositionalMultiphaseFlow::AssembleAccumulationTerms( DomainPartition * const domain,
                                                             EpetraBlockSystem * const blockSystem,
                                                             real64 const time_n,
                                                             real64 const dt )
{
//***** extract data required for assembly of system *****
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  ConstitutiveManager const * constitutiveManager =
    domain->GetGroup<ConstitutiveManager>( keys::ConstitutiveManager );

  ElementRegionManager const * elemManager = mesh->getElemManager();

  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::compositionalBlock,
                                                                BlockIDs::compositionalBlock );
  Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );

  auto elemGhostRank =
    elemManager->ConstructViewAccessor<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

  auto blockLocalDofNumber =
    elemManager->ConstructViewAccessor<array1d<globalIndex >>( viewKeysCompMultiphaseFlow.blockLocalDofNumber.Key() );

  auto pres      = elemManager->ConstructViewAccessor<array1d<real64>>( viewKeysCompMultiphaseFlow.pressure.Key() );
  auto dPres     = elemManager->ConstructViewAccessor<array1d<real64>>( viewKeysCompMultiphaseFlow.deltaPressure.Key() );
  auto volume    = elemManager->ConstructViewAccessor<array1d<real64>>( CellBlock::viewKeyStruct::elementVolumeString );

  auto phaseVolFrac        = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.phaseVolumeFraction.Key() );
  auto dPhaseVolFrac_dPres = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.dPhaseVolumeFraction_dPressure.Key() );
  auto dPhaseVolFrac_dComp = elemManager->ConstructViewAccessor<array3d<real64>>( viewKeysCompMultiphaseFlow.dPhaseVolumeFraction_dGlobalCompDensity.Key() );

  auto phaseDensOld =
    elemManager->ConstructViewAccessor<array2d<real64>>(viewKeysCompMultiphaseFlow.phaseDensityOld.Key());

  auto phaseDens =
    elemManager->ConstructMaterialViewAccessor< array3d<real64> >( MultiFluidBase::
                                                                   viewKeyStruct::
                                                                   phaseDensityString,
                                                                   constitutiveManager );

  auto dPhaseDens_dPres =
    elemManager->ConstructMaterialViewAccessor< array3d<real64> >( MultiFluidBase::
                                                                   viewKeyStruct::
                                                                   dPhaseDensity_dPressureString,
                                                                   constitutiveManager );

  auto dPhaseDens_dComp =
    elemManager->ConstructMaterialViewAccessor< array4d<real64> >( MultiFluidBase::
                                                                   viewKeyStruct::
                                                                   dPhaseDensity_dGlobalCompFractionString,
                                                                   constitutiveManager );

  auto phaseVolFracOld =
    elemManager->ConstructViewAccessor<array2d<real64>>(viewKeysCompMultiphaseFlow.phaseVolumeFractionOld.Key() );

  auto phaseCompFracOld =
    elemManager->ConstructViewAccessor<array3d<real64>>(viewKeysCompMultiphaseFlow.phaseComponentFractionOld.Key());

  auto phaseCompFrac =
    elemManager->ConstructMaterialViewAccessor< array4d<real64> >( MultiFluidBase::
                                                                   viewKeyStruct::
                                                                   phaseCompFractionString,
                                                                   constitutiveManager );

  auto dPhaseCompFrac_dPres =
    elemManager->ConstructMaterialViewAccessor< array4d<real64> >( MultiFluidBase::
                                                                   viewKeyStruct::
                                                                   dPhaseCompFraction_dPressureString,
                                                                   constitutiveManager );

  auto dPhaseCompFrac_dComp =
    elemManager->ConstructMaterialViewAccessor< array5d<real64> >( MultiFluidBase::
                                                                   viewKeyStruct::
                                                                   dPhaseCompFraction_dGlobalCompFractionString,
                                                                   constitutiveManager );

  auto porosityOld =
    elemManager->ConstructViewAccessor<array1d<real64>>( viewKeysCompMultiphaseFlow.porosityOld.Key() );

  auto porosityRef =
    elemManager->ConstructViewAccessor<array1d<real64>>( viewKeysCompMultiphaseFlow.referencePorosity.Key() );

  auto pVMult =
    elemManager->ConstructMaterialViewAccessor< array2d<real64> >( ConstitutiveBase::
                                                                   viewKeyStruct::
                                                                   poreVolumeMultiplierString,
                                                                   constitutiveManager );

  auto dPVMult_dPres =
    elemManager->ConstructMaterialViewAccessor< array2d<real64> >( ConstitutiveBase::
                                                                   viewKeyStruct::
                                                                   dPVMult_dPresString,
                                                                   constitutiveManager );

  auto dCompFrac_dCompDens =
    elemManager->ConstructViewAccessor<array3d<real64>>( viewKeysCompMultiphaseFlow.dGlobalCompFraction_dGlobalCompDensity.Key() );

  localIndex const NC   = m_numComponents;
  localIndex const NP   = m_numPhases;
  localIndex const NDOF = m_numDofPerCell;

  // using Epetra types
  array1d<long long> localAccumDOF( NDOF );
  array1d<double> localAccum( NC );
  array2d<double> localAccumJacobian( NC, NDOF );

  //***** Loop over all elements and assemble the change in volume/density terms *****
//  forAllElemsInMesh( mesh, [=] ( localIndex const er,
//                                 localIndex const esr,
//                                 localIndex const ei ) -> void

  // temporary work arrays
  array1d<real64> dPhaseAmount_dC( NC );

  for (localIndex er = 0; er < elemManager->numRegions(); ++er)
  {
    ElementRegion const * const elemRegion = elemManager->GetRegion(er);
    for (localIndex esr = 0; esr < elemRegion->numSubRegions(); ++esr)
    {
      CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion(esr);

      // set up array views to reduce amount of indexing
      arrayView1d<globalIndex const> const dofNumber = blockLocalDofNumber[er][esr].get();
      arrayView2d<real64 const> const pVMultSub = pVMult[er][esr][m_solidIndex].get();
      arrayView2d<real64 const> const dPVMult_dPresSub = dPVMult_dPres[er][esr][m_solidIndex].get();

      arrayView2d<real64 const> const phaseDensOldSub = phaseDensOld[er][esr].get();
      arrayView3d<real64 const> const phaseDensSub = phaseDens[er][esr][m_fluidIndex].get();
      arrayView3d<real64 const> const dPhaseDens_dPresSub = dPhaseDens_dPres[er][esr][m_fluidIndex].get();
      arrayView4d<real64 const> const dPhaseDens_dCompSub = dPhaseDens_dComp[er][esr][m_fluidIndex].get();

      arrayView2d<real64 const> const phaseVolFracOldSub = phaseVolFracOld[er][esr].get();
      arrayView2d<real64 const> const phaseVolFracSub = phaseVolFrac[er][esr].get();
      arrayView2d<real64 const> const dPhaseVolFrac_dPresSub = dPhaseVolFrac_dPres[er][esr].get();
      arrayView3d<real64 const> const dPhaseVolFrac_dCompSub = dPhaseVolFrac_dComp[er][esr].get();

      arrayView3d<real64 const> const phaseCompFracOldSub = phaseCompFracOld[er][esr].get();
      arrayView4d<real64 const> const phaseCompFracSub = phaseCompFrac[er][esr][m_fluidIndex].get();
      arrayView4d<real64 const> const dPhaseCompFrac_dPresSub = dPhaseCompFrac_dPres[er][esr][m_fluidIndex].get();
      arrayView5d<real64 const> const dPhaseCompFrac_dCompSub = dPhaseCompFrac_dComp[er][esr][m_fluidIndex].get();

      arrayView3d<real64 const> const dCompFrac_dCompDensSub = dCompFrac_dCompDens[er][esr].get();

      for (localIndex ei = 0; ei < subRegion->size(); ++ei)
      {
        if (elemGhostRank[er][esr][ei] < 0)
        {
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
          real64 const volNew   = volume[er][esr][ei];
          real64 const volOld   = volume[er][esr][ei];
          real64 const dVol_dP  = 0.0; // used in poroelastic solver

          real64 const poroNew  = porosityRef[er][esr][ei] * pVMultSub[ei][0];
          real64 const poroOld  = porosityOld[er][esr][ei];
          real64 const dPoro_dP = porosityRef[er][esr][ei] * dPVMult_dPresSub[ei][0];

          real64 const poreVolNew = volNew * poroNew;
          real64 const poreVolOld = volOld * poroOld;
          real64 const dPoreVol_dP = dVol_dP * poroNew + volNew * dPoro_dP;

          // sum contributions to component accumulation from each phase
          for (localIndex ip = 0; ip < NP; ++ip)
          {
            real64 const phaseAmountNew = poreVolNew * phaseVolFracSub[ei][ip] * phaseDensSub[ei][0][ip];
            real64 const phaseAmountOld = poreVolOld * phaseVolFracOldSub[ei][ip] * phaseDensOldSub[ei][ip];

            real64 const dPhaseAmount_dP = dPoreVol_dP * phaseVolFracSub[ei][ip] * phaseDensSub[ei][0][ip]
                                          + poreVolNew * (dPhaseVolFrac_dPresSub[ei][ip] * phaseDensSub[ei][0][ip]
                                                       + phaseVolFracSub[ei][ip] * dPhaseDens_dPresSub[ei][0][ip]);

            // assemble density dependence
            applyChainRule( NC, dCompFrac_dCompDensSub[ei], dPhaseDens_dCompSub[ei][0][ip], dPhaseAmount_dC );
            for (localIndex jc = 0; jc < NC; ++jc)
            {
              dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc]     * phaseVolFracSub[ei][ip]
                                  + phaseDensSub[ei][0][ip] * dPhaseVolFrac_dCompSub[ei][ip][jc];
              dPhaseAmount_dC[jc] *= poreVolNew;
            }

            // ic - index of component whose conservation equation is assembled
            // (i.e. row number in local matrix)
            for (localIndex ic = 0; ic < NC; ++ic)
            {
              real64 const phaseCompAmountNew = phaseAmountNew * phaseCompFracSub[ei][0][ip][ic];
              real64 const phaseCompAmountOld = phaseAmountOld * phaseCompFracOldSub[ei][ip][ic];

              real64 const dPhaseCompAmount_dP = dPhaseAmount_dP * phaseCompFracSub[ei][0][ip][ic]
                                               + phaseAmountNew  * dPhaseCompFrac_dPresSub[ei][0][ip][ic];

              localAccum[ic] += phaseCompAmountNew - phaseCompAmountOld;
              localAccumJacobian[ic][0] = dPhaseCompAmount_dP;

              // jc - index of component w.r.t. whose compositional var the derivative is being taken
              // (i.e. col number in local matrix)

              // assemble phase composition dependence
              applyChainRule(NC, dCompFrac_dCompDensSub[ei], dPhaseCompFrac_dCompSub[ei][0][ip][ic], &localAccumJacobian[ic][1]);
              for (localIndex jc = 0; jc < NC; ++jc)
              {
                localAccumJacobian[ic][jc+1] = localAccumJacobian[ic][jc+1]    * phaseAmountNew
                                             + phaseCompFracSub[ei][0][ip][ic] * dPhaseAmount_dC[jc];
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
        }
      }
    }
  }//)
}

void CompositionalMultiphaseFlow::AssembleFluxTerms( DomainPartition * const domain,
                                                     EpetraBlockSystem * const blockSystem,
                                                     real64 const time_n,
                                                     real64 const dt )
{
  //***** extract data required for assembly of system *****
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  ConstitutiveManager const * constitutiveManager =
    domain->GetGroup<ConstitutiveManager>( keys::ConstitutiveManager );

  ElementRegionManager const * elemManager = mesh->getElemManager();

  NumericalMethodsManager const * numericalMethodManager =
    domain->getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager const * fvManager =
    numericalMethodManager->GetGroup<FiniteVolumeManager>(keys::finiteVolumeManager);

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation( m_discretizationName );
  FluxApproximationBase::CellStencil const & stencilCollection = fluxApprox->getStencil();

  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::compositionalBlock,
                                                                BlockIDs::compositionalBlock );
  Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );

  auto elemGhostRank =
    elemManager->ConstructViewAccessor<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

  auto blockLocalDofNumber =
    elemManager->ConstructViewAccessor<array1d<globalIndex >>( viewKeysCompMultiphaseFlow.blockLocalDofNumber.Key() );

  auto pres      = elemManager->ConstructViewAccessor<array1d<real64>>( viewKeysCompMultiphaseFlow.pressure.Key() );
  auto dPres     = elemManager->ConstructViewAccessor<array1d<real64>>( viewKeysCompMultiphaseFlow.deltaPressure.Key() );
  auto gravDepth = elemManager->ConstructViewAccessor<array1d<real64>>( viewKeysCompMultiphaseFlow.gravityDepth.Key() );

  auto phaseVolFrac        = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.phaseVolumeFraction.Key() );
  auto dPhaseVolFrac_dPres = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.dPhaseVolumeFraction_dPressure.Key() );
  auto dPhaseVolFrac_dComp = elemManager->ConstructViewAccessor<array3d<real64>>( viewKeysCompMultiphaseFlow.dPhaseVolumeFraction_dGlobalCompDensity.Key() );

  auto phaseDens =
    elemManager->ConstructMaterialViewAccessor< array3d<real64> >( MultiFluidBase::
                                                                   viewKeyStruct::
                                                                   phaseDensityString,
                                                                   constitutiveManager );

  auto dPhaseDens_dPres =
    elemManager->ConstructMaterialViewAccessor< array3d<real64> >( MultiFluidBase::
                                                                   viewKeyStruct::
                                                                   dPhaseDensity_dPressureString,
                                                                   constitutiveManager );

  auto dPhaseDens_dComp =
    elemManager->ConstructMaterialViewAccessor< array4d<real64> >( MultiFluidBase::
                                                                   viewKeyStruct::
                                                                   dPhaseDensity_dGlobalCompFractionString,
                                                                   constitutiveManager );

  auto phaseVisc =
    elemManager->ConstructMaterialViewAccessor< array3d<real64> >( MultiFluidBase::
                                                                   viewKeyStruct::
                                                                   phaseViscosityString,
                                                                   constitutiveManager );

  auto dPhaseVisc_dPres =
    elemManager->ConstructMaterialViewAccessor< array3d<real64> >( MultiFluidBase::
                                                                   viewKeyStruct::
                                                                   dPhaseViscosity_dPressureString,
                                                                   constitutiveManager );

  auto dPhaseVisc_dComp =
    elemManager->ConstructMaterialViewAccessor< array4d<real64> >( MultiFluidBase::
                                                                   viewKeyStruct::
                                                                   dPhaseViscosity_dGlobalCompFractionString,
                                                                   constitutiveManager );

  auto phaseCompFrac =
    elemManager->ConstructMaterialViewAccessor< array4d<real64> >( MultiFluidBase::
                                                                   viewKeyStruct::
                                                                   phaseCompFractionString,
                                                                   constitutiveManager );

  auto dPhaseCompFrac_dPres =
    elemManager->ConstructMaterialViewAccessor< array4d<real64> >( MultiFluidBase::
                                                                   viewKeyStruct::
                                                                   dPhaseCompFraction_dPressureString,
                                                                   constitutiveManager );

  auto dPhaseCompFrac_dComp =
    elemManager->ConstructMaterialViewAccessor< array5d<real64> >( MultiFluidBase::
                                                                   viewKeyStruct::
                                                                   dPhaseCompFraction_dGlobalCompFractionString,
                                                                   constitutiveManager );

  auto dCompFrac_dCompDens =
    elemManager->ConstructViewAccessor<array3d<real64>>( viewKeysCompMultiphaseFlow.dGlobalCompFraction_dGlobalCompDensity.Key() );

  localIndex const NC   = m_numComponents;
  localIndex const NP   = m_numPhases;
  localIndex const NDOF = m_numDofPerCell;

  constexpr localIndex numElems = 2; // number of connected elements
  array1d<long long> eqnRowIndices( numElems * NC );
  array1d<long long> dofColIndices( numElems * NDOF ); // to be resized for stencil size
  array1d<double> localFlux( numElems * NC );
  array2d<double> localFluxJacobian( numElems * NC, numElems * NDOF ); // to be resized for stencil size

  // temporary working arrays

  real64 densWeight[numElems] = { 0.5, 0.5 };

  // these arrays have constant size

  array1d<real64> compFlux( NC );
  array1d<real64> dRelPerm_dC( NC );

  array1d<real64> dDens_dC( NC );
  array1d<real64> dVisc_dC( NC );

  array1d<real64> mobility( numElems );
  array1d<real64> dMobility_dP( numElems );
  array2d<real64> dMobility_dC( numElems, NC );

  array1d<real64> dDensMean_dP( numElems );
  array2d<real64> dDensMean_dC( numElems, NC );

  // the arrays below are resized for each cell's stencil size

  array1d<real64> dPhaseFlux_dP( numElems );
  array2d<real64> dPhaseFlux_dC( numElems, NC );

  array2d<real64> dCompFlux_dP( numElems, NC );
  array3d<real64> dCompFlux_dC( numElems, NC, NC );

  stencilCollection.forAll<RAJA::seq_exec>([=] (StencilCollection<CellDescriptor, real64>::Accessor stencil) mutable -> void
  {
    localIndex const stencilSize = stencil.size();

    // reset the local values
    compFlux = 0.0;
    dCompFlux_dP = 0.0;
    dCompFlux_dC = 0.0;

    localFlux = 0.0;
    localFluxJacobian = 0.0;

    // resize local working arrays that are stencil-dependent
    dDensMean_dP.resize( stencilSize ); // doesn't need to be that large, but it's convenient
    dDensMean_dC.resizeDimension<0>( stencilSize );
    dPhaseFlux_dP.resize( stencilSize );
    dPhaseFlux_dC.resizeDimension<0>( stencilSize );
    dCompFlux_dP.resizeDimension<0>( stencilSize );
    dCompFlux_dC.resizeDimension<0>( stencilSize );

    // resize local matrices and vectors
    dofColIndices.resize( stencilSize * NDOF );
    localFluxJacobian.resizeDimension<1>( stencilSize * NDOF );

    // set equation indices for both connected cells
    stencil.forConnected( [&] ( auto const & cell,
                                localIndex i ) -> void
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

      // calculate quantities on primary connected cells
      stencil.forConnected( [&] ( auto const & cell,
                                 localIndex i ) -> void
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
        real64 const relPerm = 1.0; // TODO
        real64 dRelPerm_dP = 0.0;
        dRelPerm_dC = 0.0;
        for (localIndex jp = 0; jp < NP; ++jp)
        {
          real64 const dRelPerm_dS = 0.0; // TODO
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
      real64 potDif = 0.0;
      stencil.forAll( [&] ( CellDescriptor cell,
                            real64 w,
                            localIndex i ) -> void
      {
        localIndex const er  = cell.region;
        localIndex const esr = cell.subRegion;
        localIndex const ei  = cell.index;

        globalIndex const offset = NDOF * blockLocalDofNumber[er][esr][ei];
        for (localIndex jdof = 0; jdof < NDOF; ++jdof)
        {
          dofColIndices[i * NDOF + jdof] = offset + jdof;
        }

        real64 const gravD = gravDepth[er][esr][ei];
        real64 const gravTerm = m_gravityFlag ? densMean * gravD : 0.0;
        real64 const dGrav_dP = m_gravityFlag ? dDensMean_dP[i] * gravD : 0.0;

        potDif += w * (pres[er][esr][ei] + dPres[er][esr][ei] + gravTerm);
        dPhaseFlux_dP[i] = w * (1.0 + dGrav_dP);

        for (localIndex jc = 0; jc < NC; ++jc)
        {
          real64 const dGrav_dC = m_gravityFlag ? dDensMean_dC[i][jc] * gravD : 0.0;
          dPhaseFlux_dC[i][jc] = w * dGrav_dC;
        }
      });

      // upwinding of fluid properties
      localIndex const k_up = (potDif >= 0) ? 0 : 1;

      // compute the phase flux and derivatives
      real64 const phaseFlux = mobility[k_up] * potDif;
      for (localIndex ke = 0; ke < stencilSize; ++ke)
      {
        dPhaseFlux_dP[ke] *= mobility[k_up];
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dPhaseFlux_dC[ke][jc] *= mobility[k_up];
        }
      }
      dPhaseFlux_dP[k_up] += dMobility_dP[k_up] * potDif;
      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dPhaseFlux_dC[k_up][jc] += dMobility_dC[k_up][jc] * potDif;
      }

      CellDescriptor cell_up = stencil.connectedIndex( k_up );
      localIndex er_up  = cell_up.region;
      localIndex esr_up = cell_up.subRegion;
      localIndex ei_up  = cell_up.index;

      // compute component fluxes and derivatives
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        real64 const ycp = phaseCompFrac[er_up][esr_up][m_fluidIndex][ei_up][0][ip][ic];
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

        // additional derivatives stemming from upwinding of phase composition
        dCompFlux_dP[k_up][ic] += phaseFlux * dPhaseCompFrac_dPres[er_up][esr_up][m_fluidIndex][ei_up][0][ip][ic];
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dCompFlux_dC[k_up][ic][jc] += phaseFlux * dPhaseCompFrac_dComp[er_up][esr_up][m_fluidIndex][ei_up][0][ip][ic][jc];
        }
      }
    }

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

void CompositionalMultiphaseFlow::AssembleVolumeBalanceTerms( DomainPartition * const domain,
                                                              EpetraBlockSystem * const blockSystem,
                                                              real64 const time_n,
                                                              real64 const dt )
{
//***** extract data required for assembly of system *****
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  ConstitutiveManager const * constitutiveManager =
    domain->GetGroup<ConstitutiveManager>( keys::ConstitutiveManager );

  ElementRegionManager const * elemManager = mesh->getElemManager();

  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::compositionalBlock,
                                                                BlockIDs::compositionalBlock );
  Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );

  auto elemGhostRank =
    elemManager->ConstructViewAccessor<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

  auto blockLocalDofNumber =
    elemManager->ConstructViewAccessor<array1d<globalIndex >>( viewKeysCompMultiphaseFlow.blockLocalDofNumber.Key() );

  auto pres      = elemManager->ConstructViewAccessor<array1d<real64>>( viewKeysCompMultiphaseFlow.pressure.Key() );
  auto dPres     = elemManager->ConstructViewAccessor<array1d<real64>>( viewKeysCompMultiphaseFlow.deltaPressure.Key() );
  auto volume    = elemManager->ConstructViewAccessor<array1d<real64>>( CellBlock::viewKeyStruct::elementVolumeString );

  auto phaseVolFrac = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.phaseVolumeFraction.Key() );
  auto dPhaseVolFrac_dPres = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.dPhaseVolumeFraction_dPressure.Key() );
  auto dPhaseVolFrac_dComp = elemManager->ConstructViewAccessor<array3d<real64>>( viewKeysCompMultiphaseFlow.dPhaseVolumeFraction_dGlobalCompDensity.Key() );

  auto porosityRef =
    elemManager->ConstructViewAccessor<array1d<real64>>( viewKeysCompMultiphaseFlow.referencePorosity.Key() );

  auto pVMult =
    elemManager->ConstructMaterialViewAccessor< array2d<real64> >( ConstitutiveBase::
                                                                   viewKeyStruct::
                                                                   poreVolumeMultiplierString,
                                                                   constitutiveManager );

  auto dPVMult_dPres =
    elemManager->ConstructMaterialViewAccessor< array2d<real64> >( ConstitutiveBase::
                                                                   viewKeyStruct::
                                                                   dPVMult_dPresString,
                                                                   constitutiveManager );

  localIndex const NC   = m_numComponents;
  localIndex const NP   = m_numPhases;
  localIndex const NDOF = m_numDofPerCell;

  // using Epetra types
  array1d<long long> localVolBalanceDOF( NDOF );
  array1d<double> localVolBalanceJacobian( NDOF );

  //***** Loop over all elements and assemble the change in volume/density terms *****
//  forAllElemsInMesh( mesh, [=] ( localIndex const er,
//                                 localIndex const esr,
//                                 localIndex const ei ) -> void

  for (localIndex er = 0; er < elemManager->numRegions(); ++er)
  {
    ElementRegion const * const elemRegion = elemManager->GetRegion(er);
    for (localIndex esr = 0; esr < elemRegion->numSubRegions(); ++esr)
    {
      CellBlockSubRegion const * const cellBlockSubRegion = elemRegion->GetSubRegion(esr);

      // set up array views to reduce amount of indexing
      arrayView1d<globalIndex const> const dofNumber = blockLocalDofNumber[er][esr].get();

      arrayView2d<real64 const> const pVMultSub = pVMult[er][esr][m_solidIndex].get();
      arrayView2d<real64 const> const dPVMult_dPresSub = dPVMult_dPres[er][esr][m_solidIndex].get();

      arrayView2d<real64 const> const phaseVolFracSub = phaseVolFrac[er][esr].get();
      arrayView2d<real64 const> const dPhaseVolFrac_dPresSub = dPhaseVolFrac_dPres[er][esr].get();
      arrayView3d<real64 const> const dPhaseVolFrac_dCompSub = dPhaseVolFrac_dComp[er][esr].get();

      for (localIndex ei = 0; ei < cellBlockSubRegion->size(); ++ei)
      {
        if (elemGhostRank[er][esr][ei] < 0)
        {
          // compute pore volume
          real64 const vol      = volume[er][esr][ei];
          real64 const dVol_dP  = 0.0; // used in poroelastic solver

          real64 const poro     = porosityRef[er][esr][ei] * pVMultSub[ei][0];
          real64 const dPoro_dP = porosityRef[er][esr][ei] * dPVMult_dPresSub[ei][0];

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
            localVolBalance -= phaseVolFracSub[ei][ip];
            localVolBalanceJacobian[0] -= dPhaseVolFrac_dPresSub[ei][ip];

            for (localIndex jc = 0; jc < NC; ++jc)
            {
              localVolBalanceJacobian[jc+1] -= dPhaseVolFrac_dCompSub[ei][ip][jc];
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
        }
      }
    }
  }//)
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

  jacobian->GlobalAssemble(true);
  residual->GlobalAssemble();

  if (verboseLevel() >= 2)
  {
    jacobian->Print( std::cout );
    residual->Print( std::cout );
  }
}

void
CompositionalMultiphaseFlow::ApplyDirichletBC_implicit( DomainPartition * domain,
                                                        real64 const time_n, real64 const dt,
                                                        EpetraBlockSystem * const blockSystem )
{
  BoundaryConditionManager * bcManager = BoundaryConditionManager::get();

  // 1. apply pressure Dirichlet BCs
  bcManager->ApplyBoundaryCondition( time_n + dt,
                                     domain,
                                     "ElementRegions",
                                     viewKeysCompMultiphaseFlow.pressure.Key(),
                                     [&]( BoundaryConditionBase const * const bc,
                                          string const &,
                                          set<localIndex> const & targetSet,
                                          ManagedGroup * subRegion,
                                          string const & ) -> void
  {
    // 1.1. Apply BC to set the field values
    bc->ApplyBoundaryConditionToField<BcEqual>( targetSet,
                                                time_n + dt,
                                                subRegion,
                                                viewKeysCompMultiphaseFlow.bcPressure.Key() );
  });

  // TODO how to check consistency between pressure and composition BC?

  // 2. Apply composition BC (global component fraction) and store them for constitutive call
  bcManager->ApplyBoundaryCondition( time_n + dt,
                                     domain,
                                     "ElementRegions",
                                     viewKeysCompMultiphaseFlow.globalCompFraction.Key(),
                                     [&] ( BoundaryConditionBase const * const bc,
                                           string const &,
                                           set<localIndex> const & targetSet,
                                           ManagedGroup * subRegion,
                                           string const & ) -> void
  {
    // 2.1. Apply BC to set the field values
    bc->ApplyBoundaryConditionToField<BcEqual>( targetSet,
                                                time_n + dt,
                                                subRegion,
                                                viewKeysCompMultiphaseFlow.globalCompFraction.Key() );
  });

  // 3. Call constitutive update on boundary set cells to get total density
  bcManager->ApplyBoundaryCondition( time_n + dt,
                                     domain,
                                     "ElementRegions",
                                     viewKeysCompMultiphaseFlow.pressure.Key(),
                                     [&] ( BoundaryConditionBase const * const bc,
                                           string const & setName,
                                           set<localIndex> const & targetSet,
                                           ManagedGroup * subRegion,
                                           string const & ) -> void
  {
    UpdateFluidModel( subRegion, targetSet,
                      viewKeyStruct::bcPressureString,
                      "",
                      viewKeyStruct::globalCompFractionString );
  });

  // 4. Back-calculate target global component densities and apply to the system
  bcManager->ApplyBoundaryCondition( time_n + dt,
                                     domain,
                                     "ElementRegions",
                                     viewKeysCompMultiphaseFlow.pressure.Key(),
                                     [&] ( BoundaryConditionBase const * const bc,
                                           string const & setName,
                                           set<localIndex> const & targetSet,
                                           ManagedGroup * targetGroup,
                                           string const & ) -> void
  {
    // TODO for now can only be a subregion; this will change in the future
    CellBlockSubRegion * subRegion = targetGroup->group_cast<CellBlockSubRegion *>();

    GEOS_ERROR_IF( subRegion == nullptr, "Target data group is not a subregion" );

    ManagedGroup * const constitutiveGroup = subRegion->GetConstitutiveModels();
    MultiFluidBase * const fluid = constitutiveGroup->GetGroup<MultiFluidBase>( m_fluidIndex ); // TODO could be incorrect

    GEOS_ERROR_IF( fluid == nullptr, "Fluid model does not exist in subregion " << subRegion->getName() );

    auto const & pres      = subRegion->getReference<array1d<real64>>( viewKeysCompMultiphaseFlow.pressure.Key() );
    auto const & dPres     = subRegion->getReference<array1d<real64>>( viewKeysCompMultiphaseFlow.deltaPressure.Key() );
    auto const & bcPres    = subRegion->getReference<array1d<real64>>( viewKeysCompMultiphaseFlow.bcPressure.Key() );
    auto const & compFrac  = subRegion->getReference<array2d<real64>>( viewKeysCompMultiphaseFlow.globalCompFraction.Key() );
    auto const & compDens  = subRegion->getReference<array2d<real64>>( viewKeysCompMultiphaseFlow.globalCompDensity.Key() );
    auto const & dCompDens = subRegion->getReference<array2d<real64>>( viewKeysCompMultiphaseFlow.deltaGlobalCompDensity.Key() );

    auto const & totalDens = fluid->getReference<array2d<real64>>( MultiFluidBase::viewKeyStruct::totalDensityString );

    auto const & dofNumber = subRegion->getReference<array1d<globalIndex>>( viewKeysCompMultiphaseFlow.blockLocalDofNumber.Key() );

    Epetra_FEVector * const rhs = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );
    array1d<real64> rhsContribution( targetSet.size() * m_numDofPerCell );
    array1d<globalIndex> dof( targetSet.size() * m_numDofPerCell );

    integer counter = 0;

    for (localIndex a : targetSet)
    {
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

        BcEqual::ApplyBcValue( dof[counter],
                               blockSystem,
                               BlockIDs::compositionalBlock,
                               rhsContribution[counter],
                               totalDens[a][0] * compFrac[a][ic],
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
CompositionalMultiphaseFlow::ApplyFaceDirichletBC_implicit( DomainPartition * domain,
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

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  auto elemGhostRank =
    elemManager->ConstructViewAccessor<integer_array>( ObjectManagerBase::viewKeyStruct::ghostRankString );

  auto blockLocalDofNumber = elemManager->
    ConstructViewAccessor<array1d<globalIndex>>(viewKeysCompMultiphaseFlow.blockLocalDofNumber.Key());

  auto refPoro = elemManager->ConstructViewAccessor<real64_array>(viewKeysCompMultiphaseFlow.referencePorosity.Key());
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
      globalIndex const offset = m_numDofPerCell * blockLocalDofNumber[er][esr][ei];
      for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
      {
        int const lid = rowMap->LID(integer_conversion<int>(offset + idof));
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

bool
CompositionalMultiphaseFlow::CheckSystemSolution( EpetraBlockSystem const * const blockSystem,
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
    elemManager->ConstructViewAccessor<array1d<globalIndex>>( viewKeysCompMultiphaseFlow.blockLocalDofNumber.Key() );

  auto pres      = elemManager->ConstructViewAccessor<array1d<real64>>( viewKeysCompMultiphaseFlow.pressure.Key() );
  auto dPres     = elemManager->ConstructViewAccessor<array1d<real64>>( viewKeysCompMultiphaseFlow.deltaPressure.Key() );
  auto compDens  = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.globalCompDensity.Key() );
  auto dCompDens = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.deltaGlobalCompDensity.Key() );

  auto elemGhostRank =
    elemManager->ConstructViewAccessor<integer_array>( ObjectManagerBase::viewKeyStruct::ghostRankString );

  bool result = true;

  // TODO use reduction over mesh
  // loop over all elements to update incremental pressure
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei )->void
  {
    if( elemGhostRank[er][esr][ei] < 0 )
    {
      globalIndex const offset = m_numDofPerCell * blockLocalDofNumber[er][esr][ei];
      // extract solution and apply to dP
      {
        int const lid = rowMap->LID(integer_conversion<int>(offset));
        real64 const newPres = pres[er][esr][ei] + dPres[er][esr][ei] + scalingFactor * local_solution[lid];
        if (newPres < 0.0)
        {
          result = false;
        }
      }

      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        int const lid = rowMap->LID(integer_conversion<int>(offset + ic + 1));
        real64 const newDens = compDens[er][esr][ei][ic] + dCompDens[er][esr][ei][ic] + scalingFactor * local_solution[lid];
        if (newDens < 0.0)
        {
          result = false;
        }
      }
    }
  });

  return result;
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
    elemManager->ConstructViewAccessor<array1d<globalIndex>>( viewKeysCompMultiphaseFlow.blockLocalDofNumber.Key() );

  auto dPres     = elemManager->ConstructViewAccessor<array1d<real64>>( viewKeysCompMultiphaseFlow.deltaPressure.Key() );
  auto dCompDens = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.deltaGlobalCompDensity.Key() );

  auto elemGhostRank =
    elemManager->ConstructViewAccessor<integer_array>( ObjectManagerBase::viewKeyStruct::ghostRankString );

  // loop over all elements to update incremental pressure
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei )->void
  {
    if( elemGhostRank[er][esr][ei] < 0 )
    {
      globalIndex const offset = m_numDofPerCell * blockLocalDofNumber[er][esr][ei];
      // extract solution and apply to dP
      {
        int const lid = rowMap->LID(integer_conversion<int>(offset));
        dPres[er][esr][ei] += scalingFactor * local_solution[lid];
      }

      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        int const lid = rowMap->LID(integer_conversion<int>(offset + ic + 1));
        dCompDens[er][esr][ei][ic] += scalingFactor * local_solution[lid];
      }
    }
  });

  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeysCompMultiphaseFlow.deltaPressure.Key() );
  fieldNames["elems"].push_back( viewKeysCompMultiphaseFlow.deltaGlobalCompDensity.Key() );
  CommunicationTools::SynchronizeFields( fieldNames,
                                         mesh,
                                         domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  UpdateComponentFraction( domain );
  UpdateConstitutiveModels( domain );
  UpdatePhaseVolumeFraction( domain );
}

void CompositionalMultiphaseFlow::ResetStateToBeginningOfStep(DomainPartition * const domain)
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  auto dPres     = elemManager->ConstructViewAccessor<array1d<real64>>(viewKeysCompMultiphaseFlow.deltaPressure.Key());
  auto dCompDens = elemManager->ConstructViewAccessor<array2d<real64>>(viewKeysCompMultiphaseFlow.deltaGlobalCompDensity.Key());

  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei )->void
  {
    dPres[er][esr][ei] = 0.0;
    for (localIndex ic = 0; ic < m_numComponents; ++ic)
      dCompDens[er][esr][ei][ic] = 0.0;
  });

  UpdateComponentFraction( domain );
  UpdateConstitutiveModels( domain );
  UpdatePhaseVolumeFraction( domain );
}

void CompositionalMultiphaseFlow::ImplicitStepComplete( real64 const & time,
                                                        real64 const & dt,
                                                        DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  auto pres      = elemManager->ConstructViewAccessor<array1d<real64>>(viewKeysCompMultiphaseFlow.pressure.Key());
  auto dPres     = elemManager->ConstructViewAccessor<array1d<real64>>(viewKeysCompMultiphaseFlow.deltaPressure.Key());
  auto compDens  = elemManager->ConstructViewAccessor<array2d<real64>>(viewKeysCompMultiphaseFlow.globalCompDensity.Key());
  auto dCompDens = elemManager->ConstructViewAccessor<array2d<real64>>(viewKeysCompMultiphaseFlow.deltaGlobalCompDensity.Key());

  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei )->void
  {
    pres[er][esr][ei] += dPres[er][esr][ei];
    for (localIndex ic = 0; ic < m_numComponents; ++ic)
      compDens[er][esr][ei][ic] += dCompDens[er][esr][ei][ic];
  });
}

bool CompositionalMultiphaseFlow::TestNumericalJacobian( DomainPartition * domain,
                                                         EpetraBlockSystem * blockSystem,
                                                         real64 const time_n,
                                                         real64 const dt,
                                                         double perturbParameter,
                                                         double relTol )
{
  Epetra_FECrsMatrix const * jacobian = blockSystem->GetMatrix( BlockIDs::compositionalBlock, BlockIDs::compositionalBlock );
  Epetra_FEVector const * residual = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );
  Epetra_Map      const * rowMap   = blockSystem->GetRowMap( BlockIDs::compositionalBlock );

  // get a view into local residual vector
  int localSizeInt;
  double* localResidual = nullptr;
  residual->ExtractView(&localResidual, &localSizeInt);

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  auto elemGhostRank =
    elemManager->ConstructViewAccessor<integer_array>( ObjectManagerBase::viewKeyStruct::ghostRankString );

  auto pres      = elemManager->ConstructViewAccessor<array1d<real64>>( viewKeysCompMultiphaseFlow.pressure.Key() );
  auto dPres     = elemManager->ConstructViewAccessor<array1d<real64>>( viewKeysCompMultiphaseFlow.deltaPressure.Key() );
  auto compDens  = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.globalCompDensity.Key() );
  auto dCompDens = elemManager->ConstructViewAccessor<array2d<real64>>( viewKeysCompMultiphaseFlow.deltaGlobalCompDensity.Key() );

  auto blockLocalDofNumber =
    elemManager->ConstructViewAccessor<array1d<globalIndex>>( viewKeysCompMultiphaseFlow.blockLocalDofNumber.Key() );

  // assemble the analytical residual
  ResetStateToBeginningOfStep( domain );
  AssembleSystem( domain, blockSystem, time_n, dt );

  // copy the analytical residual
  auto residualOrig = std::make_unique<Epetra_FEVector>( *residual );
  double* localResidualOrig = nullptr;
  residualOrig->ExtractView(&localResidualOrig, &localSizeInt);

  // create the numerical jacobian
  auto jacobianFD = std::make_unique<Epetra_FECrsMatrix>( Copy, jacobian->Graph() );
  jacobianFD->Scale( 0.0 );

  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei ) -> void
  {
    if (elemGhostRank[er][esr][ei] >= 0)
      return;

    globalIndex offset = blockLocalDofNumber[er][esr][ei] * m_numDofPerCell;

    real64 totalDensity = 0.0;
    for (localIndex ic = 0; ic < m_numComponents; ++ic)
    {
      totalDensity += compDens[er][esr][ei][ic];
    }

    {
      ResetStateToBeginningOfStep(domain);
      real64 const dP = perturbParameter * std::min(std::fabs(pres[er][esr][ei]), 1e5);
      dPres[er][esr][ei] = dP;
      AssembleSystem(domain, blockSystem, time_n, dt);
      long long const dofIndex = integer_conversion<long long>(offset);

      for (int lid = 0; lid < localSizeInt; ++lid)
      {
        real64 dRdP = (localResidualOrig[lid] - localResidual[lid]) / dP;
        if (std::fabs(dRdP) > 0.0)
        {
          long long gid = rowMap->GID64(lid);
          jacobianFD->ReplaceGlobalValues(gid, 1, &dRdP, &dofIndex);
        }
      }
    }

    for (localIndex ic = 0; ic < m_numComponents; ++ic)
    {
      ResetStateToBeginningOfStep(domain);
      real64 const dRho = perturbParameter * totalDensity;
      dCompDens[er][esr][ei][ic] = dRho;
      AssembleSystem(domain, blockSystem, time_n, dt);
      long long const dofIndex = integer_conversion<long long>(offset + ic + 1);

      for (int lid = 0; lid < localSizeInt; ++lid)
      {
        real64 dRdRho = (localResidualOrig[lid] - localResidual[lid]) / dRho;
        if (std::fabs(dRdRho) > 0.0)
        {
          long long gid = rowMap->GID64(lid);
          jacobianFD->ReplaceGlobalValues(gid, 1, &dRdRho, &dofIndex);
        }
      }
    }
  });

  jacobianFD->GlobalAssemble(true);

  // assemble the analytical jacobian
  ResetStateToBeginningOfStep( domain );
  AssembleSystem( domain, blockSystem, time_n, dt );

#if 1
  jacobian->Print(std::cout);
  jacobianFD->Print(std::cout);
#endif

  bool result = true;
  double * row = nullptr;
  double * rowFD = nullptr;
  int numEntries, numEntriesFD;
  int * indices = nullptr;
  int* indicesFD = nullptr;

  // check the accuracy across local rows
  for (int i = 0; i < localSizeInt; ++i)
  {
    jacobian->ExtractMyRowView( i, numEntries, row, indices );
    jacobianFD->ExtractMyRowView( i, numEntriesFD, rowFD, indicesFD );

    if (numEntries != numEntriesFD)
    {
      GEOS_LOG( "Number of entries in local row " << i
                << " of analytical and numerical jacobians does not match: "
                << numEntries << " != " << numEntriesFD);
      result = false;
    }
    for (int j = 0, jFD = 0; j < numEntries && jFD < numEntriesFD; ++j, ++jFD)
    {
      while (j < numEntries && jFD < numEntriesFD && indices[j] != indices[jFD])
      {
        while (j < numEntries && indices[j] < indicesFD[jFD])
        {
          GEOS_LOG( "Entry (" << i << ", " << indices[j] << ") in analytical jacobian does not have a match" );
          result = false;
          j++;
        }
        while (jFD < numEntriesFD && indicesFD[jFD] < indices[j])
        {
          GEOS_LOG( "Entry (" << i << ", " << indicesFD[jFD] << ") in numerical jacobian does not have a match" );
          result = false;
          jFD++;
        }
      }
      if (j < numEntries && jFD < numEntriesFD)
      {
        double const delta = std::fabs(row[j] - rowFD[j]);
        double const value = std::fmax(std::fabs(row[j]), std::fabs(rowFD[jFD]));
        if (value > 0 && delta / value > relTol)
        {
          GEOS_LOG( "Entry (" << i << ", " << indices[j] << ") mismatch in rel tol: " << delta/value << " > " << relTol );
          result = false;
        }
      }
    }
  }

  return result;
}


REGISTER_CATALOG_ENTRY(SolverBase, CompositionalMultiphaseFlow, string const &, ManagedGroup * const)
}// namespace geosx
