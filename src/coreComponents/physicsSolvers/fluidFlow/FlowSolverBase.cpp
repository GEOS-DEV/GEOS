/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FlowSolverBase.cpp
 */

#include "FlowSolverBase.hpp"

#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

FlowSolverBase::FlowSolverBase( std::string const & name,
                                Group * const parent ):
  SolverBase( name, parent ),
  m_fluidName(),
  m_solidName(),
  m_fluidIndex(),
  m_solidIndex(),
  m_poroElasticFlag(0),
  m_coupledWellsFlag(0),
  m_numDofPerCell(0),
  m_derivativeFluxResidual_dAperture(),
  m_fluxEstimate(),
  m_elemGhostRank(),
  m_volume(),
  m_gravCoef(),
  m_porosityRef(),
  m_elementArea(),
  m_elementAperture0(),
  m_elementAperture()
{
  this->registerWrapper( viewKeyStruct::discretizationString, &m_discretizationName, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of discretization object to use for this solver.");

  this->registerWrapper( viewKeyStruct::fluidNameString,  &m_fluidName,  false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of fluid constitutive object to use for this solver.");

  this->registerWrapper( viewKeyStruct::solidNameString,  &m_solidName,  false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of solid constitutive object to use for this solver");

  this->registerWrapper( viewKeyStruct::fluidIndexString, &m_fluidIndex, false );
  this->registerWrapper( viewKeyStruct::solidIndexString, &m_solidIndex, false );
  
  this->registerWrapper( viewKeyStruct::inputFluxEstimateString,  &m_fluxEstimate,  false )->
    setApplyDefaultValue(1.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Initial estimate of the input flux used only for residual scaling. This should be "
                   "essentially equivalent to the input flux * dt.");
}

void FlowSolverBase::RegisterDataOnMesh( Group * const MeshBodies )
{
  SolverBase::RegisterDataOnMesh( MeshBodies );

  for( auto & subgroup : MeshBodies->GetSubGroups() )
  {
    MeshBody * const meshBody = subgroup.second->group_cast<MeshBody *>();
    MeshLevel * const mesh = meshBody->getMeshLevel(0);

    applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
    {
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::referencePorosityString )->setPlotLevel(PlotLevel::LEVEL_0);
      subRegion->registerWrapper< array1d<R1Tensor> >( viewKeyStruct::permeabilityString )->setPlotLevel(PlotLevel::LEVEL_0);
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::gravityCoefString )->setApplyDefaultValue( 0.0 );
    });

    ElementRegionManager * const elemManager = mesh->getElemManager();

    elemManager->forElementRegions<FaceElementRegion>( [&] ( FaceElementRegion * const region )
    {
      region->forElementSubRegions<FaceElementSubRegion>( [&]( FaceElementSubRegion * const subRegion )
      {
        subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::referencePorosityString )->
          setApplyDefaultValue( 1.0 );

        subRegion->registerWrapper< array1d<R1Tensor> >( viewKeyStruct::permeabilityString )->setPlotLevel(PlotLevel::LEVEL_0);
        subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::gravityCoefString )->setApplyDefaultValue( 0.0 );
        subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::aperture0String )->
          setDefaultValue( region->getDefaultAperture() );
        subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::effectiveApertureString )->
          setApplyDefaultValue( subRegion->getWrapper<array1d<real64>>( FaceElementSubRegion::
                                                                        viewKeyStruct::
                                                                        elementApertureString)->getDefaultValue() )->
          setPlotLevel(PlotLevel::LEVEL_0);
      });
    });

    FaceManager * const faceManager = mesh->getFaceManager();
    faceManager->registerWrapper< array1d<real64> >( viewKeyStruct::gravityCoefString )->setApplyDefaultValue( 0.0 );
  }
}

void FlowSolverBase::InitializePreSubGroups(Group * const rootGroup)
{
  SolverBase::InitializePreSubGroups(rootGroup);

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  ConstitutiveManager * const cm = domain->getConstitutiveManager();

  ConstitutiveBase const * fluid  = cm->GetConstitutiveRelation<ConstitutiveBase>( m_fluidName );
  GEOSX_ERROR_IF( fluid == nullptr, "Fluid model " + m_fluidName + " not found" );
  m_fluidIndex = fluid->getIndexInParent();

  ConstitutiveBase const * solid  = cm->GetConstitutiveRelation<ConstitutiveBase>( m_solidName );
  GEOSX_ERROR_IF( solid == nullptr, "Solid model " + m_solidName + " not found" );
  m_solidIndex = solid->getIndexInParent();

  // fill stencil targetRegions
  NumericalMethodsManager * const
  numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager * const
  fvManager = numericalMethodManager->GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  FluxApproximationBase * const fluxApprox = fvManager->getFluxApproximation( m_discretizationName );
  array1d<string> & stencilTargetRegions = fluxApprox->targetRegions();
  std::set<string> stencilTargetRegionsSet( stencilTargetRegions.begin(), stencilTargetRegions.end() );
  for( auto const & targetRegion : m_targetRegions )
  {
    stencilTargetRegionsSet.insert(targetRegion);
  }

  stencilTargetRegions.clear();
  for( auto const & targetRegion : stencilTargetRegionsSet )
  {
    stencilTargetRegions.push_back( targetRegion );
  }




}

void FlowSolverBase::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  SolverBase::InitializePostInitialConditions_PreSubGroups(rootGroup);

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  ResetViews( domain );

  // Precompute solver-specific constant data (e.g. gravity-depth)
  PrecomputeData(domain);
}

void FlowSolverBase::PrecomputeData( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody(0)->getMeshLevel(0);
  FaceManager * const faceManager = mesh->getFaceManager();

  R1Tensor const gravVector = gravityVector();
  
  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
  {
    arrayView1d<R1Tensor const> const & elemCenter =
      subRegion->getReference<array1d<R1Tensor>>( CellBlock::viewKeyStruct::elementCenterString );

    arrayView1d<real64> const & gravityCoef =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::gravityCoefString );

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex a )
    {
      gravityCoef[a] = Dot( elemCenter[a], gravVector );
    } );
  } );

  {
    arrayView1d<R1Tensor const> const & faceCenter =
      faceManager->getReference<array1d<R1Tensor>>(FaceManager::viewKeyStruct::faceCenterString);

    arrayView1d<real64> const & gravityCoef =
      faceManager->getReference<array1d<real64>>(viewKeyStruct::gravityCoefString);

    forall_in_range<serialPolicy>( 0, faceManager->size(), GEOSX_LAMBDA ( localIndex a )
    {
      gravityCoef[a] = Dot( faceCenter[a], gravVector );
    } );
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
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( ElementSubRegionBase::viewKeyStruct::elementVolumeString );
  m_gravCoef =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::gravityCoefString );
  m_porosityRef =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::referencePorosityString );

  m_elementArea =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( FaceElementSubRegion::viewKeyStruct::elementAreaString );
  m_elementAperture =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( FaceElementSubRegion::viewKeyStruct::elementApertureString );
  m_elementAperture0 =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::aperture0String );

  m_effectiveAperture =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::effectiveApertureString );

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  m_elementSeparationCoefficient =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( FaceElementSubRegion::viewKeyStruct::separationCoeffString );

  m_element_dSeparationCoefficient_dAperture =
      elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( FaceElementSubRegion::viewKeyStruct::dSeparationCoeffdAperString );
#endif
}


} // namespace geosx
