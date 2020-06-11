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

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

FlowSolverBase::FlowSolverBase( std::string const & name,
                                Group * const parent ):
  SolverBase( name, parent ),
  m_fluidModelNames(),
  m_solidModelNames(),
  m_poroElasticFlag( 0 ),
  m_coupledWellsFlag( 0 ),
  m_numDofPerCell( 0 ),
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
  this->registerWrapper( viewKeyStruct::discretizationString, &m_discretizationName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of discretization object to use for this solver." );

  this->registerWrapper( viewKeyStruct::fluidNamesString, &m_fluidModelNames )->
    setInputFlag( InputFlags::REQUIRED )->
    setSizedFromParent( 0 )->
    setDescription( "Names of fluid constitutive models for each region." );

  this->registerWrapper( viewKeyStruct::solidNamesString, &m_solidModelNames )->
    setInputFlag( InputFlags::REQUIRED )->
    setSizedFromParent( 0 )->
    setDescription( "Names of solid constitutive models for each region." );

  this->registerWrapper( viewKeyStruct::inputFluxEstimateString, &m_fluxEstimate )->
    setApplyDefaultValue( 1.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Initial estimate of the input flux used only for residual scaling. This should be "
                    "essentially equivalent to the input flux * dt." );

  this->registerWrapper( viewKeyStruct::meanPermCoeffString, &m_meanPermCoeff )->
    setApplyDefaultValue( 1.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Coefficient to move between harmonic mean (1.0) and arithmetic mean (0.0) for the "
                    "calculation of permeability between elements." );

}

void FlowSolverBase::RegisterDataOnMesh( Group * const MeshBodies )
{
  SolverBase::RegisterDataOnMesh( MeshBodies );

  for( auto & subgroup : MeshBodies->GetSubGroups() )
  {
    MeshBody & meshBody = *subgroup.second->group_cast< MeshBody * >();
    MeshLevel & mesh = *meshBody.getMeshLevel( 0 );

    forTargetSubRegions( mesh, [&]( localIndex const,
                                    ElementSubRegionBase & subRegion )
    {
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::referencePorosityString )->setPlotLevel( PlotLevel::LEVEL_0 );
      subRegion.registerWrapper< array1d< R1Tensor > >( viewKeyStruct::permeabilityString )->setPlotLevel( PlotLevel::LEVEL_0 );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::gravityCoefString )->setApplyDefaultValue( 0.0 );
    } );

    ElementRegionManager * const elemManager = mesh.getElemManager();

    elemManager->forElementSubRegionsComplete< FaceElementSubRegion >( [&]( localIndex const,
                                                                            localIndex const,
                                                                            ElementRegionBase & region,
                                                                            FaceElementSubRegion & subRegion )
    {
      FaceElementRegion & faceRegion = dynamicCast< FaceElementRegion & >( region );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::referencePorosityString )->
        setApplyDefaultValue( 1.0 );

      subRegion.registerWrapper< array1d< R1Tensor > >( viewKeyStruct::permeabilityString )->setPlotLevel( PlotLevel::LEVEL_0 );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::gravityCoefString )->setApplyDefaultValue( 0.0 );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::aperture0String )->
        setDefaultValue( faceRegion.getDefaultAperture() );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::effectiveApertureString )->
        setApplyDefaultValue( subRegion.getWrapper< array1d< real64 > >( FaceElementSubRegion::
                                                                           viewKeyStruct::
                                                                           elementApertureString )->getDefaultValue() )->
        setPlotLevel( PlotLevel::LEVEL_0 );
    } );

    FaceManager * const faceManager = mesh.getFaceManager();
    faceManager->registerWrapper< array1d< real64 > >( viewKeyStruct::gravityCoefString )->setApplyDefaultValue( 0.0 );
  }
}

void FlowSolverBase::PostProcessInput()
{
  SolverBase::PostProcessInput();
  CheckModelNames( m_fluidModelNames, viewKeyStruct::fluidNamesString );
  CheckModelNames( m_solidModelNames, viewKeyStruct::solidNamesString );
}

void FlowSolverBase::InitializePreSubGroups( Group * const rootGroup )
{
  SolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup< DomainPartition >( keys::domain );

  // Validate solid models in regions (fluid models are validated by derived classes)
  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel & meshLevel = *Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );
    ValidateModelMapping( *meshLevel.getElemManager(), m_solidModelNames );
  }

  // fill stencil targetRegions
  NumericalMethodsManager & numericalMethodManager = domain->getNumericalMethodManager();

  FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();

  FluxApproximationBase * const fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  array1d< string > & stencilTargetRegions = fluxApprox->targetRegions();
  std::set< string > stencilTargetRegionsSet( stencilTargetRegions.begin(), stencilTargetRegions.end() );
  for( auto const & targetRegion : targetRegionNames() )
  {
    stencilTargetRegionsSet.insert( targetRegion );
  }

  stencilTargetRegions.clear();
  for( auto const & targetRegion : stencilTargetRegionsSet )
  {
    stencilTargetRegions.push_back( targetRegion );
  }
}

void FlowSolverBase::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  SolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup< DomainPartition >( keys::domain );

  ResetViews( domain );

  // Precompute solver-specific constant data (e.g. gravity-depth)
  PrecomputeData( domain );
}

void FlowSolverBase::PrecomputeData( DomainPartition * const domain )
{
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  FaceManager & faceManager = *mesh.getFaceManager();

  R1Tensor const gravVector = gravityVector();

  forTargetSubRegions( mesh, [&]( localIndex const,
                                  ElementSubRegionBase & subRegion )
  {
    arrayView2d< real64 const > const & elemCenter = subRegion.getElementCenter();

    arrayView1d< real64 > const & gravityCoef =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex a )
    {
      gravityCoef[ a ] = elemCenter( a, 0 ) * gravVector[ 0 ] + elemCenter( a, 1 ) * gravVector[ 1 ] + elemCenter( a, 2 ) * gravVector[ 2 ];
    } );
  } );

  {
    arrayView2d< real64 const > const & faceCenter = faceManager.faceCenter();

    arrayView1d< real64 > const & gravityCoef =
      faceManager.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

    forAll< serialPolicy >( faceManager.size(), [=] ( localIndex a )
    {
      // TODO change to LvArray::tensorOps::AiBi once gravVector is a c-array.
      gravityCoef[ a ] = faceCenter[ a ][ 0 ] * gravVector[ 0 ];
      gravityCoef[ a ] += faceCenter[ a ][ 1 ] * gravVector[ 1 ];
      gravityCoef[ a ] += faceCenter[ a ][ 2 ] * gravVector[ 2 ];
    } );
  }
}

FlowSolverBase::~FlowSolverBase() = default;

void FlowSolverBase::ResetViews( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();

  m_elemGhostRank =
    elemManager->ConstructViewAccessor< array1d< integer >, arrayView1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );
  m_volume =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString );
  m_gravCoef =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::gravityCoefString );
  m_porosityRef =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::referencePorosityString );

  m_elementArea =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( FaceElementSubRegion::viewKeyStruct::elementAreaString );
  m_elementAperture =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( FaceElementSubRegion::viewKeyStruct::elementApertureString );
  m_elementAperture0 =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::aperture0String );

  m_effectiveAperture =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::effectiveApertureString );

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  m_elementSeparationCoefficient =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( FaceElementSubRegion::viewKeyStruct::separationCoeffString );

  m_element_dSeparationCoefficient_dAperture =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( FaceElementSubRegion::viewKeyStruct::dSeparationCoeffdAperString );
#endif
}

} // namespace geosx
