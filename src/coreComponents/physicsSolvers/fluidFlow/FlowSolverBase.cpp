/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FlowSolverBase.cpp
 */

#include "FlowSolverBase.hpp"

#include "constitutive/ConstitutivePassThru.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "physicsSolvers/fluidFlow/FluxKernelsHelper.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

template< typename POROUSWRAPPER_TYPE >
void execute1( POROUSWRAPPER_TYPE porousWrapper,
               CellElementSubRegion & subRegion,
               arrayView1d< real64 const > const & pressure,
               arrayView1d< real64 const > const & deltaPressure )
{
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < porousWrapper.numGauss(); ++q )
    {
      porousWrapper.updateStateFromPressure( k, q,
                                             pressure[k],
                                             deltaPressure[k] );
    }
  } );
}

template< typename POROUSWRAPPER_TYPE >
void execute2( POROUSWRAPPER_TYPE porousWrapper,
               SurfaceElementSubRegion & subRegion,
               arrayView1d< real64 const > const & pressure,
               arrayView1d< real64 const > const & deltaPressure,
               arrayView1d< real64 const > const & oldHydraulicAperture,
               arrayView1d< real64 const > const & newHydraulicAperture )
{
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < porousWrapper.numGauss(); ++q )
    {
      porousWrapper.updateStateFromPressureAndAperture( k, q,
                                                        pressure[k],
                                                        deltaPressure[k],
                                                        oldHydraulicAperture[k],
                                                        newHydraulicAperture[k] );
    }
  } );
}

FlowSolverBase::FlowSolverBase( string const & name,
                                Group * const parent ):
  SolverBase( name, parent ),
  m_poroElasticFlag( 0 ),
  m_coupledWellsFlag( 0 ),
  m_numDofPerCell( 0 ),
  m_fluxEstimate(),
  m_elemGhostRank(),
  m_volume(),
  m_gravCoef()
{
  this->registerWrapper( viewKeyStruct::discretizationString(), &m_discretizationName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of discretization object to use for this solver." );

  this->registerWrapper( viewKeyStruct::inputFluxEstimateString(), &m_fluxEstimate ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Initial estimate of the input flux used only for residual scaling. This should be "
                    "essentially equivalent to the input flux * dt." );
}

void FlowSolverBase::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  forMeshTargets( meshBodies, [&] ( string const &,
                                    MeshLevel & mesh,
                                    arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::gravityCoefString() ).
        setApplyDefaultValue( 0.0 );

      subRegion.registerWrapper< string >( viewKeyStruct::fluidNamesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setSizedFromParent(0);

      subRegion.registerWrapper< string >( viewKeyStruct::solidNamesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setSizedFromParent(0);

      subRegion.registerWrapper< string >( viewKeyStruct::permeabilityNamesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setSizedFromParent(0);
    } );

    elemManager.forElementSubRegionsComplete< SurfaceElementSubRegion >( [&]( localIndex const,
                                                                              localIndex const,
                                                                              ElementRegionBase & region,
                                                                              SurfaceElementSubRegion & subRegion )
    {
      SurfaceElementRegion & faceRegion = dynamicCast< SurfaceElementRegion & >( region );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::gravityCoefString() ).
        setApplyDefaultValue( 0.0 );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::aperture0String() ).
        setDefaultValue( faceRegion.getDefaultAperture() );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::effectiveApertureString() ).
        setApplyDefaultValue( faceRegion.getDefaultAperture() ).
        setPlotLevel( PlotLevel::LEVEL_0 );
    } );

    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerWrapper< array1d< real64 > >( viewKeyStruct::gravityCoefString() ).setApplyDefaultValue( 0.0 );
  } );
}

void FlowSolverBase::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d<string const> const & regionNames )
  {
    ElementRegionManager & elementRegionManager = mesh.getElemManager();
    elementRegionManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                       [&]( localIndex const,
                                                                            ElementSubRegionBase & subRegion )
    {
      string & solidName = subRegion.getReference<string>( viewKeyStruct::solidNamesString() );
      solidName = getConstitutiveName< CoupledSolidBase >( subRegion );
      GEOSX_ERROR_IF( solidName.empty(), GEOSX_FMT( "Solid model not found on subregion {}", subRegion.getName() ) );

      string & permName = subRegion.getReference<string>( viewKeyStruct::permeabilityNamesString() );
      permName = getConstitutiveName< PermeabilityBase >( subRegion );
      GEOSX_ERROR_IF( solidName.empty(), GEOSX_FMT( "Permeability model not found on subregion {}", subRegion.getName() ) );
    } );
  } );

  // fill stencil targetRegions
  NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();

  if( fvManager.hasGroup< FluxApproximationBase >( m_discretizationName ) )
  {
    forMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                  MeshLevel & mesh,
                                                  arrayView1d< string const > const & regionNames )
    {
      FluxApproximationBase & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
      array1d< string > & stencilTargetRegions = fluxApprox.targetRegions( meshBodyName );
      array1d< string > & stencilCoeffModelNames = fluxApprox.coefficientModelNames( meshBodyName );

      std::set< string > stencilTargetRegionsSet( stencilTargetRegions.begin(), stencilTargetRegions.end() );
      map< string, string > coeffModelNames;

      mesh.getElemManager().forElementSubRegionsComplete< ElementSubRegionBase >( regionNames,
                                                                                  [&]( localIndex const,
                                                                                       localIndex const,
                                                                                       localIndex const,
                                                                                       ElementRegionBase const & region,
                                                                                       ElementSubRegionBase const & subRegion )
      {
        string const & permName = subRegion.getReference<string>( viewKeyStruct::permeabilityNamesString() );
        coeffModelNames[region.getName()] = permName;
        stencilTargetRegionsSet.insert( region.getName() );
      } );

      stencilTargetRegions.clear();
      stencilCoeffModelNames.clear();
      for( auto const & targetRegion: stencilTargetRegionsSet )
      {
        stencilTargetRegions.emplace_back( targetRegion );
        stencilCoeffModelNames.emplace_back( coeffModelNames[targetRegion] );
      }
    } );
  }
}

void FlowSolverBase::initializePostInitialConditionsPreSubGroups()
{
  SolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    resetViews( mesh ); // TODO: won't work with multiple meshes
    precomputeData( mesh, regionNames );
  } );
}

void FlowSolverBase::precomputeData( MeshLevel & mesh,
                                     arrayView1d< string const > const & regionNames )
{
  FaceManager & faceManager = mesh.getFaceManager();
  real64 const gravVector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

  mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames, [&]( localIndex const,
                                                                                        ElementSubRegionBase & subRegion )
  {
    arrayView2d< real64 const > const elemCenter = subRegion.getElementCenter();

    arrayView1d< real64 > const gravityCoef =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString() );

    forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const ei )
    {
      gravityCoef[ ei ] = LvArray::tensorOps::AiBi< 3 >( elemCenter[ ei ], gravVector );
    } );
  } );

  {
    arrayView2d< real64 const > const faceCenter = faceManager.faceCenter();

    arrayView1d< real64 > const gravityCoef =
      faceManager.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString() );

    forAll< parallelHostPolicy >( faceManager.size(), [=] ( localIndex const kf )
    {
      gravityCoef[ kf ] = LvArray::tensorOps::AiBi< 3 >( faceCenter[ kf ], gravVector );
    } );
  }
}

void FlowSolverBase::updatePorosityAndPermeability( CellElementSubRegion & subRegion ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const & pressure = subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
  arrayView1d< real64 const > const & deltaPressure = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );

  string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
  CoupledSolidBase & porousSolid = subRegion.template getConstitutiveModel< CoupledSolidBase >( solidName );

  constitutive::ConstitutivePassThru< CoupledSolidBase >::execute( porousSolid, [=, &subRegion] ( auto & castedPorousSolid )
  {
    typename TYPEOFREF( castedPorousSolid ) ::KernelWrapper porousWrapper = castedPorousSolid.createKernelUpdates();

    execute1( porousWrapper, subRegion, pressure, deltaPressure );
  } );
}

void FlowSolverBase::updatePorosityAndPermeability( SurfaceElementSubRegion & subRegion ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const & pressure = subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
  arrayView1d< real64 const > const & deltaPressure = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );

  arrayView1d< real64 const > const newHydraulicAperture = subRegion.getReference< array1d< real64 > >( viewKeyStruct::effectiveApertureString() );
  arrayView1d< real64 const > const oldHydraulicAperture = subRegion.getReference< array1d< real64 > >( viewKeyStruct:: viewKeyStruct::aperture0String() );

  string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
  CoupledSolidBase & porousSolid = subRegion.getConstitutiveModel< CoupledSolidBase >( solidName );

  constitutive::ConstitutivePassThru< CompressibleSolidBase >::execute( porousSolid, [=, &subRegion] ( auto & castedPorousSolid )
  {
    typename TYPEOFREF( castedPorousSolid ) ::KernelWrapper porousWrapper = castedPorousSolid.createKernelUpdates();

    execute2( porousWrapper, subRegion, pressure, deltaPressure, oldHydraulicAperture, newHydraulicAperture );

  } );
}

void FlowSolverBase::resetViews( MeshLevel & mesh )
{
  ElementRegionManager const & elemManager = mesh.getElemManager();

  m_pressure.clear();
  m_pressure = elemManager.constructArrayViewAccessor< real64, 1 >( viewKeyStruct::pressureString() );
  m_pressure.setName( getName() + "/accessors/" + viewKeyStruct::pressureString() );

  m_deltaPressure.clear();
  m_deltaPressure = elemManager.constructArrayViewAccessor< real64, 1 >( viewKeyStruct::deltaPressureString() );
  m_deltaPressure.setName( getName() + "/accessors/" + viewKeyStruct::deltaPressureString() );

  m_elemGhostRank.clear();
  m_elemGhostRank = elemManager.constructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
  m_elemGhostRank.setName( getName() + "/accessors/" + ObjectManagerBase::viewKeyStruct::ghostRankString() );

  m_volume.clear();
  m_volume = elemManager.constructArrayViewAccessor< real64, 1 >( ElementSubRegionBase::viewKeyStruct::elementVolumeString() );
  m_volume.setName( getName() + "/accessors/" + ElementSubRegionBase::viewKeyStruct::elementVolumeString() );

  m_gravCoef.clear();
  m_gravCoef = elemManager.constructArrayViewAccessor< real64, 1 >( viewKeyStruct::gravityCoefString() );
  m_gravCoef.setName( getName() + "/accessors/" + viewKeyStruct::gravityCoefString() );

  {
    using keys = PermeabilityBase::viewKeyStruct;

    m_permeability.clear();
    m_permeability = elemManager.constructMaterialArrayViewAccessor< PermeabilityBase, real64, 3 >( keys::permeabilityString() );
    m_permeability.setName( getName() + "/accessors/" + keys::permeabilityString() );

    m_dPerm_dPressure.clear();
    m_dPerm_dPressure = elemManager.constructMaterialArrayViewAccessor< PermeabilityBase, real64, 3 >( keys::dPerm_dPressureString() );
    m_dPerm_dPressure.setName( getName() + "/accessors/" + keys::dPerm_dPressureString() );
  }

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  m_elementSeparationCoefficient.clear();
  m_elementSeparationCoefficient = elemManager.constructArrayViewAccessor< real64, 1 >( FaceElementSubRegion::viewKeyStruct::separationCoeffString() );
  m_elementSeparationCoefficient.setName( getName() + "/accessors/" + FaceElementSubRegion::viewKeyStruct::separationCoeffString() );

  m_element_dSeparationCoefficient_dAperture.clear();
  m_element_dSeparationCoefficient_dAperture = elemManager.constructArrayViewAccessor< real64, 1 >(
    FaceElementSubRegion::viewKeyStruct::dSeparationCoeffdAperString() );
  m_element_dSeparationCoefficient_dAperture.setName( getName() + "/accessors/" + FaceElementSubRegion::viewKeyStruct::dSeparationCoeffdAperString() );
#endif
}

void FlowSolverBase::saveAquiferConvergedState( real64 const & time,
                                                real64 const & dt,
                                                DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  // This step requires three passes:
  //    - First we count the number of individual aquifers
  //    - Second we loop over all the stencil entries to compute the sum of aquifer influxes
  //    - Third we loop over the aquifers to save the sums of each individual aquifer

  // Step 1: count individual aquifers

  std::map< string, localIndex > aquiferNameToAquiferId;
  localIndex aquiferCounter = 0;

  fsManager.forSubGroups< AquiferBoundaryCondition >( [&] ( AquiferBoundaryCondition const & bc )
  {
    aquiferNameToAquiferId[bc.getName()] = aquiferCounter;
    aquiferCounter++;
  } );

  // Step 2: sum the aquifer fluxes for each individual aquifer

  array1d< real64 > globalSumFluxes( aquiferNameToAquiferId.size() );
  array1d< real64 > localSumFluxes( aquiferNameToAquiferId.size() );

  fsManager.apply< AquiferBoundaryCondition >( time + dt,
                                               domain,
                                               "faceManager",
                                               AquiferBoundaryCondition::catalogName(),
                                               [&] ( AquiferBoundaryCondition const & bc,
                                                     string const & setName,
                                                     SortedArrayView< localIndex const > const &,
                                                     Group &,
                                                     string const & )
  {
    BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );
    if( stencil.size() == 0 )
    {
      return;
    }

    AquiferBoundaryCondition::KernelWrapper aquiferBCWrapper = bc.createKernelWrapper();

    real64 const targetSetSumFluxes =
      FluxKernelsHelper::AquiferBCKernel::sumFluxes( stencil,
                                                     aquiferBCWrapper,
                                                     m_pressure.toNestedViewConst(),
                                                     m_deltaPressure.toNestedViewConst(),
                                                     m_gravCoef.toNestedViewConst(),
                                                     time,
                                                     dt );

    localIndex const aquiferIndex = aquiferNameToAquiferId.at( bc.getName() );
    localSumFluxes[aquiferIndex] += targetSetSumFluxes;
  } );

  MpiWrapper::allReduce( localSumFluxes.data(),
                         globalSumFluxes.data(),
                         localSumFluxes.size(),
                         MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                         MPI_COMM_GEOSX );

  // Step 3: we are ready to save the summed fluxes for each individual aquifer

  fsManager.forSubGroups< AquiferBoundaryCondition >( [&] ( AquiferBoundaryCondition & bc )
  {
    localIndex const aquiferIndex = aquiferNameToAquiferId.at( bc.getName() );
    bc.saveConvergedState( dt * globalSumFluxes[aquiferIndex] );
  } );
}


} // namespace geosx
