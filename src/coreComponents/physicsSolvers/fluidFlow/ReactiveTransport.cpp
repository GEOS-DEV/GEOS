/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ProppantTransport.cpp
 */

#include "ReactiveTransport.hpp"

#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"

#include "physicsSolvers/fluidFlow/ReactiveTransportKernels.hpp"
#include "constitutive/fluid/ThermoDatabases/KineticReactionsBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace ReactiveTransportKernels;

ReactiveTransport::ReactiveTransport( const std::string & name,
                                      Group * const parent ):
  FlowSolverBase( name, parent )
{

  this->registerWrapper( viewKeyStruct::viscosityString, &m_viscosity )->setApplyDefaultValue( 0.001 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Fluid viscosity" );

  this->registerWrapper( viewKeyStruct::densityString, &m_density )->setApplyDefaultValue( 988.52 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Fluid density" );


  this->registerWrapper( GeochemicalModel::viewKeyStruct::reactiveFluidNamesString, &m_reactiveFluidNames )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of chemical system constitutive objects to use each target region." );

  this->registerWrapper( viewKeyStruct::permPoroPowerString, &m_permPoroPower )->setApplyDefaultValue( 0.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Perm-Porosity power" );

  this->registerWrapper( viewKeyStruct::updatePorosityString, &m_updatePorosity )->setApplyDefaultValue( 1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Update porosity" );

  this->registerWrapper( viewKeyStruct::maxChangeString, &m_maxChange )->setApplyDefaultValue( 3.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Max change" );


}

void ReactiveTransport::PostProcessInput()
{
  FlowSolverBase::PostProcessInput();

}

void ReactiveTransport::RegisterDataOnMesh( Group * const MeshBodies )
{
  FlowSolverBase::RegisterDataOnMesh( MeshBodies );

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel & meshLevel = *Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    forTargetSubRegions< CellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                                 CellElementSubRegion & subRegion )
    {

      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::bcComponentConcentrationString )->setDefaultValue( 0.0 );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::initialPorosityString );

      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::initialMineralSurfaceAreaString );

      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::initialMineralVolumeFractionString );

      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::mineralVolumeFractionString )->setPlotLevel( PlotLevel::LEVEL_0 );

    } );
  }

}

void ReactiveTransport::InitializePreSubGroups( Group * const rootGroup )
{
  FlowSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup< DomainPartition >( keys::domain );

  MeshLevel & meshLevel = *(domain->getMeshBody( 0 )->getMeshLevel( 0 ));
  this->forTargetSubRegions( meshLevel,
                             [&]
                               ( localIndex const targetRegionIndex,
                               ElementSubRegionBase const & subRegion )
  {
    string const & fluidName = m_reactiveFluidNames[targetRegionIndex];
    ReactiveFluidBase const & reactiveFluid = *(subRegion.getConstitutiveModel< ReactiveFluidBase >( fluidName ) );

    m_numComponents = reactiveFluid.numBasisSpecies();
    m_numDofPerCell = m_numComponents;
    m_numKineticReactions = reactiveFluid.numKineticReaction();

  } );

  ResizeFields( meshLevel );

}

void ReactiveTransport::ResizeFields( MeshLevel & mesh )
{

  localIndex const NC = m_numComponents;
  localIndex const NR = m_numKineticReactions;

  forTargetSubRegions( mesh,
                       [&]
                         ( localIndex const,
                         ElementSubRegionBase & subRegion )
  {

    subRegion.getReference< array2d< real64 > >( viewKeyStruct::bcComponentConcentrationString ).resizeDimension< 1 >( NC );

    subRegion.getReference< array2d< real64 > >( viewKeyStruct::initialMineralSurfaceAreaString ).resizeDimension< 1 >( NR );
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::initialMineralVolumeFractionString ).resizeDimension< 1 >( NR );
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::mineralVolumeFractionString ).resizeDimension< 1 >( NR );

  } );

}

void ReactiveTransport::UpdateReactiveFluidModel( Group * const dataGroup, localIndex const targetIndex )const
{
  GEOSX_MARK_FUNCTION;

  ReactiveFluidBase & reactiveFluid = GetConstitutiveModel< ReactiveFluidBase >( *dataGroup, m_reactiveFluidNames[targetIndex] );

  arrayView1d< real64 const > const & pres = dataGroup->getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dPres = dataGroup->getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  arrayView1d< real64 const > const & temp = dataGroup->getReference< array1d< real64 > >( GeochemicalModel::viewKeyStruct::temperatureString );
  arrayView1d< real64 const > const & dTemp = dataGroup->getReference< array1d< real64 > >( GeochemicalModel::viewKeyStruct::deltaTemperatureString );

  arrayView2d< real64 const > const & conc = dataGroup->getReference< array2d< real64 > >( GeochemicalModel::viewKeyStruct::concentrationString );
  arrayView2d< real64 const > const & dConc = dataGroup->getReference< array2d< real64 > >( GeochemicalModel::viewKeyStruct::deltaConcentrationString );
  arrayView2d< real64 > & concNew = dataGroup->getReference< array2d< real64 > >( GeochemicalModel::viewKeyStruct::concentrationNewString );

  arrayView1d< real64 const > const & porosity0 = dataGroup->getReference< array1d< real64 > >( viewKeyStruct::initialPorosityString );

  arrayView1d< real64 const > const & porosity = dataGroup->getReference< array1d< real64 > >( viewKeyStruct::referencePorosityString );

  arrayView2d< real64 const > const & A0 = dataGroup->getReference< array2d< real64 > >( viewKeyStruct::initialMineralSurfaceAreaString );

  arrayView2d< real64 const > const & theta0 = dataGroup->getReference< array2d< real64 > >( viewKeyStruct::initialMineralVolumeFractionString );

  arrayView2d< real64 const > const & theta = dataGroup->getReference< array2d< real64 > >( viewKeyStruct::mineralVolumeFractionString );

  forAll< serialPolicy >( dataGroup->size(), [&] ( localIndex const a )
  {
    for( localIndex ic = 0; ic < m_numComponents; ++ic )
    {
      concNew[a][ic] = conc[a][ic] + dConc[a][ic];
    }

    reactiveFluid.PointUpdateChemistry( pres[a] + dPres[a], temp[a] + dTemp[a], concNew[a], A0[a], theta0[a], theta[a], porosity0[a], porosity[a], a );
  } );

}


void ReactiveTransport::UpdateState( Group * dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  UpdateReactiveFluidModel( dataGroup, targetIndex );

}

void ReactiveTransport::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition & domain = *rootGroup->GetGroup< DomainPartition >( keys::domain );
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( GeochemicalModel::viewKeyStruct::deltaConcentrationString ) );

  CommunicationTools::SynchronizeFields( fieldNames, &mesh, domain.getNeighbors(), true );

  ResetViews( mesh );

}

real64 ReactiveTransport::SolverStep( real64 const & time_n,
                                      real64 const & dt,
                                      const int cycleNumber,
                                      DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  ImplicitStepSetup( time_n, dt, domain );

  if( cycleNumber == 0 )
  {
    FieldSpecificationManager const & boundaryConditionManager = FieldSpecificationManager::get();

    boundaryConditionManager.ApplyInitialConditions( &domain );

  }

  PreStepUpdate( time_n, dt, domain );

  // currently the only method is implicit time integration
  real64 const dtReturn= this->NonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.

  ImplicitStepComplete( time_n, dtReturn, domain );

  PostStepUpdate( time_n, dtReturn, domain );

  return dtReturn;
}


void ReactiveTransport::PreStepUpdate( real64 const & GEOSX_UNUSED_PARAM( time ),
                                       real64 const & GEOSX_UNUSED_PARAM( dt ),
                                       DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  FlowSolverBase::PrecomputeData( mesh );

  localIndex const NC = m_numComponents;

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion )
  {

    arrayView2d< real64 > const deltaConc = subRegion.getReference< array2d< real64 > >( GeochemicalModel::viewKeyStruct::deltaConcentrationString );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      for( localIndex c = 0; c < NC; ++c )
      {
        deltaConc[ei][c] = 0.0;
      }
    } );
  } );

}

void ReactiveTransport::PostStepUpdate( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                        real64 const & dt_return,
                                        DomainPartition & domain )
{

  UpdateRockProperties( dt_return, domain );

}

void ReactiveTransport::ImplicitStepSetup( real64 const & time_n,
                                           real64 const & GEOSX_UNUSED_PARAM( dt ),
                                           DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  ResetViews( mesh );

  forTargetSubRegionsComplete( mesh,
                               [&]
                                 ( localIndex const targetRegionIndex,
                                 localIndex const er,
                                 localIndex const esr,
                                 ElementRegionBase & GEOSX_UNUSED_PARAM( region ),
                                 ElementSubRegionBase & subRegion )
  {

    arrayView1d< real64 > const & dPres   = m_deltaPressure[er][esr];
    arrayView1d< real64 > const & dTemp   = m_deltaTemperature[er][esr];

    arrayView2d< real64 > const & dConc = m_deltaConcentration[er][esr];

    arrayView2d< real64 > const & componentConc = m_componentConcentration[er][esr];
    arrayView2d< real64 const > const totalConc = m_totalConc[er][esr];

    forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
    {
      dPres[ei] = 0.0;
      dTemp[ei] = 0.0;
      for( localIndex ic = 0; ic < m_numComponents; ++ic )
      {
        dConc[ei][ic] = 0.0;
      }

    } );

    forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
    {

      for( localIndex ic = 0; ic < m_numComponents; ++ic )
      {

        componentConc[ei][ic] = totalConc[ei][ic];

      }

    } );


    UpdateState( &subRegion, targetRegionIndex );

    forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
    {

      for( localIndex ic = 0; ic < m_numComponents; ++ic )
      {

        componentConc[ei][ic] = totalConc[ei][ic];

      }

    } );


  } );

  if( time_n <= 0.0 )
  {

    forTargetSubRegionsComplete( mesh,
                                 [&]
                                   ( localIndex const,
                                   localIndex const,
                                   localIndex const,
                                   ElementRegionBase & GEOSX_UNUSED_PARAM( region ),
                                   ElementSubRegionBase & subRegion )
    {

      arrayView1d< R1Tensor > const transTMult =
        subRegion.getReference< array1d< R1Tensor > >( SinglePhaseBase::viewKeyStruct::transTMultString );

      forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
      {
        transTMult[ei][0] = 1.0;
        transTMult[ei][1] = 1.0;
        transTMult[ei][2] = 1.0;

      } );



    } );

  }

}

void ReactiveTransport::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                              real64 const & GEOSX_UNUSED_PARAM( dt ),
                                              DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  localIndex const NC = m_numComponents;

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion )
  {

    arrayView2d< real64 > const conc =
      subRegion.getReference< array2d< real64 > >( GeochemicalModel::viewKeyStruct::concentrationString );
    arrayView2d< real64 const > const dConc =
      subRegion.getReference< array2d< real64 > >( GeochemicalModel::viewKeyStruct::deltaConcentrationString );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {

      for( localIndex c = 0; c < NC; ++c )
      {
        conc[ei][c] += dConc[ei][c];

      }
    } );
  } );


  forTargetSubRegionsComplete( mesh,
                               [&]
                                 ( localIndex const targetRegionIndex,
                                 localIndex const er,
                                 localIndex const esr,
                                 ElementRegionBase & GEOSX_UNUSED_PARAM( region ),
                                 ElementSubRegionBase & subRegion )
  {


    ReactiveFluidBase & reactiveFluid  = GetConstitutiveModel< ReactiveFluidBase >( subRegion, m_reactiveFluidNames[targetRegionIndex] );

    arrayView1d< bool > const isHplus = reactiveFluid.IsHplus();

    arrayView2d< real64 > const & conc = m_concentration[er][esr];

    arrayView2d< real64 > const & concOut =
      subRegion.getReference< array2d< real64 > >( GeochemicalModel::viewKeyStruct::concentrationOutString );

    forAll< serialPolicy >( subRegion.size(),
                            [=]
                              ( localIndex const ei )
    {
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        if( isHplus[ic] )
          concOut[ei][ic] = -conc[ei][ic];
        else
          concOut[ei][ic] = pow( 10.0, conc[ei][ic] );

      }
    } );

  } );


}

void ReactiveTransport::SetupDofs( DomainPartition const & domain,
                                   DofManager & dofManager ) const
{

  dofManager.addField( viewKeyStruct::reactiveTransportModelString,
                       DofManager::Location::Elem,
                       m_numDofPerCell,
                       targetRegionNames() );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  dofManager.addCoupling( viewKeyStruct::reactiveTransportModelString, fluxApprox );

}


void ReactiveTransport::AssembleSystem( real64 const time,
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  if( 1 )
    AssembleAccumulationTerms( dt,
                               domain,
                               dofManager,
                               localMatrix,
                               localRhs );
  if( 1 )
    AssembleFluxTerms( time,
                       dt,
                       domain,
                       dofManager,
                       localMatrix,
                       localRhs );


}

void ReactiveTransport::AssembleAccumulationTerms( real64 const dt,
                                                   DomainPartition const & domain,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegionsComplete( mesh,
                               [&]
                                 ( localIndex const,
                                 localIndex const er,
                                 localIndex const esr,
                                 ElementRegionBase const & GEOSX_UNUSED_PARAM( region ),
                                 ElementSubRegionBase const & subRegion )
  {

    string const dofKey = dofManager.getKey( viewKeyStruct::reactiveTransportModelString );
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const & volume = subRegion.getElementVolume();

    arrayView2d< real64 const > const componentConc = m_componentConcentration[er][esr];

    arrayView2d< real64 const > const totalConc = m_totalConc[er][esr];
    arrayView3d< real64 const > const dTotalConc_dConc = m_dTotalConc_dConc[er][esr];

    arrayView2d< real64 const > const kineticSpeciesReactionRate = m_kineticSpeciesReactionRate[er][esr];
    arrayView3d< real64 const > const dKineticSpeciesReactionRate_dConc = m_dKineticSpeciesReactionRate_dConc[er][esr];

    arrayView1d< real64 const > const & porosity = subRegion.getReference< array1d< real64 > >( viewKeyStruct::referencePorosityString );

    AccumulationKernel::Launch( subRegion.size(),
                                m_numComponents,
                                m_numDofPerCell,
                                dofManager.rankOffset(),
                                dofNumber,
                                elemGhostRank,
                                componentConc,
                                totalConc,
                                dTotalConc_dConc,
                                kineticSpeciesReactionRate,
                                dKineticSpeciesReactionRate_dConc,
                                porosity,
                                volume,
                                m_density,
                                dt,
                                localMatrix,
                                localRhs );
  } );
}


void ReactiveTransport::AssembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                           real64 const dt,
                                           DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                           arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const & elemManager = *mesh.getElemManager();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const dofKey = dofManager.getKey( viewKeyStruct::reactiveTransportModelString );

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  dofNumberAccessor = elemManager.ConstructViewAccessor< array1d< globalIndex >, arrayView1d< globalIndex const > >( dofKey );

  FluxKernel::ElementViewConst< arrayView1d< globalIndex const > > const dofNumber = dofNumberAccessor.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const pres  = m_pressure.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const dPres = m_deltaPressure.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const totalConc   = m_totalConc.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView3d< real64 const > > const dTotalConc_dConc  = m_dTotalConc_dConc.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView1d< integer const > > const elemGhostRank = m_elemGhostRank.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView1d< R1Tensor const > > const transTMult   = m_transTMult.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView1d< R1Tensor const > > const permeability   = m_permeability.toNestedViewConst();

  fluxApprox.forStencils< CellElementStencilTPFA >( mesh, [&]( auto const & stencil )
  {

    FluxKernel::Launch( stencil,
                        m_numDofPerCell,
                        dt,
                        dofManager.rankOffset(),
                        dofNumber,
                        elemGhostRank,
                        pres,
                        dPres,
                        totalConc,
                        dTotalConc_dConc,
                        transTMult,
                        permeability,
                        m_viscosity,
                        localMatrix,
                        localRhs );
  } );
}

void ReactiveTransport::ApplyBoundaryConditions( real64 const time_n,
                                                 real64 const dt,
                                                 DomainPartition & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();
  string const dofKey = dofManager.getKey( viewKeyStruct::reactiveTransportModelString );
  globalIndex const rankOffset = dofManager.rankOffset();

  //  Apply Dirichlet BC for component concentration

  localIndex const NC = m_numComponents;

  if( NC > 0 )
  {
    map< string, map< string, array1d< bool > > > bcStatusMap; // map to check consistent application of BC

    fsManager.Apply( time_n + dt,
                     &domain,
                     "ElementRegions",
                     GeochemicalModel::viewKeyStruct::concentrationString,
                     [&]( FieldSpecificationBase const * const GEOSX_UNUSED_PARAM( fs ),
                          string const & setName,
                          SortedArrayView< localIndex const > const & GEOSX_UNUSED_PARAM( targetSet ),
                          Group * const subRegion,
                          string const & )
    {

      string const & subRegionName = subRegion->getName();
      bcStatusMap[subRegionName][setName].resize( NC );
      bcStatusMap[subRegionName][setName].setValues< serialPolicy >( false );

    } );

    fsManager.Apply( time_n + dt,
                     &domain,
                     "ElementRegions",
                     GeochemicalModel::viewKeyStruct::concentrationString,
                     [&] ( FieldSpecificationBase const * const fs,
                           string const & setName,
                           SortedArrayView< localIndex const > const & targetSet,
                           Group * const subRegion,
                           string const & )
    {

      string const & subRegionName = subRegion->getName();
      localIndex const comp = fs->GetComponent();

      GEOSX_ERROR_IF( bcStatusMap[subRegionName][setName][comp], "Conflicting composition[" << comp << "] boundary conditions on set '" << setName << "'" );
      bcStatusMap[subRegionName][setName][comp] = true;

      fs->ApplyFieldValue< FieldSpecificationEqual >( targetSet,
                                                      time_n + dt,
                                                      subRegion,
                                                      viewKeyStruct::bcComponentConcentrationString );

    } );

    bool bcConsistent = true;
    for( auto const & bcStatusEntryOuter : bcStatusMap )
    {
      for( auto const & bcStatusEntryInner : bcStatusEntryOuter.second )
      {
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          bcConsistent &= bcStatusEntryInner.second[ic];
          GEOSX_WARNING_IF( !bcConsistent, "Composition boundary condition not applied to component " << ic
                                                                                                      << " on region '" << bcStatusEntryOuter.first << "',"
                                                                                                      << " set '" << bcStatusEntryInner.first << "'" );
        }
      }
    }

    GEOSX_ERROR_IF( !bcConsistent, "Inconsistent composition boundary conditions" );

    fsManager.Apply( time_n + dt,
                     &domain,
                     "ElementRegions",
                     GeochemicalModel::viewKeyStruct::concentrationString,
                     [&] ( FieldSpecificationBase const * const GEOSX_UNUSED_PARAM( bc ),
                           string const & GEOSX_UNUSED_PARAM( setName ),
                           SortedArrayView< localIndex const > const & targetSet,
                           Group * const subRegion,
                           string const & )
    {
      arrayView1d< integer const > const ghostRank =
        subRegion->getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );
      arrayView1d< globalIndex const > const dofNumber = subRegion->getReference< array1d< globalIndex > >( dofKey );

      arrayView2d< real64 const > const conc =
        subRegion->getReference< array2d< real64 > >( GeochemicalModel::viewKeyStruct::concentrationString );
      arrayView2d< real64 const > const deltaConc =
        subRegion->getReference< array2d< real64 > >( GeochemicalModel::viewKeyStruct::deltaConcentrationString );
      arrayView2d< real64 const > const bcCompConc =
        subRegion->getReference< array2d< real64 > >( viewKeyStruct::bcComponentConcentrationString );

      forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
      {
        localIndex const ei = targetSet[a];
        if( ghostRank[ei] >= 0 )
          return;

        globalIndex const dofIndex = dofNumber[ei];
        localIndex const localRow = dofIndex - rankOffset;
        real64 rhsValue;

        for( localIndex ic = 0; ic < NC; ++ic )
        {
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + ic,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      bcCompConc[ei][ic],
                                                      conc[ei][ic] + deltaConc[ei][ic] );
          localRhs[localRow + ic] = rhsValue;
        }
      } );
    } );
  }
}

real64
ReactiveTransport::CalculateResidualNorm( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                          DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                          arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( localRhs ) )
{

  //  need to comput global maxDSol for parallel run

  return m_maxDSol;

}

void ReactiveTransport::ApplySystemSolution( DofManager const & dofManager,
                                             arrayView1d< real64 const > const & localSolution,
                                             real64 const,
                                             DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  localIndex const NDOF = m_numDofPerCell;

  localIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::reactiveTransportModelString );

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion )
  {
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();

    arrayView2d< real64 > const & dConc =
      subRegion.getReference< array2d< real64 > >( GeochemicalModel::viewKeyStruct::deltaConcentrationString );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( elemGhostRank[ei] < 0 )
      {
        localIndex const lid = dofNumber[ei] - rankOffset;
        for( localIndex idof = 0; idof < NDOF; ++idof )
        {
          real64 val = localSolution[lid + idof];

          if( fabs( val ) > m_maxChange )
            val = val / fabs( val ) * m_maxChange;


          dConc[ei][idof] += val;
        }
      }
    } );

  } );


  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( GeochemicalModel::viewKeyStruct::deltaConcentrationString ) );

  CommunicationTools::SynchronizeFields( fieldNames, &mesh, domain.getNeighbors(), true );

  this->forTargetSubRegions( mesh, [&] ( localIndex const targetRegionIndex,
                                         ElementSubRegionBase & subRegion )
  {
    UpdateState( &subRegion, targetRegionIndex );
  } );


  forTargetSubRegionsComplete( mesh,
                               [&]
                                 ( localIndex const,
                                 localIndex const,
                                 localIndex const,
                                 ElementRegionBase & GEOSX_UNUSED_PARAM( region ),
                                 ElementSubRegionBase & subRegion )
  {

    localIndex num = 0.0;

    m_maxDSol = -1e10;

    for( localIndex ie = 0; ie < subRegion.size(); ie++ )
    {

      for( localIndex ic = 0; ic < m_numComponents; ++ic )
      {

        if( fabs( localSolution[num] ) > m_maxDSol )
          m_maxDSol = fabs( localSolution[num] );

        num++;

      }

    }

  } );

}

void ReactiveTransport::SolveSystem( DofManager const & dofManager,
                                     ParallelMatrix & matrix,
                                     ParallelVector & rhs,
                                     ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );
}

void ReactiveTransport::ResetStateToBeginningOfStep( DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  localIndex const NC = m_numComponents;

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion )
  {
    arrayView2d< real64 > const & dComponentConc =
      subRegion.getReference< array2d< real64 > >( GeochemicalModel::viewKeyStruct::deltaConcentrationString );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      for( localIndex c = 0; c < NC; ++c )
      {
        dComponentConc[ei][c] = 0.0;
      }
    } );

  } );

}

void ReactiveTransport::ResetViews( MeshLevel & mesh )
{
  FlowSolverBase::ResetViews( mesh );
  ElementRegionManager & elemManager = *mesh.getElemManager();

  m_pressure.clear();
  m_pressure = elemManager.ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::pressureString );

  m_deltaPressure.clear();
  m_deltaPressure = elemManager.ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::deltaPressureString );

  m_temperature.clear();
  m_temperature = elemManager.ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( GeochemicalModel::viewKeyStruct::temperatureString );

  m_deltaTemperature.clear();
  m_deltaTemperature = elemManager.ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( GeochemicalModel::viewKeyStruct::deltaTemperatureString );

  m_concentration.clear();
  m_concentration = elemManager.ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( GeochemicalModel::viewKeyStruct::concentrationString );

  m_deltaConcentration.clear();
  m_deltaConcentration = elemManager.ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( GeochemicalModel::viewKeyStruct::deltaConcentrationString );

  m_componentConcentration.clear();
  m_componentConcentration = elemManager.ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( GeochemicalModel::viewKeyStruct::totalConcentrationString );

  m_transTMult.clear();
  m_transTMult = elemManager.ConstructArrayViewAccessor< R1Tensor, 1 >( SinglePhaseBase::viewKeyStruct::transTMultString );

  m_permeability.clear();
  m_permeability = elemManager.ConstructArrayViewAccessor< R1Tensor, 1 >( SinglePhaseBase::viewKeyStruct::permeabilityString );


  m_totalConc.clear();
  m_totalConc =
    elemManager.ConstructMaterialArrayViewAccessor< real64, 2 >( ReactiveFluidBase::viewKeyStruct::totalConcString,
                                                                 targetRegionNames(),
                                                                 m_reactiveFluidNames );

  m_dTotalConc_dConc.clear();
  m_dTotalConc_dConc =
    elemManager.ConstructMaterialArrayViewAccessor< real64, 3 >( ReactiveFluidBase::viewKeyStruct::dTotalConc_dConcString,
                                                                 targetRegionNames(),
                                                                 m_reactiveFluidNames );

  m_kineticSpeciesReactionRate.clear();
  m_kineticSpeciesReactionRate =
    elemManager.ConstructMaterialArrayViewAccessor< real64, 2 >( ReactiveFluidBase::viewKeyStruct::kineticSpeciesReactionRateString,
                                                                 targetRegionNames(),
                                                                 m_reactiveFluidNames );

  m_dKineticSpeciesReactionRate_dConc.clear();
  m_dKineticSpeciesReactionRate_dConc =
    elemManager.ConstructMaterialArrayViewAccessor< real64, 3 >( ReactiveFluidBase::viewKeyStruct::dKineticSpeciesReactionRate_dConcString,
                                                                 targetRegionNames(),
                                                                 m_reactiveFluidNames );

  m_kineticReactionRate.clear();
  m_kineticReactionRate =
    elemManager.ConstructMaterialArrayViewAccessor< real64, 2 >( ReactiveFluidBase::viewKeyStruct::kineticReactionRateString,
                                                                 targetRegionNames(),
                                                                 m_reactiveFluidNames );

}

void ReactiveTransport::UpdateRockProperties( real64 const & dt,
                                              DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *(domain.getMeshBody( 0 )->getMeshLevel( 0 ));

  forTargetSubRegionsComplete( mesh,
                               [&]
                                 ( localIndex const targetRegionIndex,
                                 localIndex const er,
                                 localIndex const esr,
                                 ElementRegionBase & GEOSX_UNUSED_PARAM( region ),
                                 ElementSubRegionBase & subRegion )
  {

    ReactiveFluidBase & reactiveFluid  = GetConstitutiveModel< ReactiveFluidBase >( subRegion, m_reactiveFluidNames[targetRegionIndex] );

    arrayView1d< real64 > const & porosity = subRegion.getReference< array1d< real64 > >( viewKeyStruct::referencePorosityString );

    arrayView2d< real64 > const & theta = subRegion.getReference< array2d< real64 > >( viewKeyStruct::mineralVolumeFractionString );

    arrayView1d< R1Tensor > const transTMultiplier =
      subRegion.getReference< array1d< R1Tensor > >( SinglePhaseBase::viewKeyStruct::transTMultString );

    arrayView1d< real64 const > const & porosity0 = subRegion.getReference< array1d< real64 > >( viewKeyStruct::initialPorosityString );

    arrayView2d< real64 const > const & kineticReactionRate = m_kineticReactionRate[er][esr];

    const array1d< KineticReaction > & kineticReactions = reactiveFluid.GetKineticReactions();

    forAll< serialPolicy >( subRegion.size(), [&] ( localIndex const a )
    {

      real64 totalChange = 0.0;
      real64 dTheta;

      for( localIndex ir = 0; ir < kineticReactions.size(); ++ir )
      {

        const KineticReaction & kineticReaction = kineticReactions[ir];


        dTheta = kineticReactionRate[a][ir] * kineticReaction.MW / kineticReaction.density * dt;

        real64 thetaNew = theta[a][ir] + dTheta;
        if( thetaNew < 0.0 )
          dTheta = -theta[a][ir];

        totalChange += dTheta;
        theta[a][ir] += dTheta;

      }

      if( m_updatePorosity == 1 )
      {
        porosity[a] -= totalChange;
        if( porosity[a] < 0.0 )
          porosity[a] = 0.0;
        if( porosity[a] > 1.0 )
          porosity[a] = 1.0;
      }

      real64 TransMult = pow( porosity[a] / porosity0[a], m_permPoroPower );

      transTMultiplier[a] = TransMult;

    } );

  } );

}


REGISTER_CATALOG_ENTRY( SolverBase, ReactiveTransport, std::string const &, Group * const )
} /* namespace geosx */
