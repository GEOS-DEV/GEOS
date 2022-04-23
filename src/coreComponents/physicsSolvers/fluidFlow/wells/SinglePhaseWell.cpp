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
 * @file SinglePhaseWell.cpp
 */

#include "SinglePhaseWell.hpp"

#include "dataRepository/Group.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "constitutive/fluid/SingleFluidExtrinsicData.hpp"
#include "constitutive/fluid/singleFluidSelector.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "mesh/PerforationData.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWellExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace singlePhaseWellKernels;

SinglePhaseWell::SinglePhaseWell( const string & name,
                                  Group * const parent ):
  WellSolverBase( name, parent )
{
  m_numDofPerWellElement = 2;
}

void SinglePhaseWell::registerDataOnMesh( Group & meshBodies )
{
  WellSolverBase::registerDataOnMesh( meshBodies );

  // loop over the wells
  forMeshTargets( meshBodies, [&] ( string const &,
                                    MeshLevel & mesh,
                                    arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      subRegion.registerExtrinsicData< extrinsicMeshData::well::pressure >( getName() );
      subRegion.registerExtrinsicData< extrinsicMeshData::well::deltaPressure >( getName() ).
        setRestartFlags( RestartFlags::WRITE_AND_READ );

      subRegion.registerExtrinsicData< extrinsicMeshData::well::connectionRate >( getName() );
      subRegion.registerExtrinsicData< extrinsicMeshData::well::deltaConnectionRate >( getName() );

      subRegion.registerExtrinsicData< extrinsicMeshData::well::densityOld >( getName() );

      PerforationData & perforationData = *subRegion.getPerforationData();
      perforationData.registerExtrinsicData< extrinsicMeshData::well::perforationRate >( getName() );
      perforationData.registerExtrinsicData< extrinsicMeshData::well::dPerforationRate_dPres >( getName() ).
        reference().resizeDimension< 1 >( 2 );

      WellControls & wellControls = getWellControls( subRegion );
      wellControls.registerWrapper< real64 >( viewKeyStruct::currentBHPString() );
      wellControls.registerWrapper< real64 >( viewKeyStruct::dCurrentBHP_dPresString() );

      wellControls.registerWrapper< real64 >( viewKeyStruct::currentVolRateString() );
      wellControls.registerWrapper< real64 >( viewKeyStruct::dCurrentVolRate_dPresString() );
      wellControls.registerWrapper< real64 >( viewKeyStruct::dCurrentVolRate_dRateString() );

      string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      fluidName = getConstitutiveName< SingleFluidBase >( subRegion );
      GEOSX_ERROR_IF( fluidName.empty(), GEOSX_FMT( "Fluid model not found on subregion {}", subRegion.getName() ) );

    } );
  } );
}

void SinglePhaseWell::initializePostSubGroups()
{

  WellSolverBase::initializePostSubGroups();
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      validateWellConstraints( subRegion );
    } );
  } );
}

string SinglePhaseWell::resElementDofName() const
{
  return extrinsicMeshData::flow::pressure::key();
}

void SinglePhaseWell::validateWellConstraints( WellElementSubRegion const & subRegion ) const
{
  WellControls const & wellControls = getWellControls( subRegion );
  WellControls::Control const currentControl = wellControls.getControl();
  real64 const targetTotalRate = wellControls.getTargetTotalRate( m_currentTime + m_currentDt );
  real64 const targetPhaseRate = wellControls.getTargetPhaseRate( m_currentTime + m_currentDt );
  GEOSX_THROW_IF( currentControl == WellControls::Control::PHASEVOLRATE,
                  "WellControls named " << wellControls.getName() <<
                  ": Phase rate control is not available for SinglePhaseWell",
                  InputError );
  // The user always provides positive rates, but these rates are later multiplied by -1 internally for producers
  GEOSX_THROW_IF( ( ( wellControls.isInjector() && targetTotalRate < 0.0 ) ||
                    ( wellControls.isProducer() && targetTotalRate > 0.0) ),
                  "WellControls named " << wellControls.getName() <<
                  ": Target total rate cannot be negative",
                  InputError );
  GEOSX_THROW_IF( !isZero( targetPhaseRate ),
                  "WellControls named " << wellControls.getName() <<
                  ": Target phase rate cannot be used for SinglePhaseWell",
                  InputError );
}

void SinglePhaseWell::updateBHPForConstraint( WellElementSubRegion & subRegion )
{
  GEOSX_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.isLocallyOwned() )
  {
    return;
  }

  localIndex const iwelemRef = subRegion.getTopWellElementIndex();

  // subRegion data

  arrayView1d< real64 const > const pres =
    subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
  arrayView1d< real64 const > const dPres =
    subRegion.getExtrinsicData< extrinsicMeshData::well::deltaPressure >();

  arrayView1d< real64 const > const wellElemGravCoef =
    subRegion.getExtrinsicData< extrinsicMeshData::well::gravityCoefficient >();

  // fluid data
  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  SingleFluidBase & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );
  arrayView2d< real64 const > const & dens = fluid.density();
  arrayView2d< real64 const > const & dDens_dPres = fluid.dDensity_dPressure();

  // control data

  WellControls & wellControls = getWellControls( subRegion );

  real64 const & refGravCoef = wellControls.getReferenceGravityCoef();

  real64 & currentBHP =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentBHPString() );
  real64 & dCurrentBHP_dPres =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::dCurrentBHP_dPresString() );

  // bring everything back to host, capture the scalars by reference
  forAll< serialPolicy >( 1, [pres,
                              dPres,
                              dens,
                              dDens_dPres,
                              wellElemGravCoef,
                              &currentBHP,
                              &dCurrentBHP_dPres,
                              &iwelemRef,
                              &refGravCoef] ( localIndex const )
  {
    currentBHP = pres[iwelemRef] + dPres[iwelemRef] + dens[iwelemRef][0] * ( refGravCoef - wellElemGravCoef[iwelemRef] );
    dCurrentBHP_dPres = 1.0 + dDens_dPres[iwelemRef][0] * ( refGravCoef - wellElemGravCoef[iwelemRef] );
  } );
}

void SinglePhaseWell::updateVolRateForConstraint( WellElementSubRegion & subRegion )
{
  GEOSX_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.isLocallyOwned() )
  {
    return;
  }

  localIndex const iwelemRef = subRegion.getTopWellElementIndex();

  // subRegion data

  arrayView1d< real64 const > const pres =
    subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
  arrayView1d< real64 const > const dPres =
    subRegion.getExtrinsicData< extrinsicMeshData::well::deltaPressure >();

  arrayView1d< real64 const > const & connRate =
    subRegion.getExtrinsicData< extrinsicMeshData::well::connectionRate >();
  arrayView1d< real64 const > const & dConnRate =
    subRegion.getExtrinsicData< extrinsicMeshData::well::deltaConnectionRate >();

  // fluid data

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  SingleFluidBase & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );
  arrayView2d< real64 const > const & dens = fluid.density();
  arrayView2d< real64 const > const & dDens_dPres = fluid.dDensity_dPressure();

  // control data

  WellControls & wellControls = getWellControls( subRegion );

  integer const useSurfaceConditions = wellControls.useSurfaceConditions();
  real64 const & surfacePres = wellControls.getSurfacePressure();

  real64 & currentVolRate =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentVolRateString() );
  real64 & dCurrentVolRate_dPres =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::dCurrentVolRate_dPresString() );
  real64 & dCurrentVolRate_dRate =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::dCurrentVolRate_dRateString() );

  constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    // bring everything back to host, capture the scalars by reference
    forAll< serialPolicy >( 1, [fluidWrapper,
                                pres,
                                dPres,
                                connRate,
                                dConnRate,
                                dens,
                                dDens_dPres,
                                &useSurfaceConditions,
                                &surfacePres,
                                &currentVolRate,
                                &dCurrentVolRate_dPres,
                                &dCurrentVolRate_dRate,
                                &iwelemRef] ( localIndex const )
    {
      //    We need to evaluate the density as follows:
      //      - Surface conditions: using the surface pressure provided by the user
      //      - Reservoir conditions: using the pressure in the top element

      if( useSurfaceConditions )
      {
        // we need to compute the surface density
        fluidWrapper.update( iwelemRef, 0, surfacePres );
      }
      else
      {
        real64 const refPres = pres[iwelemRef] + dPres[iwelemRef];
        fluidWrapper.update( iwelemRef, 0, refPres );
      }

      real64 const densInv = 1.0 / dens[iwelemRef][0];
      currentVolRate = ( connRate[iwelemRef] + dConnRate[iwelemRef] ) * densInv;
      dCurrentVolRate_dPres = -( useSurfaceConditions ==  0 ) * dDens_dPres[iwelemRef][0] * currentVolRate * densInv;
      dCurrentVolRate_dRate = densInv;

    } );
  } );
}

void SinglePhaseWell::updateFluidModel( WellElementSubRegion & subRegion ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const pres = subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
  arrayView1d< real64 const > const dPres = subRegion.getExtrinsicData< extrinsicMeshData::well::deltaPressure >();

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  SingleFluidBase & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );

  constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    singlePhaseBaseKernels::FluidUpdateKernel::launch( fluidWrapper, pres, dPres );
  } );
}

void SinglePhaseWell::updateSubRegionState( MeshLevel const & meshLevel, WellElementSubRegion & subRegion )
{
  // update volumetric rates for the well constraints
  // Warning! This must be called before updating the fluid model
  updateVolRateForConstraint( subRegion );

  // update density in the well elements
  updateFluidModel( subRegion );

  // update the current BHP
  updateBHPForConstraint( subRegion );

  // update perforation rates
  computePerforationRates( meshLevel, subRegion );
}

void SinglePhaseWell::initializeWells( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  // loop over the wells
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & meshLevel,
                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = meshLevel.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      WellControls const & wellControls = getWellControls( subRegion );
      PerforationData const & perforationData = *subRegion.getPerforationData();

      // get the info stored on well elements
      arrayView1d< real64 const > const wellElemGravCoef =
        subRegion.getExtrinsicData< extrinsicMeshData::well::gravityCoefficient >();

      // get well primary variables on well elements
      arrayView1d< real64 > const wellElemPressure =
        subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
      arrayView1d< real64 > const connRate =
        subRegion.getExtrinsicData< extrinsicMeshData::well::connectionRate >();

      // get the element region, subregion, index
      arrayView1d< localIndex const > const resElementRegion =
        perforationData.getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString() );
      arrayView1d< localIndex const > const resElementSubRegion =
        perforationData.getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString() );
      arrayView1d< localIndex const > const resElementIndex =
        perforationData.getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString() );

      // TODO: change the way we access the flowSolver here
      SinglePhaseBase const & flowSolver = getParent().getGroup< SinglePhaseBase >( getFlowSolverName() );
      PresInitializationKernel::SinglePhaseFlowAccessors resSinglePhaseFlowAccessors( meshLevel.getElemManager(), flowSolver.getName() );
      PresInitializationKernel::SingleFluidAccessors resSingleFluidAccessors( meshLevel.getElemManager(), flowSolver.getName() );

      // 1) Loop over all perforations to compute an average density
      // 2) Initialize the reference pressure
      // 3) Estimate the pressures in the well elements using the average density
      PresInitializationKernel::launch( perforationData.size(),
                                        subRegion.size(),
                                        perforationData.getNumPerforationsGlobal(),
                                        wellControls,
                                        0.0, // initialization done at t = 0
                                        resSinglePhaseFlowAccessors.get( extrinsicMeshData::flow::pressure{} ),
                                        resSingleFluidAccessors.get( extrinsicMeshData::singlefluid::density{} ),
                                        resElementRegion,
                                        resElementSubRegion,
                                        resElementIndex,
                                        wellElemGravCoef,
                                        wellElemPressure );

      // 4) Recompute the pressure-dependent properties
      // Note: I am leaving that here because I would like to use the perforationRates (computed in UpdateState)
      //       to better initialize the rates
      updateSubRegionState( meshLevel, subRegion );

      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      SingleFluidBase & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );
      arrayView2d< real64 const > const & wellElemDens = fluid.density();

      // 5) Estimate the well rates
      RateInitializationKernel::launch( subRegion.size(),
                                        wellControls,
                                        0.0, // initialization done at t = 0
                                        wellElemDens,
                                        connRate );
    } );

  } );
}

void SinglePhaseWell::assembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                         real64 const dt,
                                         DomainPartition const & domain,
                                         DofManager const & dofManager,
                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                         arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;


  // loop over the wells
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {
      // get a reference to the degree-of-freedom numbers
      string const wellDofKey = dofManager.getKey( wellElementDofName() );
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );
      arrayView1d< localIndex const > const & nextWellElemIndex =
        subRegion.getReference< array1d< localIndex > >( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString() );

      // get a reference to the primary variables on well elements
      arrayView1d< real64 const > const connRate =
        subRegion.getExtrinsicData< extrinsicMeshData::well::connectionRate >();
      arrayView1d< real64 const > const dConnRate =
        subRegion.getExtrinsicData< extrinsicMeshData::well::deltaConnectionRate >();

      FluxKernel::launch( subRegion.size(),
                          dofManager.rankOffset(),
                          wellElemDofNumber,
                          nextWellElemIndex,
                          connRate,
                          dConnRate,
                          dt,
                          localMatrix,
                          localRhs );
    } );

  } );
}

void SinglePhaseWell::assemblePressureRelations( DomainPartition const & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();


    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {

      WellControls & wellControls = getWellControls( subRegion );

      // get the degrees of freedom numbers, depth, next well elem index
      string const wellDofKey = dofManager.getKey( wellElementDofName() );
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );
      arrayView1d< real64 const > const & wellElemGravCoef =
        subRegion.getExtrinsicData< extrinsicMeshData::well::gravityCoefficient >();
      arrayView1d< localIndex const > const & nextWellElemIndex =
        subRegion.getReference< array1d< localIndex > >( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString() );

      // get primary variables on well elements
      arrayView1d< real64 const > const & wellElemPressure =
        subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
      arrayView1d< real64 const > const & dWellElemPressure =
        subRegion.getExtrinsicData< extrinsicMeshData::well::deltaPressure >();

      // get well constitutive data
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      SingleFluidBase const & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );
      arrayView2d< real64 const > const & wellElemDensity = fluid.density();
      arrayView2d< real64 const > const & dWellElemDensity_dPres = fluid.dDensity_dPressure();

      localIndex const controlHasSwitched =
        PressureRelationKernel::launch( subRegion.size(),
                                        dofManager.rankOffset(),
                                        subRegion.isLocallyOwned(),
                                        subRegion.getTopWellElementIndex(),
                                        wellControls,
                                        m_currentTime + m_currentDt, // controls evaluated with BHP/rate of the end of the time interval
                                        wellElemDofNumber,
                                        wellElemGravCoef,
                                        nextWellElemIndex,
                                        wellElemPressure,
                                        dWellElemPressure,
                                        wellElemDensity,
                                        dWellElemDensity_dPres,
                                        localMatrix,
                                        localRhs );

      if( controlHasSwitched == 1 )
      {
        // Note: if BHP control is not viable, we switch to TOTALVOLRATE
        //       if TOTALVOLRATE is not viable, we switch to BHP

        real64 const timeAtEndOfStep = m_currentTime + m_currentDt;

        if( wellControls.getControl() == WellControls::Control::BHP )
        {
          wellControls.switchToTotalRateControl( wellControls.getTargetTotalRate( timeAtEndOfStep ) );
          GEOSX_LOG_LEVEL_RANK_0( 1, "Control switch for well " << subRegion.getName()
                                                                << " from BHP constraint to rate constraint" );
        }
        else
        {
          wellControls.switchToBHPControl( wellControls.getTargetBHP( timeAtEndOfStep ) );
          GEOSX_LOG_LEVEL_RANK_0( 1, "Control switch for well " << subRegion.getName()
                                                                << " from rate constraint to BHP constraint" );
        }
      }

    } );
  } );
}

void SinglePhaseWell::assembleAccumulationTerms( DomainPartition const & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {
      // for now, we do not want to model storage effects in the wells (unless the well is shut)
      WellControls const & wellControls = getWellControls( subRegion );
      if( wellControls.isWellOpen( m_currentTime + m_currentDt ) )
      {
        return;
      }

      // get a reference to the degree-of-freedom numbers
      string const wellElemDofKey = dofManager.getKey( wellElementDofName() );
      arrayView1d< globalIndex const > const wellElemDofNumber = subRegion.getReference< array1d< globalIndex > >( wellElemDofKey );
      arrayView1d< integer const > const wellElemGhostRank = subRegion.ghostRank();

      arrayView1d< real64 const > const wellElemVolume = subRegion.getElementVolume();

      arrayView1d< real64 const > const wellElemDensityOld =
        subRegion.getExtrinsicData< extrinsicMeshData::well::densityOld >();

      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      SingleFluidBase const & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );
      arrayView2d< real64 const > const wellElemDensity = fluid.density();
      arrayView2d< real64 const > const dWellElemDensity_dPres = fluid.dDensity_dPressure();

      AccumulationKernel::launch( subRegion.size(),
                                  dofManager.rankOffset(),
                                  wellElemDofNumber,
                                  wellElemGhostRank,
                                  wellElemVolume,
                                  wellElemDensity,
                                  dWellElemDensity_dPres,
                                  wellElemDensityOld,
                                  localMatrix,
                                  localRhs );

    } );
  } );

}

void SinglePhaseWell::assembleVolumeBalanceTerms( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                                  DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                                  CRSMatrixView< real64, globalIndex const > const & GEOSX_UNUSED_PARAM( localMatrix ),
                                                  arrayView1d< real64 > const & GEOSX_UNUSED_PARAM( localRhs ) )
{
  // not implemented for single phase flow
}

void SinglePhaseWell::computePerforationRates( MeshLevel const & meshLevel,
                                               WellElementSubRegion & subRegion )
{
  GEOSX_MARK_FUNCTION;

  // if the well is shut, we neglect reservoir-well flow that may occur despite the zero rate
  // therefore, we do not want to compute perforation rates and we simply assume they are zero
  WellControls const & wellControls = getWellControls( subRegion );
  if( !wellControls.isWellOpen( m_currentTime + m_currentDt ) )
  {
    return;
  }

  // get the well data
  PerforationData * const perforationData = subRegion.getPerforationData();

  // get the degrees of freedom and depth
  arrayView1d< real64 const > const wellElemGravCoef =
    subRegion.getExtrinsicData< extrinsicMeshData::well::gravityCoefficient >();

  // get well primary variables on well elements
  arrayView1d< real64 const > const wellElemPressure =
    subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
  arrayView1d< real64 const > const dWellElemPressure =
    subRegion.getExtrinsicData< extrinsicMeshData::well::deltaPressure >();

  // get well constitutive data
  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  SingleFluidBase const & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );  arrayView2d< real64 const > const wellElemDensity = fluid.density();
  arrayView2d< real64 const > const dWellElemDensity_dPres = fluid.dDensity_dPressure();
  arrayView2d< real64 const > const wellElemViscosity = fluid.viscosity();
  arrayView2d< real64 const > const dWellElemViscosity_dPres = fluid.dViscosity_dPressure();

  // get well variables on perforations
  arrayView1d< real64 const > const perfGravCoef =
    perforationData->getExtrinsicData< extrinsicMeshData::well::gravityCoefficient >();
  arrayView1d< localIndex const > const perfWellElemIndex =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::wellElementIndexString() );
  arrayView1d< real64 const > const perfTransmissibility =
    perforationData->getReference< array1d< real64 > >( PerforationData::viewKeyStruct::wellTransmissibilityString() );

  arrayView1d< real64 > const perfRate =
    perforationData->getExtrinsicData< extrinsicMeshData::well::perforationRate >();
  arrayView2d< real64 > const dPerfRate_dPres =
    perforationData->getExtrinsicData< extrinsicMeshData::well::dPerforationRate_dPres >();

  // get the element region, subregion, index
  arrayView1d< localIndex const > const resElementRegion =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString() );
  arrayView1d< localIndex const > const resElementSubRegion =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString() );
  arrayView1d< localIndex const > const resElementIndex =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString() );

  // TODO: change the way we access the flowSolver here
  SinglePhaseBase const & flowSolver = getParent().getGroup< SinglePhaseBase >( getFlowSolverName() );
  PerforationKernel::SinglePhaseFlowAccessors resSinglePhaseFlowAccessors( meshLevel.getElemManager(), flowSolver.getName() );
  PerforationKernel::SingleFluidAccessors resSingleFluidAccessors( meshLevel.getElemManager(), flowSolver.getName() );

  PerforationKernel::launch( perforationData->size(),
                             resSinglePhaseFlowAccessors.get( extrinsicMeshData::flow::pressure{} ),
                             resSinglePhaseFlowAccessors.get( extrinsicMeshData::flow::deltaPressure{} ),
                             resSingleFluidAccessors.get( extrinsicMeshData::singlefluid::density{} ),
                             resSingleFluidAccessors.get( extrinsicMeshData::singlefluid::dDensity_dPressure{} ),
                             resSingleFluidAccessors.get( extrinsicMeshData::singlefluid::viscosity{} ),
                             resSingleFluidAccessors.get( extrinsicMeshData::singlefluid::dViscosity_dPressure{} ),
                             wellElemGravCoef,
                             wellElemPressure,
                             dWellElemPressure,
                             wellElemDensity,
                             dWellElemDensity_dPres,
                             wellElemViscosity,
                             dWellElemViscosity_dPres,
                             perfGravCoef,
                             perfWellElemIndex,
                             perfTransmissibility,
                             resElementRegion,
                             resElementSubRegion,
                             resElementIndex,
                             perfRate,
                             dPerfRate_dPres );
}


real64
SinglePhaseWell::calculateResidualNorm( DomainPartition const & domain,
                                        DofManager const & dofManager,
                                        arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  real64 localResidualNorm = 0;
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {
      string const wellDofKey = dofManager.getKey( wellElementDofName() );
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );
      arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

      arrayView1d< real64 const > const wellElemVolume = subRegion.getElementVolume();
      arrayView1d< real64 const > const wellElemDensityOld =
        subRegion.getExtrinsicData< extrinsicMeshData::well::densityOld >();

      WellControls const & wellControls = getWellControls( subRegion );

      ResidualNormKernel::launch< parallelDevicePolicy<> >( localRhs,
                                                            dofManager.rankOffset(),
                                                            subRegion.isLocallyOwned(),
                                                            subRegion.getTopWellElementIndex(),
                                                            wellControls,
                                                            wellElemDofNumber,
                                                            wellElemGhostRank,
                                                            wellElemVolume,
                                                            wellElemDensityOld,
                                                            m_currentTime + m_currentDt, // residual normalized with rate of the end of the
                                                            // time interval
                                                            m_currentDt,
                                                            &localResidualNorm );

    } );
  } );
  // compute global residual norm
  return sqrt( MpiWrapper::sum( localResidualNorm, MPI_COMM_GEOSX ) );
}

bool SinglePhaseWell::checkSystemSolution( DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           arrayView1d< real64 const > const & localSolution,
                                           real64 const scalingFactor )
{
  GEOSX_MARK_FUNCTION;

  localIndex localCheck = 1;
  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )

    {
      // get the degree of freedom numbers on well elements
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );
      arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

      // get a reference to the primary variables on well elements
      arrayView1d< real64 const > const & wellElemPressure =
        subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
      arrayView1d< real64 const > const & dWellElemPressure =
        subRegion.getExtrinsicData< extrinsicMeshData::well::deltaPressure >();

      // here we can reuse the flow solver kernel checking that pressures are positive
      localIndex const subRegionSolutionCheck =
        singlePhaseWellKernels::
          SolutionCheckKernel::launch< parallelDevicePolicy<> >( localSolution,
                                                                 dofManager.rankOffset(),
                                                                 wellElemDofNumber,
                                                                 wellElemGhostRank,
                                                                 wellElemPressure,
                                                                 dWellElemPressure,
                                                                 scalingFactor );

      if( subRegionSolutionCheck == 0 )
      {
        localCheck = 0;
      }
    } );
  } );

  return MpiWrapper::min( localCheck );
}

void
SinglePhaseWell::applySystemSolution( DofManager const & dofManager,
                                      arrayView1d< real64 const > const & localSolution,
                                      real64 const scalingFactor,
                                      DomainPartition & domain )
{
  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               extrinsicMeshData::well::deltaPressure::key(),
                               scalingFactor,
                               { m_numDofPerWellElement, 0, 1 } );

  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               extrinsicMeshData::well::deltaConnectionRate::key(),
                               scalingFactor,
                               { m_numDofPerWellElement, 1, m_numDofPerWellElement } );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    FieldIdentifiers const fieldsToBeSync{
      SyncFieldsID( FieldLocation::Elem,
                    { extrinsicMeshData::well::deltaPressure::key(), extrinsicMeshData::well::deltaConnectionRate::key() },
                    regionNames )
    };

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );

  } );

  // update properties
  updateState( domain );
}

void SinglePhaseWell::resetStateToBeginningOfStep( DomainPartition & domain )
{

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      // get a reference to the primary variables on well elements
      arrayView1d< real64 > const & dWellElemPressure =
        subRegion.getExtrinsicData< extrinsicMeshData::well::deltaPressure >();
      arrayView1d< real64 > const & dConnRate =
        subRegion.getExtrinsicData< extrinsicMeshData::well::deltaConnectionRate >();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
      {
        dWellElemPressure[iwelem] = 0;
        dConnRate[iwelem] = 0;
      } );
    } );

  } );

  // call constitutive models
  updateState( domain );
}


void SinglePhaseWell::backupFields( MeshLevel & mesh, arrayView1d< string const > const & regionNames ) const
{

  ElementRegionManager & elemManager = mesh.getElemManager();

  elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                            [&]( localIndex const,
                                                                 WellElementSubRegion & subRegion )
  {
    string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
    SingleFluidBase const & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );    arrayView2d< real64 const > const wellElemDensity = fluid.density();

    arrayView1d< real64 > const & wellElemDensityOld =
      subRegion.getExtrinsicData< extrinsicMeshData::well::densityOld >();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      wellElemDensityOld[ei] = wellElemDensity[ei][0];
    } );
  } );
}

void SinglePhaseWell::implicitStepSetup( real64 const & time,
                                         real64 const & dt,
                                         DomainPartition & domain )
{
  WellSolverBase::implicitStepSetup( time, dt, domain );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      validateWellConstraints( subRegion );
    } );
  } );
}


void SinglePhaseWell::implicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time ),
                                            real64 const & GEOSX_UNUSED_PARAM( dt ),
                                            DomainPartition & domain )
{
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      // get a reference to the primary variables on well elements
      arrayView1d< real64 > const wellElemPressure = subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
      arrayView1d< real64 const > const dWellElemPressure = subRegion.getExtrinsicData< extrinsicMeshData::well::deltaPressure >();
      arrayView1d< real64 > const connRate = subRegion.getExtrinsicData< extrinsicMeshData::well::connectionRate >();
      arrayView1d< real64 const > const dConnRate = subRegion.getExtrinsicData< extrinsicMeshData::well::deltaConnectionRate >();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
      {
        wellElemPressure[iwelem] += dWellElemPressure[iwelem];
        connRate[iwelem]         += dConnRate[iwelem];
      } );
    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseWell, string const &, Group * const )
}// namespace geosx
