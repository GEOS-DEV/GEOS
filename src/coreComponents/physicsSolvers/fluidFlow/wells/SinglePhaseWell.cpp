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
 * @file SinglePhaseWell.cpp
 */

#include "SinglePhaseWell.hpp"

#include "mpiCommunications/CommunicationTools.hpp"
#include "dataRepository/Group.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "constitutive/fluid/singleFluidSelector.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "meshUtilities/PerforationData.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace SinglePhaseWellKernels;

SinglePhaseWell::SinglePhaseWell( const string & name,
                                  Group * const parent ):
  WellSolverBase( name, parent )
{
  m_numDofPerWellElement = 2;
}

void SinglePhaseWell::PostProcessInput()
{
  WellSolverBase::PostProcessInput();

  SinglePhaseBase const * const flowSolver = getParent()->GetGroup< SinglePhaseBase >( GetFlowSolverName() );
  GEOSX_ERROR_IF( flowSolver == nullptr,
                  "Flow solver " << GetFlowSolverName() << " not found or incompatible type "
                                                           "(referenced from well solver " << getName() << ")" );
}

void SinglePhaseWell::RegisterDataOnMesh( Group * const meshBodies )
{
  WellSolverBase::RegisterDataOnMesh( meshBodies );

  MeshLevel & meshLevel = *meshBodies->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {
    subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::pressureString )->setPlotLevel( PlotLevel::LEVEL_0 );
    subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaPressureString );
    subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::connRateString )->setPlotLevel( PlotLevel::LEVEL_0 );
    subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaConnRateString );

    PerforationData & perforationData = *subRegion.GetPerforationData();
    perforationData.registerWrapper< array1d< real64 > >( viewKeyStruct::perforationRateString );
    perforationData.registerWrapper< array2d< real64 > >( viewKeyStruct::dPerforationRate_dPresString )->
      reference().resizeDimension< 1 >( 2 );

    WellControls & wellControls = GetWellControls( subRegion );
    wellControls.registerWrapper< real64 >( viewKeyStruct::currentBHPString );
    wellControls.registerWrapper< real64 >( viewKeyStruct::dCurrentBHP_dPresString );

    wellControls.registerWrapper< real64 >( viewKeyStruct::currentVolRateString );
    wellControls.registerWrapper< real64 >( viewKeyStruct::dCurrentVolRate_dPresString );
    wellControls.registerWrapper< real64 >( viewKeyStruct::dCurrentVolRate_dRateString );

  } );
}

void SinglePhaseWell::InitializePreSubGroups( Group * const rootGroup )
{

  WellSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup< DomainPartition >( keys::domain );
  MeshLevel & meshLevel = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  ValidateModelMapping< SingleFluidBase >( *meshLevel.getElemManager(), m_fluidModelNames );
  ValidateWellConstraints( meshLevel );
}

void SinglePhaseWell::ValidateWellConstraints( MeshLevel const & meshLevel ) const
{
  // now that we know we are single-phase, we can check a few things in the constraints
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {
    WellControls const & wellControls = GetWellControls( subRegion );
    WellControls::Control const currentControl = wellControls.GetControl();
    real64 const targetTotalRate = wellControls.GetTargetTotalRate();
    real64 const targetPhaseRate = wellControls.GetTargetPhaseRate();
    GEOSX_ERROR_IF( currentControl == WellControls::Control::PHASEVOLRATE,
                    "Phase rate control is not available for SinglePhaseWell" );
    GEOSX_ERROR_IF( isZero( targetTotalRate ),
                    "Target total rate cannot be equal to zero" );
    GEOSX_ERROR_IF( !isZero( targetPhaseRate ),
                    "Target phase rate cannot be used for SinglePhaseWell" );
  } );
}

void SinglePhaseWell::UpdateBHPForConstraint( WellElementSubRegion & subRegion, localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.IsLocallyOwned() )
  {
    return;
  }

  localIndex const iwelemRef = subRegion.GetTopWellElementIndex();

  // subRegion data

  arrayView1d< real64 const > const pres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const dPres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  arrayView1d< real64 const > const wellElemGravCoef =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

  // fluid data

  SingleFluidBase & fluid = GetConstitutiveModel< SingleFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
  arrayView2d< real64 const > const & dens = fluid.density();
  arrayView2d< real64 const > const & dDens_dPres = fluid.dDensity_dPressure();

  // control data

  WellControls & wellControls = GetWellControls( subRegion );

  real64 const & refGravCoef = wellControls.GetReferenceGravityCoef();

  real64 & currentBHP =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentBHPString );
  real64 & dCurrentBHP_dPres =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::dCurrentBHP_dPresString );

  // bring everything back to host, capture the scalars by reference
  forAll< serialPolicy >( 1, [pres,
                              dPres,
                              dens,
                              dDens_dPres,
                              wellElemGravCoef,
                              &currentBHP,
                              &dCurrentBHP_dPres,
                              &iwelemRef,
                              &refGravCoef] GEOSX_HOST_DEVICE ( localIndex const )
  {
    currentBHP = pres[iwelemRef] + dPres[iwelemRef] + dens[iwelemRef][0] * ( refGravCoef - wellElemGravCoef[iwelemRef] );
    dCurrentBHP_dPres = 1.0 + dDens_dPres[iwelemRef][0] * ( refGravCoef - wellElemGravCoef[iwelemRef] );
  } );
}

void SinglePhaseWell::UpdateVolRateForConstraint( WellElementSubRegion & subRegion, localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.IsLocallyOwned() )
  {
    return;
  }

  localIndex const iwelemRef = subRegion.GetTopWellElementIndex();

  // subRegion data

  arrayView1d< real64 const > const pres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const dPres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  arrayView1d< real64 const > const & connRate =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::connRateString );
  arrayView1d< real64 const > const & dConnRate =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaConnRateString );

  // fluid data

  SingleFluidBase & fluid = GetConstitutiveModel< SingleFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
  arrayView2d< real64 const > const & dens = fluid.density();
  arrayView2d< real64 const > const & dDens_dPres = fluid.dDensity_dPressure();

  // control data

  WellControls & wellControls = GetWellControls( subRegion );

  integer const useSurfaceConditions = wellControls.UseSurfaceConditions();

  real64 & currentVolRate =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentVolRateString );
  real64 & dCurrentVolRate_dPres =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::dCurrentVolRate_dPresString );
  real64 & dCurrentVolRate_dRate =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::dCurrentVolRate_dRateString );

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
                                &currentVolRate,
                                &dCurrentVolRate_dPres,
                                &dCurrentVolRate_dRate,
                                &iwelemRef] GEOSX_HOST_DEVICE ( localIndex const )
    {
      //    We need to evaluate the density as follows:
      //      - Surface conditions: using the surface pressure provided by the user
      //      - Reservoir conditions: using the pressure in the top element

      if( useSurfaceConditions )
      {
        real64 const surfacePres = 101325.0;
        // we need to compute the surface density
        fluidWrapper.Update( iwelemRef, 0, surfacePres );
      }
      else
      {
        real64 const refPres = pres[iwelemRef] + dPres[iwelemRef];
        fluidWrapper.Update( iwelemRef, 0, refPres );
      }

      real64 const densInv = 1.0 / dens[iwelemRef][0];
      currentVolRate = ( connRate[iwelemRef] + dConnRate[iwelemRef] ) * densInv;
      dCurrentVolRate_dPres = -( useSurfaceConditions ==  0 ) * dDens_dPres[iwelemRef][0] * currentVolRate * densInv;
      dCurrentVolRate_dRate = densInv;

    } );
  } );
}

void SinglePhaseWell::UpdateFluidModel( WellElementSubRegion & subRegion, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const pres = subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const dPres = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  SingleFluidBase & fluid = GetConstitutiveModel< SingleFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

  constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    SinglePhaseBaseKernels::FluidUpdateKernel::Launch( fluidWrapper, pres, dPres );
  } );
}

void SinglePhaseWell::UpdateState( WellElementSubRegion & subRegion, localIndex const targetIndex )
{
  // update volumetric rates for the well constraints
  // Warning! This must be called before updating the fluid model
  UpdateVolRateForConstraint( subRegion, targetIndex );

  // update density in the well elements
  UpdateFluidModel( subRegion, targetIndex );

  // update the current BHP
  UpdateBHPForConstraint( subRegion, targetIndex );

  // update perforation rates
  ComputePerforationRates( subRegion, targetIndex );
}

void SinglePhaseWell::InitializeWells( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
                                                               WellElementSubRegion & subRegion )
  {
    WellControls const & wellControls = GetWellControls( subRegion );
    PerforationData const & perforationData = *subRegion.GetPerforationData();

    // get the info stored on well elements
    arrayView1d< real64 const > const wellElemGravCoef =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

    // get well primary variables on well elements
    arrayView1d< real64 > const wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 > const connRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::connRateString );

    // get the element region, subregion, index
    arrayView1d< localIndex const > const resElementRegion =
      perforationData.getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString );
    arrayView1d< localIndex const > const resElementSubRegion =
      perforationData.getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString );
    arrayView1d< localIndex const > const resElementIndex =
      perforationData.getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString );

    // 1) Loop over all perforations to compute an average density
    // 2) Initialize the reference pressure
    // 3) Estimate the pressures in the well elements using the average density
    PresInitializationKernel::Launch< parallelDevicePolicy<> >( perforationData.size(),
                                                                subRegion.size(),
                                                                perforationData.GetNumPerforationsGlobal(),
                                                                wellControls,
                                                                m_resPressure.toNestedViewConst(),
                                                                m_resDensity.toNestedViewConst(),
                                                                resElementRegion,
                                                                resElementSubRegion,
                                                                resElementIndex,
                                                                wellElemGravCoef,
                                                                wellElemPressure );

    // 4) Recompute the pressure-dependent properties
    // Note: I am leaving that here because I would like to use the perforationRates (computed in UpdateState)
    //       to better initialize the rates
    UpdateState( subRegion, targetIndex );

    SingleFluidBase const & fluid = GetConstitutiveModel< SingleFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
    arrayView2d< real64 const > const & wellElemDens = fluid.density();

    // 5) Estimate the well rates
    RateInitializationKernel::Launch< parallelDevicePolicy<> >( subRegion.size(),
                                                                wellControls,
                                                                wellElemDens,
                                                                connRate );
  } );
}

void SinglePhaseWell::AssembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                         real64 const dt,
                                         DomainPartition const & domain,
                                         DofManager const & dofManager,
                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                         arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  // saved current dt for residual normalization
  m_currentDt = dt;

  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {
    // get a reference to the degree-of-freedom numbers
    string const wellDofKey = dofManager.getKey( WellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< localIndex const > const & nextWellElemIndex =
      subRegion.getReference< array1d< localIndex > >( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );

    // get a reference to the primary variables on well elements
    arrayView1d< real64 const > const & connRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::connRateString );
    arrayView1d< real64 const > const & dConnRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaConnRateString );

    FluxKernel::Launch< parallelDevicePolicy<> >( subRegion.size(),
                                                  dofManager.rankOffset(),
                                                  wellElemDofNumber,
                                                  nextWellElemIndex,
                                                  connRate,
                                                  dConnRate,
                                                  dt,
                                                  localMatrix,
                                                  localRhs );
  } );
}

void SinglePhaseWell::FormPressureRelations( DomainPartition const & domain,
                                             DofManager const & dofManager,
                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                             arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
                                                               WellElementSubRegion const & subRegion )
  {

    WellControls & wellControls = GetWellControls( subRegion );

    // get the degrees of freedom numbers, depth, next well elem index
    string const wellDofKey = dofManager.getKey( WellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< real64 const > const & wellElemGravCoef =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );
    arrayView1d< localIndex const > const & nextWellElemIndex =
      subRegion.getReference< array1d< localIndex > >( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );

    // get primary variables on well elements
    arrayView1d< real64 const > const & wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 const > const & dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

    // get well constitutive data
    SingleFluidBase const & fluid = GetConstitutiveModel< SingleFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
    arrayView2d< real64 const > const & wellElemDensity = fluid.density();
    arrayView2d< real64 const > const & dWellElemDensity_dPres = fluid.dDensity_dPressure();

    localIndex const controlHasSwitched =
      PressureRelationKernel::Launch< parallelDevicePolicy<>,
                                      parallelDeviceReduce >( subRegion.size(),
                                                              dofManager.rankOffset(),
                                                              subRegion.IsLocallyOwned(),
                                                              subRegion.GetTopWellElementIndex(),
                                                              wellControls,
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

      if( wellControls.GetControl() == WellControls::Control::BHP )
      {
        wellControls.SwitchToTotalRateControl( wellControls.GetTargetTotalRate() );
        GEOSX_LOG_LEVEL_RANK_0( 1, "Control switch for well " << subRegion.getName()
                                                              << " from BHP constraint to rate constraint" );
      }
      else
      {
        wellControls.SwitchToBHPControl( wellControls.GetTargetBHP() );
        GEOSX_LOG_LEVEL_RANK_0( 1, "Control switch for well " << subRegion.getName()
                                                              << " from rate constraint to BHP constraint" );
      }
    }

  } );
}

void SinglePhaseWell::AssembleVolumeBalanceTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                  real64 const GEOSX_UNUSED_PARAM( dt ),
                                                  DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                                  DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                                  CRSMatrixView< real64, globalIndex const > const & GEOSX_UNUSED_PARAM( localMatrix ),
                                                  arrayView1d< real64 > const & GEOSX_UNUSED_PARAM( localRhs ) )
{
  // not implemented for single phase flow
}

void SinglePhaseWell::ComputePerforationRates( WellElementSubRegion & subRegion, localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  // get the well data
  PerforationData * const perforationData = subRegion.GetPerforationData();

  // get the degrees of freedom and depth
  arrayView1d< real64 const > const wellElemGravCoef =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

  // get well primary variables on well elements
  arrayView1d< real64 const > const
  wellElemPressure = subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const
  dWellElemPressure = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  // get well constitutive data
  SingleFluidBase const & fluid = GetConstitutiveModel< SingleFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
  arrayView2d< real64 const > const wellElemDensity = fluid.density();
  arrayView2d< real64 const > const dWellElemDensity_dPres = fluid.dDensity_dPressure();
  arrayView2d< real64 const > const wellElemViscosity = fluid.viscosity();
  arrayView2d< real64 const > const dWellElemViscosity_dPres = fluid.dViscosity_dPressure();

  // get well variables on perforations
  arrayView1d< real64 const > const perfGravCoef =
    perforationData->getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );
  arrayView1d< localIndex const > const perfWellElemIndex =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::wellElementIndexString );
  arrayView1d< real64 const > const perfTransmissibility =
    perforationData->getReference< array1d< real64 > >( PerforationData::viewKeyStruct::wellTransmissibilityString );

  arrayView1d< real64 > const perfRate =
    perforationData->getReference< array1d< real64 > >( viewKeyStruct::perforationRateString );
  arrayView2d< real64 > const dPerfRate_dPres =
    perforationData->getReference< array2d< real64 > >( viewKeyStruct::dPerforationRate_dPresString );

  // get the element region, subregion, index
  arrayView1d< localIndex const > const resElementRegion =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString );
  arrayView1d< localIndex const > const resElementSubRegion =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString );
  arrayView1d< localIndex const > const resElementIndex =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString );

  PerforationKernel::Launch< parallelDevicePolicy<> >( perforationData->size(),
                                                       m_resPressure.toNestedViewConst(),
                                                       m_deltaResPressure.toNestedViewConst(),
                                                       m_resDensity.toNestedViewConst(),
                                                       m_dResDens_dPres.toNestedViewConst(),
                                                       m_resViscosity.toNestedViewConst(),
                                                       m_dResVisc_dPres.toNestedViewConst(),
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
SinglePhaseWell::CalculateResidualNorm( DomainPartition const & domain,
                                        DofManager const & dofManager,
                                        arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  real64 localResidualNorm = 0;
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
                                                               WellElementSubRegion const & subRegion )
  {
    string const wellDofKey = dofManager.getKey( WellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

    SingleFluidBase const & fluid = GetConstitutiveModel< SingleFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
    arrayView2d< real64 const > const & wellElemDensity = fluid.density();

    WellControls const & wellControls = GetWellControls( subRegion );

    ResidualNormKernel::Launch< parallelDevicePolicy<>,
                                parallelDeviceReduce >( localRhs,
                                                        dofManager.rankOffset(),
                                                        subRegion.IsLocallyOwned(),
                                                        subRegion.GetTopWellElementIndex(),
                                                        wellControls,
                                                        wellElemDofNumber,
                                                        wellElemGhostRank,
                                                        wellElemDensity,
                                                        m_currentDt,
                                                        &localResidualNorm );

  } );

  // compute global residual norm
  return sqrt( MpiWrapper::Sum( localResidualNorm, MPI_COMM_GEOSX ) );
}

bool SinglePhaseWell::CheckSystemSolution( DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           arrayView1d< real64 const > const & localSolution,
                                           real64 const scalingFactor )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  localIndex localCheck = 1;

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {
    // get the degree of freedom numbers on well elements
    string const wellDofKey = dofManager.getKey( WellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

    // get a reference to the primary variables on well elements
    arrayView1d< real64 const > const & wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 const > const & dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

    // here we can reuse the flow solver kernel checking that pressures are positive
    localIndex const subRegionSolutionCheck =
      SinglePhaseWellKernels::SolutionCheckKernel::Launch< parallelDevicePolicy<>,
                                                           parallelDeviceReduce >( localSolution,
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

  return MpiWrapper::Min( localCheck );
}

void
SinglePhaseWell::ApplySystemSolution( DofManager const & dofManager,
                                      arrayView1d< real64 const > const & localSolution,
                                      real64 const scalingFactor,
                                      DomainPartition & domain )
{
  dofManager.addVectorToField( localSolution,
                               WellElementDofName(),
                               viewKeyStruct::deltaPressureString,
                               scalingFactor,
                               0, 1 );

  dofManager.addVectorToField( localSolution,
                               WellElementDofName(),
                               viewKeyStruct::deltaConnRateString,
                               scalingFactor,
                               1, m_numDofPerWellElement );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaPressureString ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaConnRateString ) );
  CommunicationTools::SynchronizeFields( fieldNames,
                                         domain.getMeshBody( 0 )->getMeshLevel( 0 ),
                                         domain.getNeighbors(),
                                         true );

  // update properties
  UpdateStateAll( domain );
}

void SinglePhaseWell::ResetStateToBeginningOfStep( DomainPartition & domain )
{

  MeshLevel & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {
    // get a reference to the primary variables on well elements
    arrayView1d< real64 > const & dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );
    arrayView1d< real64 > const & dConnRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaConnRateString );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      dWellElemPressure[iwelem] = 0;
      dConnRate[iwelem] = 0;
    } );
  } );

  // call constitutive models
  UpdateStateAll( domain );
}


void SinglePhaseWell::ResetViews( DomainPartition & domain )
{
  WellSolverBase::ResetViews( domain );

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager & elemManager = *mesh.getElemManager();

  SinglePhaseBase & flowSolver = *getParent()->GetGroup< SinglePhaseBase >( GetFlowSolverName() );

  {
    using keys = SinglePhaseBase::viewKeyStruct;

    m_resPressure.clear();
    m_resPressure = elemManager.ConstructArrayViewAccessor< real64, 1 >( keys::pressureString );
    m_resPressure.setName( getName() + "/accessors/" + keys::pressureString );

    m_deltaResPressure.clear();
    m_deltaResPressure = elemManager.ConstructArrayViewAccessor< real64, 1 >( keys::deltaPressureString );
    m_deltaResPressure.setName( getName() + "/accessors/" + keys::deltaPressureString );

  }
  {
    using keys = SingleFluidBase::viewKeyStruct;

    m_resDensity.clear();
    m_resDensity = elemManager.ConstructMaterialArrayViewAccessor< real64, 2 >( keys::densityString,
                                                                                flowSolver.targetRegionNames(),
                                                                                flowSolver.fluidModelNames() );
    m_resDensity.setName( getName() + "/accessors/" + keys::densityString );

    m_dResDens_dPres.clear();
    m_dResDens_dPres = elemManager.ConstructMaterialArrayViewAccessor< real64, 2 >( keys::dDens_dPresString,
                                                                                    flowSolver.targetRegionNames(),
                                                                                    flowSolver.fluidModelNames() );
    m_dResDens_dPres.setName( getName() + "/accessors/" + keys::dDens_dPresString );

    m_resViscosity.clear();
    m_resViscosity = elemManager.ConstructMaterialArrayViewAccessor< real64, 2 >( keys::viscosityString,
                                                                                  flowSolver.targetRegionNames(),
                                                                                  flowSolver.fluidModelNames() );
    m_resViscosity.setName( getName() + "/accessors/" + keys::viscosityString );

    m_dResVisc_dPres.clear();
    m_dResVisc_dPres = elemManager.ConstructMaterialArrayViewAccessor< real64, 2 >( keys::dVisc_dPresString,
                                                                                    flowSolver.targetRegionNames(),
                                                                                    flowSolver.fluidModelNames() );
    m_dResVisc_dPres.setName( getName() + "/accessors/" + keys::dVisc_dPresString );

  }
}

void SinglePhaseWell::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time ),
                                            real64 const & GEOSX_UNUSED_PARAM( real64 const & dt ),
                                            DomainPartition & domain )
{
  MeshLevel & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager & elemManager = *meshLevel.getElemManager();

  elemManager.forElementSubRegions< WellElementSubRegion >( [&]( WellElementSubRegion & subRegion )
  {
    // get a reference to the primary variables on well elements
    arrayView1d< real64 > const wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 const > const dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );
    arrayView1d< real64 > const connRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::connRateString );
    arrayView1d< real64 const > const dConnRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaConnRateString );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      wellElemPressure[iwelem] += dWellElemPressure[iwelem];
      connRate[iwelem]         += dConnRate[iwelem];
    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseWell, string const &, Group * const )
}// namespace geosx
