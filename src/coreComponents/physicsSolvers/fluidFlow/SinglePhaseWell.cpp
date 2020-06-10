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
#include "wells/PerforationData.hpp"
#include "wells/WellElementSubRegion.hpp"
#include "wells/WellControls.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseWellKernels.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace SinglePhaseWellKernels;

SinglePhaseWell::SinglePhaseWell( const string & name,
                                  Group * const parent )
  :
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

    PerforationData * const perforationData = subRegion.GetPerforationData();
    perforationData->registerWrapper< array1d< real64 > >( viewKeyStruct::perforationRateString );
    perforationData->registerWrapper< array2d< real64 > >( viewKeyStruct::dPerforationRate_dPresString );
  } );
}

void SinglePhaseWell::InitializePreSubGroups( Group * const rootGroup )
{

  WellSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup< DomainPartition >( keys::domain );
  MeshLevel & meshLevel = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  ValidateModelMapping< SingleFluidBase >( *meshLevel.getElemManager(), m_fluidModelNames );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {
    PerforationData & perforationData = *subRegion.GetPerforationData();
    perforationData.getReference< array2d< real64 > >( viewKeyStruct::dPerforationRate_dPresString ).resizeDimension< 1 >( 2 );
  } );
}

void SinglePhaseWell::UpdateFluidModel( WellElementSubRegion & subRegion, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const & pres = subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dPres = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  SingleFluidBase & fluid = GetConstitutiveModel< SingleFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

  bool const success =
    constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    SinglePhaseBaseKernels::FluidUpdateKernel::Launch( fluidWrapper, pres, dPres );
  } );
  GEOSX_ERROR_IF( !success, "Kernel not launched due to unknown fluid type" );
}

void SinglePhaseWell::UpdateState( WellElementSubRegion & subRegion, localIndex const targetIndex )
{
  // update density in the well elements
  UpdateFluidModel( subRegion, targetIndex );

  // update perforation rates
  ComputePerforationRates( subRegion, targetIndex );
}

void SinglePhaseWell::InitializeWells( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const & resPressure = m_resPressure;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const & resDensity = m_resDensity;

  MeshLevel & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
                                                               WellElementSubRegion & subRegion )
  {
    WellControls const & wellControls = GetWellControls( subRegion );
    PerforationData const * const perforationData = subRegion.GetPerforationData();

    // get the info stored on well elements
    arrayView1d< real64 const > const & wellElemGravCoef =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

    // get well primary variables on well elements
    arrayView1d< real64 > const & wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 > const & connRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::connRateString );

    // get the element region, subregion, index
    arrayView1d< localIndex const > const & resElementRegion =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString );
    arrayView1d< localIndex const > const & resElementSubRegion =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString );
    arrayView1d< localIndex const > const & resElementIndex =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString );

    // define a reservoir pressure used for initialization
    real64 resPres = ( wellControls.GetType() == WellControls::Type::PRODUCER )
                     ? 1e20 : 0;

    // 1) Loop over all perforations to compute an average density
    real64 avgDensity = 0;
    for( localIndex iperf = 0; iperf < perforationData->size(); ++iperf )
    {
      // get the reservoir (sub)region and element indices
      localIndex const er = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei = resElementIndex[iperf];

      avgDensity += resDensity[er][esr][ei][0];

      // save min pressure for producer
      if( wellControls.GetType() == WellControls::Type::PRODUCER &&
          resPres > resPressure[er][esr][ei] )
      {
        resPres = resPressure[er][esr][ei];
      }
      // save max pressure for injector
      else if( wellControls.GetType() == WellControls::Type::INJECTOR &&
               resPres < resPressure[er][esr][ei] )
      {
        resPres = resPressure[er][esr][ei];
      }
    }

    // communicate the pressures to the ranks without perforations
    // this will be used to initialize the pressure, starting by the owner rank
    if( wellControls.GetType() == WellControls::Type::PRODUCER )
    {
      resPres = MpiWrapper::Min( resPres );
    }
    else if( wellControls.GetType() == WellControls::Type::INJECTOR )
    {
      resPres = MpiWrapper::Max( resPres );
    }

    avgDensity = MpiWrapper::Sum( avgDensity );

    globalIndex const numPerforationsGlobal = perforationData->GetNumPerforationsGlobal();
    avgDensity /= numPerforationsGlobal;

    real64 pressureControl = 0.0;
    real64 gravCoefControl = 0.0;
    if( subRegion.IsLocallyOwned() )
    {

      // get the reference data for this well
      localIndex const iwelemControl = wellControls.GetReferenceWellElementIndex();
      gravCoefControl = wellElemGravCoef[iwelemControl];

      // 2) Initialize the reference pressure
      real64 const & targetBHP = wellControls.GetTargetBHP();
      if( wellControls.GetControl() == WellControls::Control::BHP )
      {
        // if pressure constraint, set the ref pressure at the constraint
        pressureControl = targetBHP;
      }
      else // rate control
      {
        // if rate constraint, set the ref pressure slightly
        // above/below the target pressure depending on well type
        pressureControl =
          ( wellControls.GetType() == WellControls::Type::PRODUCER )
          ? 0.5 * resPres // hard-coded values come from personal communication with Hui
          : 2.0 * resPres;
      }

      wellElemPressure[iwelemControl] = pressureControl;
    }

    // TODO optimize
    MpiWrapper::Broadcast( pressureControl, subRegion.GetTopRank() );
    MpiWrapper::Broadcast( gravCoefControl, subRegion.GetTopRank() );

    GEOSX_ERROR_IF( pressureControl <= 0, "Invalid well initialization: negative pressure was found" );

    // 3) Estimate the pressures in the well elements using this avgDensity
    forAll< serialPolicy >( subRegion.size(), [=]( localIndex const iwelem )
    {
      wellElemPressure[iwelem] = pressureControl
                                 + avgDensity * ( wellElemGravCoef[iwelem] - gravCoefControl );
    } );

    // 4) Recompute the pressure-dependent properties
    UpdateState( subRegion, targetIndex );

    real64 const targetRate = wellControls.GetTargetRate();
    WellControls::Control const control = wellControls.GetControl();
    WellControls::Type const wellType = wellControls.GetType();

    // 5) Estimate the connection rates based on the min/max pressure
    forAll< serialPolicy >( subRegion.size(), [=]( localIndex const iwelem )
    {
      if( control == WellControls::Control::BHP )
      {
        // if BHP constraint set rate below the absolute max rate
        // with the appropriate sign (negative for prod, positive for inj)
        connRate[iwelem] = ( wellType == WellControls::Type::PRODUCER )
                           ? std::max( 0.1 * targetRate, -1e3 )  // hard-coded values come from personal communication
                                                                 // with Hui
                           : std::min( 0.1 * targetRate, 1e3 );
      }
      else
      {
        connRate[iwelem] = targetRate;
      }
    } );
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

    FluxKernel::Launch< serialPolicy >( subRegion.size(),
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


void SinglePhaseWell::FormControlEquation( DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                           arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {

    if( !subRegion.IsLocallyOwned() )
    {
      return;
    }

    WellControls const & wellControls = GetWellControls( subRegion );

    // get the degrees of freedom
    string const wellDofKey = dofManager.getKey( WellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );

    // get primary variables on well elements
    arrayView1d< real64 const > const & wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 const > const & dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );
    arrayView1d< real64 const > const & connRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::connRateString );
    arrayView1d< real64 const > const & dConnRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaConnRateString );

    // get the index of the element where the control is enforced
    localIndex const iwelemControl = wellControls.GetReferenceWellElementIndex();

    ControlEquationHelper::ComputeJacobianEntry( dofManager.rankOffset(),
                                                 wellControls,
                                                 wellElemDofNumber[iwelemControl],
                                                 wellElemPressure[iwelemControl],
                                                 dWellElemPressure[iwelemControl],
                                                 connRate[iwelemControl],
                                                 dConnRate[iwelemControl],
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

    WellControls const & wellControls = GetWellControls( subRegion );

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

    PressureRelationKernel::Launch< serialPolicy >( subRegion.size(),
                                                    dofManager.rankOffset(),
                                                    wellElemDofNumber,
                                                    wellElemGravCoef,
                                                    nextWellElemIndex,
                                                    wellElemPressure,
                                                    dWellElemPressure,
                                                    wellElemDensity,
                                                    dWellElemDensity_dPres,
                                                    wellControls.GetTargetBHP(),
                                                    localMatrix,
                                                    localRhs );
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


void SinglePhaseWell::CheckWellControlSwitch( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {

    if( !subRegion.IsLocallyOwned() )
    {
      return;
    }

    WellControls & wellControls = GetWellControls( subRegion );

    // get the primary variables
    arrayView1d< real64 const > const & wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 const > const & dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );
    arrayView1d< real64 const > const & connRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::connRateString );
    arrayView1d< real64 const > const & dConnRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaConnRateString );

    // if isViable is true at the end of the following checks, no need to switch
    bool controlIsViable = false;

    // get well control and type
    WellControls::Control const currentControl = wellControls.GetControl();
    WellControls::Type const type = wellControls.GetType();

    // again we assume here that the first well element is on this MPI rank
    localIndex const iwelemControl = wellControls.GetReferenceWellElementIndex();

    real64 const refRate = connRate[iwelemControl] + dConnRate[iwelemControl];
    real64 const refPressure = wellElemPressure[iwelemControl] + dWellElemPressure[iwelemControl];

    // BHP control
    if( currentControl == WellControls::Control::BHP )
    {
      // the control is viable if the reference rate is below the max rate
      real64 const & maxRate = wellControls.GetTargetRate();
      controlIsViable = ( fabs( refRate ) <= fabs( maxRate ) );
    }
    else // rate control
    {
      // the control is viable if the reference pressure is below/above the max/min pressure
      if( type == WellControls::Type::PRODUCER )
      {
        // targetBHP specifies a min pressure here
        real64 const & minPressure = wellControls.GetTargetBHP();
        controlIsViable = ( refPressure >= minPressure );
      }
      else
      {
        // targetBHP specifies a max pressure here
        real64 const & maxPressure = wellControls.GetTargetBHP();
        controlIsViable = ( refPressure <= maxPressure );
      }
    }

    if( !controlIsViable )
    {
      if( currentControl == WellControls::Control::BHP )
      {
        wellControls.SetControl( WellControls::Control::LIQUIDRATE,
                                 wellControls.GetTargetRate() );
        GEOSX_LOG_LEVEL_RANK_0( 1, "Control switch for well " << subRegion.getName()
                                                              << " from BHP constraint to rate constraint" );
      }
      else // rate control
      {
        wellControls.SetControl( WellControls::Control::BHP,
                                 wellControls.GetTargetBHP() );

        // Debug information for logLevel >= 1
        GEOSX_LOG_LEVEL_RANK_0( 1, "Control switch for well " << subRegion.getName()
                                                              << " from rate constraint to BHP constraint" );
      }
    }
  } );
}


void SinglePhaseWell::ComputePerforationRates( WellElementSubRegion & subRegion, localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  // get the reservoir data
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > >  const & resPressure        = m_resPressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > >  const & dResPressure       = m_deltaResPressure;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const & resDensity          = m_resDensity;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const & dResDensity_dPres   = m_dResDens_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const & resViscosity        = m_resViscosity;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const & dResViscosity_dPres = m_dResVisc_dPres;

  // get the well data
  PerforationData * const perforationData = subRegion.GetPerforationData();

  // get the degrees of freedom and depth
  arrayView1d< real64 const > const & wellElemGravCoef =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

  // get well primary variables on well elements
  arrayView1d< real64 const > const &
  wellElemPressure = subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const &
  dWellElemPressure = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  // get well constitutive data
  SingleFluidBase const & fluid = GetConstitutiveModel< SingleFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
  arrayView2d< real64 const > const & wellElemDensity = fluid.density();
  arrayView2d< real64 const > const & dWellElemDensity_dPres = fluid.dDensity_dPressure();
  arrayView2d< real64 const > const & wellElemViscosity = fluid.viscosity();
  arrayView2d< real64 const > const & dWellElemViscosity_dPres = fluid.dViscosity_dPressure();

  // get well variables on perforations
  arrayView1d< real64 const > const & perfGravCoef =
    perforationData->getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );
  arrayView1d< localIndex const > const & perfWellElemIndex =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::wellElementIndexString );
  arrayView1d< real64 const > const & perfTransmissibility =
    perforationData->getReference< array1d< real64 > >( PerforationData::viewKeyStruct::wellTransmissibilityString );

  arrayView1d< real64 > const & perfRate =
    perforationData->getReference< array1d< real64 > >( viewKeyStruct::perforationRateString );
  arrayView2d< real64 > const & dPerfRate_dPres =
    perforationData->getReference< array2d< real64 > >( viewKeyStruct::dPerforationRate_dPresString );

  // get the element region, subregion, index
  arrayView1d< localIndex const > const & resElementRegion =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString );
  arrayView1d< localIndex const > const & resElementSubRegion =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString );
  arrayView1d< localIndex const > const & resElementIndex =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString );

  PerforationKernel::Launch< serialPolicy >( perforationData->size(),
                                             resPressure.toViewConst(),
                                             dResPressure.toViewConst(),
                                             resDensity.toViewConst(),
                                             dResDensity_dPres.toViewConst(),
                                             resViscosity.toViewConst(),
                                             dResViscosity_dPres.toViewConst(),
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
    arrayView1d< integer const > const & wellElemGhostRank =
      subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );
    arrayView1d< real64 const > const & wellElemVolume = subRegion.getElementVolume();

    SingleFluidBase const & fluid = GetConstitutiveModel< SingleFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
    arrayView2d< real64 const > const & wellElemDensity = fluid.density();

    ResidualNormKernel::Launch< serialPolicy, serialReduce >( localRhs,
                                                              dofManager.rankOffset(),
                                                              wellElemDofNumber,
                                                              wellElemGhostRank,
                                                              wellElemVolume,
                                                              wellElemDensity,
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
    arrayView1d< integer const > const & wellElemGhostRank =
      subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );

    // get a reference to the primary variables on well elements
    arrayView1d< real64 const > const & wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 const > const & dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

    // here we can reuse the flow solver kernel checking that pressures are positive
    localIndex const subRegionSolutionCheck =
      SinglePhaseWellKernels::SolutionCheckKernel::Launch< serialPolicy,
                                                           serialReduce >( localSolution,
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
  fieldNames["elems"].push_back( viewKeyStruct::deltaPressureString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaConnRateString );
  CommunicationTools::SynchronizeFields( fieldNames,
                                         domain.getMeshBody( 0 )->getMeshLevel( 0 ),
                                         domain.getNeighbors() );

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

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex const iwelem )
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

    m_resPressure =
      elemManager.ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( keys::pressureString );
    m_deltaResPressure =
      elemManager.ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( keys::deltaPressureString );
  }
  {
    using keys = SingleFluidBase::viewKeyStruct;

    m_resDensity =
      elemManager.ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::densityString,
                                                                                             flowSolver.targetRegionNames(),
                                                                                             flowSolver.fluidModelNames() );
    m_dResDens_dPres =
      elemManager.ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::dDens_dPresString,
                                                                                             flowSolver.targetRegionNames(),
                                                                                             flowSolver.fluidModelNames() );
    m_resViscosity =
      elemManager.ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::viscosityString,
                                                                                             flowSolver.targetRegionNames(),
                                                                                             flowSolver.fluidModelNames() );
    m_dResVisc_dPres =
      elemManager.ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::dVisc_dPresString,
                                                                                             flowSolver.targetRegionNames(),
                                                                                             flowSolver.fluidModelNames() );
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
    arrayView1d< real64 > const & wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 const > const & dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );
    arrayView1d< real64 > const & connRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::connRateString );
    arrayView1d< real64 const > const & dConnRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaConnRateString );

    forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const iwelem )
    {
      wellElemPressure[iwelem] += dWellElemPressure[iwelem];
      connRate[iwelem]         += dConnRate[iwelem];
    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseWell, string const &, Group * const )
}// namespace geosx
