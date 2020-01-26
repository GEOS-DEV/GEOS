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
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "managers/DomainPartition.hpp"
#include "wells/PerforationData.hpp"
#include "wells/WellElementSubRegion.hpp"
#include "wells/WellControls.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

SinglePhaseWell::SinglePhaseWell( const string & name,
                                  Group * const parent )
  :
  WellSolverBase( name, parent )
{
  m_numDofPerWellElement = 2;
}

void SinglePhaseWell::RegisterDataOnMesh(Group * const meshBodies)
{

  WellSolverBase::RegisterDataOnMesh(meshBodies);

  MeshLevel * const meshLevel = meshBodies->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  // loop over the wells
  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion * const subRegion )
  {
    subRegion->registerWrapper<array1d<real64>>( viewKeyStruct::pressureString )->setPlotLevel(PlotLevel::LEVEL_0);
    subRegion->registerWrapper<array1d<real64>>( viewKeyStruct::deltaPressureString );
    subRegion->registerWrapper<array1d<real64>>( viewKeyStruct::connRateString )->setPlotLevel(PlotLevel::LEVEL_0);
    subRegion->registerWrapper<array1d<real64>>( viewKeyStruct::deltaConnRateString );

    PerforationData * const perforationData = subRegion->GetPerforationData();
    perforationData->registerWrapper<array1d<real64>>( viewKeyStruct::perforationRateString );
    perforationData->registerWrapper<array2d<real64>>( viewKeyStruct::dPerforationRate_dPresString );
  });

}
  
void SinglePhaseWell::InitializePreSubGroups( Group * const rootGroup )
{

  WellSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup<DomainPartition>( keys::domain );

  MeshLevel * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  // loop over the wells
  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion * const subRegion )
  {

    PerforationData * const perforationData = subRegion->GetPerforationData();
    perforationData->getReference<array2d<real64>>( viewKeyStruct::dPerforationRate_dPresString ).resizeDimension<1>(2);

  });

}
  
void SinglePhaseWell::UpdateState( WellElementSubRegion * const subRegion )
{
  SinglePhaseBase * const flowSolver = getParent()->GetGroup<SinglePhaseBase>( GetFlowSolverName() );

  GEOSX_ERROR_IF( flowSolver == nullptr,
                 "Flow solver " << GetFlowSolverName() << " not found in well solver " << getName() );

  flowSolver->UpdateFluidModel( subRegion );
}

void SinglePhaseWell::InitializeWells( DomainPartition * const domain )
{

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & resPressure = m_resPressure;
  
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & resDensity = m_resDensity;

  MeshLevel * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  // loop over the wells
  elemManager->forElementSubRegions<WellElementSubRegion>([&]( WellElementSubRegion * const subRegion )
  {
    WellControls const * const wellControls = GetWellControls( subRegion );
    PerforationData const * const perforationData = subRegion->GetPerforationData();

    // get the info stored on well elements
    arrayView1d<real64 const> const & wellElemGravCoef =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::gravityCoefString );

    // get well primary variables on well elements
    arrayView1d<real64> const & wellElemPressure =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    arrayView1d<real64> const & connRate  =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::connRateString );
    
    // get the element region, subregion, index
    arrayView1d<localIndex const> const & resElementRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );

    arrayView1d<localIndex const> const & resElementSubRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

    arrayView1d<localIndex const> const & resElementIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementIndexString );
    
    // define a reservoir pressure used for initialization
    real64 resPres = ( wellControls->GetType() == WellControls::Type::PRODUCER )
                   ? 1e20 : 0;

    // 1) Loop over all perforations to compute an average density
    real64 avgDensity = 0;
    for (localIndex iperf = 0; iperf < perforationData->size(); ++iperf)
    {
      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];
      
      avgDensity += resDensity[er][esr][m_resFluidIndex][ei][0];

      // save min pressure for producer
      if ( wellControls->GetType() == WellControls::Type::PRODUCER &&
           resPres > resPressure[er][esr][ei] )
      { 
        resPres = resPressure[er][esr][ei];
      }
      // save max pressure for injector
      else if ( wellControls->GetType() == WellControls::Type::INJECTOR &&
                resPres < resPressure[er][esr][ei] )
      { 
        resPres = resPressure[er][esr][ei];
      }
    }

    // communicate the pressures to the ranks without perforations
    // this will be used to initialize the pressure, starting by the owner rank
    if ( wellControls->GetType() == WellControls::Type::PRODUCER )
    { 
      resPres = MpiWrapper::Min( resPres );
    }
    else if ( wellControls->GetType() == WellControls::Type::INJECTOR )
    { 
      resPres = MpiWrapper::Max( resPres );
    }
    
    avgDensity = MpiWrapper::Sum( avgDensity );

    globalIndex const numPerforationsGlobal = perforationData->GetNumPerforationsGlobal();
    avgDensity /= numPerforationsGlobal;

    real64 pressureControl = 0.0;
    real64 gravCoefControl = 0.0;
    if (subRegion->IsLocallyOwned())
    {

      // get the reference data for this well
      localIndex const iwelemControl = wellControls->GetReferenceWellElementIndex();
      gravCoefControl = wellElemGravCoef[iwelemControl];

      // 2) Initialize the reference pressure
      real64 const & targetBHP = wellControls->GetTargetBHP();
      if (wellControls->GetControl() == WellControls::Control::BHP)
      {
        // if pressure constraint, set the ref pressure at the constraint
        pressureControl = targetBHP;
      }
      else // rate control
      {
        // if rate constraint, set the ref pressure slightly 
        // above/below the target pressure depending on well type
        pressureControl = 
          (wellControls->GetType() == WellControls::Type::PRODUCER)
          ? 0.5 * resPres // hard-coded values come from personal communication with Hui
          : 2.0 * resPres;
      }

      wellElemPressure[iwelemControl] = pressureControl;
    }

    // TODO optimize
    MpiWrapper::Broadcast( pressureControl, subRegion->GetTopRank() );
    MpiWrapper::Broadcast( gravCoefControl, subRegion->GetTopRank() );

    GEOSX_ERROR_IF( pressureControl <= 0, "Invalid well initialization: negative pressure was found" );

    // 3) Estimate the pressures in the well elements using this avgDensity
    forall_in_range( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex const iwelem )
    {
      wellElemPressure[iwelem] = pressureControl
        + avgDensity * ( wellElemGravCoef[iwelem] - gravCoefControl );
    });

    // 4) Recompute the pressure-dependent properties
    UpdateState( subRegion );

    // 5) Estimate the connection rates based on the min/max pressure
    forall_in_range( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex const iwelem )
    {
      real64 const & targetRate = wellControls->GetTargetRate();
      if (wellControls->GetControl() == WellControls::Control::BHP)
      {
        // if BHP constraint set rate below the absolute max rate 
        // with the appropriate sign (negative for prod, positive for inj)
        connRate[iwelem] = (wellControls->GetType() == WellControls::Type::PRODUCER)
          ? std::max( 0.1 * targetRate, - 1e3 ) // hard-coded values come from personal communication with Hui
          : std::min( 0.1 * targetRate,   1e3 );
      }
      else
      {
        connRate[iwelem] = targetRate;
      }
    });
  });
}


void SinglePhaseWell::SetupDofs( DomainPartition const * const domain,
                                 DofManager & dofManager ) const
{
  MeshLevel const * const meshLevel = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  array1d<string> regions;
  elemManager->forElementRegions<WellElementRegion>( [&]( WellElementRegion const * const region )
  {
    regions.push_back( region->getName() );
  } );

  dofManager.addField( WellElementDofName(),
                       DofManager::Location::Elem,
                       NumDofPerWellElement(),
                       regions );

  dofManager.addCoupling( WellElementDofName(),
                          WellElementDofName(),
                          DofManager::Connectivity::Node );
}

void SinglePhaseWell::AssembleFluxTerms( real64 const GEOSX_UNUSED_ARG( time_n ),
                                         real64 const dt,
                                         DomainPartition const * const domain,
                                         DofManager const * const dofManager,
                                         ParallelMatrix * const matrix,
                                         ParallelVector * const rhs )
{

  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  string const wellDofKey = dofManager->getKey( WellElementDofName() );

  // loop over the wells
  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion const * const subRegion )
  {

    // get a reference to the degree-of-freedom numbers
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      subRegion->getReference<array1d<globalIndex>>( wellDofKey );

    arrayView1d<integer const> const & wellElemGhostRank =
      subRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

    arrayView1d<localIndex const> const & nextWellElemIndex =
      subRegion->getReference<array1d<localIndex>>( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );

    // get a reference to the primary variables on well elements
    arrayView1d<real64 const> const & connRate =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::connRateString );

    arrayView1d<real64 const> const & dConnRate =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::deltaConnRateString );

    // loop over the well elements to compute the fluxes between elements
    forall_in_range( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex const iwelem )
    {

      if ( wellElemGhostRank[iwelem] >= 0 )
      {
        return;
      }

      // 1) Compute the flux and its derivatives

      /*  currentConnRate < 0 flow from iwelem to iwelemNext
       *  currentConnRate > 0 flow from iwelemNext to iwelem
       *  With this convention, currentConnRate < 0 at the last connection for a producer
       *                        currentConnRate > 0 at the last connection for a injector
       */

      // get next well element index
      localIndex const iwelemNext = nextWellElemIndex[iwelem];

      // there is nothing to upwind for single-phase flow
      real64 const currentConnRate = connRate[iwelem] + dConnRate[iwelem];
      real64 const flux        = dt * currentConnRate;
      real64 const dFlux_dRate = dt;

      // 2) Assemble the flux into residual and Jacobian
      if ( iwelemNext < 0 )
      {
        // flux terms
        real64 const oneSidedLocalFlux = - flux;
        real64 const oneSidedLocalFluxJacobian_dRate = - dFlux_dRate;

        // jacobian indices
        globalIndex const offset = wellElemDofNumber[iwelem];
        globalIndex const oneSidedEqnRowIndex       = offset + RowOffset::MASSBAL;
        globalIndex const oneSidedDofColIndex_dRate = offset + ColOffset::DRATE;

        rhs->add( &oneSidedEqnRowIndex,
                  &oneSidedLocalFlux,
                  1 );

        matrix->add( &oneSidedEqnRowIndex,
                     &oneSidedDofColIndex_dRate,
                     &oneSidedLocalFluxJacobian_dRate,
                     1, 1 );

      }
      else
      {
        // local working variables and arrays
        stackArray1d<globalIndex, 2> eqnRowIndices( 2 );

        stackArray1d<real64, 2> localFlux( 2 );
        stackArray1d<real64, 2> localFluxJacobian_dRate( 2 );

        // flux terms
        localFlux[ElemTag::NEXT]    =   flux;
        localFlux[ElemTag::CURRENT] = - flux;

        localFluxJacobian_dRate[ElemTag::NEXT]    =   dFlux_dRate;
        localFluxJacobian_dRate[ElemTag::CURRENT] = - dFlux_dRate;

        // indices
        globalIndex const offsetCurrent = wellElemDofNumber[iwelem];
        globalIndex const offsetNext    = wellElemDofNumber[iwelemNext];
        eqnRowIndices[ElemTag::CURRENT] = offsetCurrent + RowOffset::MASSBAL;
        eqnRowIndices[ElemTag::NEXT]    = offsetNext    + RowOffset::MASSBAL;
        globalIndex const dofColIndex_dRate = offsetCurrent + ColOffset::DRATE;

        rhs->add( eqnRowIndices.data(),
                  localFlux.data(),
                  2 );

        matrix->add( eqnRowIndices.data(),
                     &dofColIndex_dRate,
                     localFluxJacobian_dRate.data(),
                     2, 1 );
      }
    });
  });
}


void SinglePhaseWell::AssemblePerforationTerms( real64 const GEOSX_UNUSED_ARG( time_n ),
                                                real64 const dt,
                                                DomainPartition const * const domain,
                                                DofManager const * const dofManager,
                                                ParallelMatrix * const matrix,
                                                ParallelVector * const rhs )
{
  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  string const wellDofKey = dofManager->getKey( WellElementDofName() );
  string const resDofKey  = dofManager->getKey( ResElementDofName() );

  ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex> > resDofNumberAccessor =
    elemManager->ConstructViewAccessor< array1d<globalIndex>, arrayView1d<globalIndex> >( resDofKey );

  ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex const> >::ViewTypeConst resDofNumber =
    resDofNumberAccessor.toViewConst();

  // loop over the wells
  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion const * const subRegion )
  {
    PerforationData const * const perforationData = subRegion->GetPerforationData();

    // compute the local rates for this well
    ComputeAllPerforationRates( subRegion );

    // get the degrees of freedom
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      subRegion->getReference<array1d<globalIndex>>( wellDofKey );

    // get well variables on perforations
    arrayView1d<real64 const> const & perfRate =
      perforationData->getReference<array1d<real64>>( viewKeyStruct::perforationRateString );

    arrayView2d<real64 const> const & dPerfRate_dPres =
      perforationData->getReference<array2d<real64>>( viewKeyStruct::dPerforationRate_dPresString );

    arrayView1d<localIndex const> const & perfWellElemIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::wellElementIndexString );

    // get the element region, subregion, index
    arrayView1d<localIndex const> const & resElementRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );

    arrayView1d<localIndex const> const & resElementSubRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

    arrayView1d<localIndex const> const & resElementIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementIndexString );

    // local working variables and arrays
    stackArray1d<globalIndex, 2> eqnRowIndices( 2 );
    stackArray1d<globalIndex, 2> dofColIndices( 2 );

    stackArray1d<real64, 2> localPerf( 2 );
    stackArray2d<real64, 4> localPerfJacobian(2, 2);

    // loop over the perforations and add the rates to the residual and jacobian
    for (localIndex iperf = 0; iperf < perforationData->size(); ++iperf)
    {
      eqnRowIndices = -1;
      dofColIndices = -1;

      localPerf = 0;
      localPerfJacobian = 0;

      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];

      // get the well element index for this perforation
      localIndex const iwelem      = perfWellElemIndex[iperf];
      globalIndex const elemOffset = wellElemDofNumber[iwelem];

      // row index on reservoir side
      eqnRowIndices[SubRegionTag::RES] = resDofNumber[er][esr][ei];

      // column index on reservoir side
      dofColIndices[SubRegionTag::RES] = resDofNumber[er][esr][ei];

      // row index on well side
      eqnRowIndices[SubRegionTag::WELL] = elemOffset + RowOffset::MASSBAL;

      // column index on well side
      dofColIndices[SubRegionTag::WELL] = elemOffset + ColOffset::DPRES;

      // populate local flux vector and derivatives
      localPerf[SubRegionTag::RES]  =  dt * perfRate[iperf];
      localPerf[SubRegionTag::WELL] = -localPerf[SubRegionTag::RES];

      for (localIndex ke = 0; ke < 2; ++ke)
      {
        localPerfJacobian[SubRegionTag::RES][ke]  = dt * dPerfRate_dPres[iperf][ke];
        localPerfJacobian[SubRegionTag::WELL][ke] = - localPerfJacobian[SubRegionTag::RES][ke];
      }

      rhs->add( eqnRowIndices.data(),
                localPerf.data(),
                2 );

      matrix->add( eqnRowIndices.data(),
                   dofColIndices.data(),
                   localPerfJacobian.data(),
                   2, 2 );
    }

  });
}


void SinglePhaseWell::FormPressureRelations( DomainPartition const * const domain,
                                             DofManager const * const dofManager,
                                             ParallelMatrix * const matrix,
                                             ParallelVector * const rhs )
{

  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  string const wellDofKey = dofManager->getKey( WellElementDofName() );

  elemManager->forElementSubRegions<WellElementSubRegion>([&]( WellElementSubRegion const * const subRegion )
  {

    WellControls const * const wellControls = GetWellControls( subRegion );

    // get the degrees of freedom numbers, depth, next well elem index
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      subRegion->getReference<array1d<globalIndex>>( wellDofKey );

    arrayView1d<integer const> const & wellElemGhostRank =
      subRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

    arrayView1d<real64 const> const & wellElemGravCoef =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::gravityCoefString );

    arrayView1d<localIndex const> const & nextWellElemIndex =
      subRegion->getReference<array1d<localIndex>>( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );

    // get primary variables on well elements
    arrayView1d<real64 const> const & wellElemPressure =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    arrayView1d<real64 const> const & dWellElemPressure =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    // get well constitutive data
    SingleFluidBase const * const fluid = GetConstitutiveModel<SingleFluidBase>( subRegion, m_fluidName );

    arrayView2d<real64 const> const & wellElemDensity =
      fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::densityString );

    arrayView2d<real64 const> const & dWellElemDensity_dPres =
      fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::dDens_dPresString );

    // loop over the well elements to compute the pressure relations between well elements
    forall_in_range( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex const iwelem )
    {

      if (wellElemGhostRank[iwelem] >= 0)
      {
        return;
      }

      localIndex const iwelemNext = nextWellElemIndex[iwelem];

      if ( iwelemNext >= 0 )  // if iwelemNext < 0, form control equation, not momentum
      {

        // local working variables and arrays
        stackArray1d<globalIndex, 2> dofColIndices( 2 );
        stackArray1d<real64, 2> localPresRelJacobian( 2 );

        dofColIndices        = -1;
        localPresRelJacobian = 0;

        // compute avg density
        real64 const avgDensity = 0.5 * ( wellElemDensity[iwelem][0] + wellElemDensity[iwelemNext][0] );
        real64 const dAvgDensity_dPresNext    = 0.5 * dWellElemDensity_dPres[iwelemNext][0];
        real64 const dAvgDensity_dPresCurrent = 0.5 * dWellElemDensity_dPres[iwelem][0];

        // compute depth diff times acceleration
        real64 const gravD = wellElemGravCoef[iwelemNext] - wellElemGravCoef[iwelem];

        // compute the current pressure in the two well elements
        real64 const pressureCurrent = wellElemPressure[iwelem]     + dWellElemPressure[iwelem];
        real64 const pressureNext    = wellElemPressure[iwelemNext] + dWellElemPressure[iwelemNext];

        // compute a coefficient to normalize the momentum equation
        real64 const & targetBHP  = wellControls->GetTargetBHP();
        real64 const normalizer = targetBHP > std::numeric_limits<real64>::epsilon()
                                ? 1.0 / targetBHP
                                : 1.0;

        // compute momentum flux and derivatives
        real64 const localPresRel = ( pressureNext - pressureCurrent - avgDensity * gravD ) * normalizer;
        localPresRelJacobian[ElemTag::NEXT]    = ( 1 - dAvgDensity_dPresNext * gravD ) * normalizer;
        localPresRelJacobian[ElemTag::CURRENT] = (-1 - dAvgDensity_dPresCurrent * gravD ) * normalizer;

        // TODO: add friction and acceleration terms

        // jacobian indices
        globalIndex const offsetNext     = wellElemDofNumber[iwelemNext];
        globalIndex const offsetCurrent  = wellElemDofNumber[iwelem];
        globalIndex const eqnRowIndex    = offsetCurrent + RowOffset::CONTROL;
        dofColIndices[ElemTag::NEXT]     = offsetNext    + ColOffset::DPRES;
        dofColIndices[ElemTag::CURRENT]  = offsetCurrent + ColOffset::DPRES;

        rhs->add( &eqnRowIndex,
                  &localPresRel,
                  1 );

        matrix->add( &eqnRowIndex,
                     dofColIndices.data(),
                     localPresRelJacobian.data(),
                     1, 2 );
      }
    });
  });
}

void SinglePhaseWell::AssembleVolumeBalanceTerms( real64 const GEOSX_UNUSED_ARG( time_n ),
                                                  real64 const GEOSX_UNUSED_ARG( dt ),
                                                  DomainPartition const * const GEOSX_UNUSED_ARG( domain ),
                                                  DofManager const * const GEOSX_UNUSED_ARG( dofManager ),
                                                  ParallelMatrix * const GEOSX_UNUSED_ARG( matrix ),
                                                  ParallelVector * const GEOSX_UNUSED_ARG( rhs ) )
{
  // not implemented for single phase flow
}


void SinglePhaseWell::CheckWellControlSwitch( DomainPartition * const domain )
{
  MeshLevel * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  // loop over the wells
  elemManager->forElementSubRegions<WellElementSubRegion>([&]( WellElementSubRegion const * const subRegion )
  {

    if (! subRegion->IsLocallyOwned())
    {
      return;
    }

    WellControls * const wellControls = GetWellControls( subRegion );

    // get the primary variables
    arrayView1d<real64 const> const & wellElemPressure =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    arrayView1d<real64 const> const & dWellElemPressure =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    arrayView1d<real64 const> const & connRate  =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::connRateString );

    arrayView1d<real64 const> const & dConnRate =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::deltaConnRateString );

    // if isViable is true at the end of the following checks, no need to switch
    bool controlIsViable = false;

    // get well control and type
    WellControls::Control const currentControl = wellControls->GetControl();
    WellControls::Type const type = wellControls->GetType();

    // again we assume here that the first well element is on this MPI rank
    localIndex const iwelemControl = wellControls->GetReferenceWellElementIndex();

    real64 const refRate = connRate[iwelemControl] + dConnRate[iwelemControl];
    real64 const refPressure = wellElemPressure[iwelemControl] + dWellElemPressure[iwelemControl];

    // BHP control
    if ( currentControl == WellControls::Control::BHP )
    {
      // the control is viable if the reference rate is below the max rate
      real64 const & maxRate = wellControls->GetTargetRate();
      controlIsViable = ( fabs(refRate) <= fabs(maxRate) );
    }
    else // rate control
    {
      // the control is viable if the reference pressure is below/above the max/min pressure
      if ( type == WellControls::Type::PRODUCER )
      {
        // targetBHP specifies a min pressure here
        real64 const & minPressure = wellControls->GetTargetBHP();
        controlIsViable = ( refPressure >= minPressure );
      }
      else
      {
        // targetBHP specifies a max pressure here
        real64 const & maxPressure = wellControls->GetTargetBHP();
        controlIsViable = ( refPressure <= maxPressure );
      }
    }

    if (!controlIsViable)
    {
      if ( currentControl == WellControls::Control::BHP )
      {
        wellControls->SetControl( WellControls::Control::LIQUIDRATE,
                                  wellControls->GetTargetRate() );
        GEOSX_LOG_LEVEL_RANK_0( 1, "Control switch for well " << subRegion->getName()
                          << " from BHP constraint to rate constraint" );
      }
      else // rate control
      {
        wellControls->SetControl( WellControls::Control::BHP,
                                  wellControls->GetTargetBHP() );

        // Debug information for logLevel >= 1 
        GEOSX_LOG_LEVEL_RANK_0( 1, "Control switch for well " << subRegion->getName()
                              << " from rate constraint to BHP constraint" );
      }
    }
  });
}


real64
SinglePhaseWell::CalculateResidualNorm( DomainPartition const * const domain,
                                        DofManager const & dofManager,
                                        ParallelVector const & rhs )
{
  // get a view into local residual vector
  real64 const * localResidual = rhs.extractLocalVector();

  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  string const wellDofKey = dofManager.getKey( WellElementDofName() );

  real64 residualNorm = 0;
  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion const * const subRegion )
  {

    // get the degree of freedom numbers
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      subRegion->getReference<array1d<globalIndex>>( wellDofKey );

    arrayView1d<integer const> const & wellElemGhostRank =
      subRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

    arrayView1d<real64 const> const & wellElemVolume =
      subRegion->getReference<array1d<real64>>( ElementSubRegionBase::viewKeyStruct::elementVolumeString );

    SingleFluidBase const * const fluid = GetConstitutiveModel<SingleFluidBase>( subRegion, m_fluidName );

    arrayView2d<real64 const> const & wellElemDensity =
      fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::densityString );

    for (localIndex iwelem = 0; iwelem < subRegion->size(); ++iwelem)
    {
      if (wellElemGhostRank[iwelem] < 0)
      {
        for (localIndex idof = 0; idof < NumDofPerWellElement(); ++idof)
        {
          real64 const normalizer = (idof == RowOffset::MASSBAL)
                                  ? wellElemDensity[iwelem][0] * wellElemVolume[iwelem]
                                  : 1;
          localIndex const lid    = rhs.getLocalRowID( wellElemDofNumber[iwelem] + idof );
          real64 const val = localResidual[lid] / normalizer;
          residualNorm += val * val;
        }
      }
    }
  });

  // compute global residual norm
  real64 globalResidualNorm;
  MpiWrapper::allReduce(&residualNorm, &globalResidualNorm, 1, MPI_SUM, MPI_COMM_GEOSX);

  return sqrt(globalResidualNorm);
}

bool
SinglePhaseWell::CheckSystemSolution( DomainPartition const * const domain,
                                      DofManager const & dofManager,
                                      ParallelVector const & solution,
                                      real64 const scalingFactor )
{
  // get the update
  real64 const * localSolution = solution.extractLocalVector();

  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  int isInvalidLocal = 0;

  string const wellDofKey = dofManager.getKey( WellElementDofName() );

  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion const * const subRegion )
  {

    // get the degree of freedom numbers on well elements
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      subRegion->getReference<array1d<globalIndex>>( wellDofKey );

    arrayView1d<integer const> const & wellElemGhostRank =
      subRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

    // get a reference to the primary variables on well elements
    arrayView1d<real64 const> const & wellElemPressure =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    arrayView1d<real64 const> const & dWellElemPressure =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    for (localIndex iwelem = 0; iwelem < subRegion->size(); ++iwelem)
    {

      if (wellElemGhostRank[iwelem] < 0)
      {
        // extract pressure solution and apply to dP
        localIndex const lid =  solution.getLocalRowID( wellElemDofNumber[iwelem] + ColOffset::DPRES );
        real64 const newPres = wellElemPressure[iwelem] + dWellElemPressure[iwelem]
                             + scalingFactor * localSolution[lid];

        if (newPres < 0.0)
        {
          isInvalidLocal = 1;
        }
      }
    }
  });

  int isInvalidGlobal;
  MpiWrapper::allReduce(&isInvalidLocal, &isInvalidGlobal, 1, MPI_SUM, MPI_COMM_GEOSX);
 
  bool isValid = (isInvalidGlobal == 0);
  return isValid;
}

void
SinglePhaseWell::ApplySystemSolution( DofManager const & dofManager,
                                      ParallelVector const & solution,
                                      real64 const scalingFactor,
                                      DomainPartition * const domain )
{
  dofManager.addVectorToField( solution,
                               WellElementDofName(),
                               viewKeyStruct::deltaPressureString,
                               scalingFactor,
                               0, 1 );

  dofManager.addVectorToField( solution,
                               WellElementDofName(),
                               viewKeyStruct::deltaConnRateString,
                               scalingFactor,
                               1, m_numDofPerWellElement );

  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::deltaPressureString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaConnRateString );
  CommunicationTools::SynchronizeFields(fieldNames,
                                        domain->getMeshBody(0)->getMeshLevel(0),
                                        domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  // update properties
  UpdateStateAll( domain );
}

void SinglePhaseWell::ResetStateToBeginningOfStep( DomainPartition * const domain )
{

  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion const * const subRegion )
  {

    // get a reference to the primary variables on well elements
    arrayView1d<real64> const & dWellElemPressure =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    arrayView1d<real64> const & dConnRate =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::deltaConnRateString );

    forall_in_range( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex const iwelem )
    {
      dWellElemPressure[iwelem] = 0;
      dConnRate[iwelem] = 0;
    });
  });

  // call constitutive models
  UpdateStateAll( domain );
}


void SinglePhaseWell::ResetViews(DomainPartition * const domain)
{
  WellSolverBase::ResetViews(domain);

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  m_resPressure =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( SinglePhaseBase::viewKeyStruct::pressureString );

  m_deltaResPressure =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( SinglePhaseBase::viewKeyStruct::deltaPressureString );

  m_resDensity =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( SingleFluidBase::viewKeyStruct::densityString,
                                                                                          constitutiveManager );
  m_dResDens_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( SingleFluidBase::viewKeyStruct::dDens_dPresString,
                                                                                          constitutiveManager );
  m_resViscosity =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( SingleFluidBase::viewKeyStruct::viscosityString,
                                                                                          constitutiveManager );
  m_dResVisc_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( SingleFluidBase::viewKeyStruct::dVisc_dPresString,
                                                                                          constitutiveManager );
}

void SinglePhaseWell::ComputeAllPerforationRates( WellElementSubRegion const * const subRegion )
{

  // get the reservoir data
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & resPressure         = m_resPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & dResPressure        = m_deltaResPressure;

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & resDensity          = m_resDensity;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dResDensity_dPres   = m_dResDens_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & resViscosity        = m_resViscosity;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dResViscosity_dPres = m_dResVisc_dPres;

  // get the well data
  PerforationData const * const perforationData = subRegion->GetPerforationData();

  // get the degrees of freedom and depth
  arrayView1d<real64 const> const & wellElemGravCoef =
    subRegion->getReference<array1d<real64>>( viewKeyStruct::gravityCoefString );

  // get well primary variables on well elements
  arrayView1d<real64 const> const & wellElemPressure =
    subRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

  arrayView1d<real64 const> const & dWellElemPressure =
    subRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

  // get well constitutive data
  SingleFluidBase const * const fluid = GetConstitutiveModel<SingleFluidBase>( subRegion, m_fluidName );

  arrayView2d<real64 const> const & wellElemDensity =
    fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::densityString );

  arrayView2d<real64 const> const & dWellElemDensity_dPres =
    fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::dDens_dPresString );

  arrayView2d<real64 const> const & wellElemViscosity =
    fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::viscosityString );

  arrayView2d<real64 const> const & dWellElemViscosity_dPres =
    fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::dVisc_dPresString );

  // get well variables on perforations
  arrayView1d<real64 const> const & perfGravCoef =
    perforationData->getReference<array1d<real64>>( viewKeyStruct::gravityCoefString );

  arrayView1d<localIndex const> const & perfWellElemIndex =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::wellElementIndexString );

  arrayView1d<real64 const> const & perfTransmissibility =
    perforationData->getReference<array1d<real64>>( PerforationData::viewKeyStruct::transmissibilityString );

  arrayView1d<real64> const & perfRate =
    perforationData->getReference<array1d<real64>>( viewKeyStruct::perforationRateString );

  arrayView2d<real64> const & dPerfRate_dPres =
    perforationData->getReference<array2d<real64>>( viewKeyStruct::dPerforationRate_dPresString );

  // get the element region, subregion, index
  arrayView1d<localIndex const> const & resElementRegion =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );

  arrayView1d<localIndex const> const & resElementSubRegion =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

  arrayView1d<localIndex const> const & resElementIndex =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementIndexString );

  // local working variables and arrays
  stackArray1d<globalIndex, 2> eqnRowIndices( 2 );
  stackArray1d<globalIndex, 2> dofColIndices( 2 );

  stackArray1d<real64, 2> pressure( 2 );
  stackArray1d<real64, 2> dPressure_dP( 2 );

  stackArray1d<localIndex, 2> multiplier( 2 );

  // loop over the perforations to compute the perforation rates
  for (localIndex iperf = 0; iperf < perforationData->size(); ++iperf)
  {
    eqnRowIndices = -1;
    dofColIndices = -1;

    pressure = 0;
    dPressure_dP = 0;

    multiplier = 0;

    // 1) Reservoir side

    // get the reservoir (sub)region and element indices
    localIndex const er  = resElementRegion[iperf];
    localIndex const esr = resElementSubRegion[iperf];
    localIndex const ei  = resElementIndex[iperf];

    // get reservoir variables
    pressure[SubRegionTag::RES] = resPressure[er][esr][ei] + dResPressure[er][esr][ei];
    dPressure_dP[SubRegionTag::RES] = 1;

    // TODO: add a buoyancy term for the reservoir side here

    // multiplier for reservoir side in the flux
    multiplier[SubRegionTag::RES] = 1;

    // 2) Well side

    // get the local index of the well element
    localIndex const iwelem = perfWellElemIndex[iperf];

    // get well variables
    pressure[SubRegionTag::WELL] = wellElemPressure[iwelem] + dWellElemPressure[iwelem];
    dPressure_dP[SubRegionTag::WELL] = 1.0;

    real64 const gravD = ( perfGravCoef[iperf] - wellElemGravCoef[iwelem] );
    pressure[SubRegionTag::WELL]     += wellElemDensity[iwelem][0] * gravD;
    dPressure_dP[SubRegionTag::WELL] += dWellElemDensity_dPres[iwelem][0] * gravD;

    // multiplier for well side in the flux
    multiplier[SubRegionTag::WELL] = -1;

    // get transmissibility at the interface
    real64 const trans = perfTransmissibility[iperf];

    // compute potential difference
    real64 potDif = 0.0;
    for (localIndex i = 0; i < 2; ++i)
    {
      potDif += multiplier[i] * trans * pressure[i];
      dPerfRate_dPres[iperf][i] = multiplier[i] * trans * dPressure_dP[i];
    }

    // choose upstream cell based on potential difference
    localIndex const k_up = (potDif >= 0) ? SubRegionTag::RES : SubRegionTag::WELL;

    // compute upstream density, viscosity, and mobility
    real64 densityUp       = 0.0;
    real64 dDensityUp_dP   = 0.0;
    real64 viscosityUp     = 0.0;
    real64 dViscosityUp_dP = 0.0;

    // upwinding the variables
    if (k_up == SubRegionTag::RES) // use reservoir vars
    {
      densityUp     = resDensity[er][esr][m_resFluidIndex][ei][0];
      dDensityUp_dP = dResDensity_dPres[er][esr][m_resFluidIndex][ei][0];

      viscosityUp     = resViscosity[er][esr][m_resFluidIndex][ei][0];
      dViscosityUp_dP = dResViscosity_dPres[er][esr][m_resFluidIndex][ei][0];
    }
    else // use well vars
    {
      densityUp = wellElemDensity[iwelem][0];
      dDensityUp_dP = dWellElemDensity_dPres[iwelem][0];

      viscosityUp = wellElemViscosity[iwelem][0];
      dViscosityUp_dP = dWellElemViscosity_dPres[iwelem][0];
    }

    // compute mobility
    real64 const mobilityUp     = densityUp / viscosityUp;
    real64 const dMobilityUp_dP = dDensityUp_dP / viscosityUp
                                - mobilityUp / viscosityUp * dViscosityUp_dP;

    perfRate[iperf] = mobilityUp * potDif;
    for (localIndex ke = 0; ke < 2; ++ke)
    {
      dPerfRate_dPres[iperf][ke] *= mobilityUp;
    }
    dPerfRate_dPres[iperf][k_up] += dMobilityUp_dP * potDif;
  }
}


void SinglePhaseWell::FormControlEquation( DomainPartition const * const domain,
                                           DofManager const * const dofManager,
                                           ParallelMatrix * const matrix,
                                           ParallelVector * const rhs )
{
  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  string const wellDofKey = dofManager->getKey( WellElementDofName() );

  elemManager->forElementSubRegions<WellElementSubRegion>([&]( WellElementSubRegion const * const subRegion )
  {

    if (! subRegion->IsLocallyOwned())
    {
      return;
    }

    WellControls const * const wellControls = GetWellControls( subRegion );

    // get the degrees of freedom
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      subRegion->getReference<array1d<globalIndex>>( wellDofKey );

    // get the index of the well element at which the control is enforced
    localIndex const iwelemControl = wellControls->GetReferenceWellElementIndex();

    // get well control
    WellControls::Control const control = wellControls->GetControl();

    // BHP control
    if (control == WellControls::Control::BHP)
    {
      // get primary variables on well elements
      arrayView1d<real64 const> const & wellElemPressure =
        subRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

      arrayView1d<real64 const> const & dWellElemPressure =
        subRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

      // get pressures and compute normalizer
      real64 const currentBHP = wellElemPressure[iwelemControl] + dWellElemPressure[iwelemControl];
      real64 const & targetBHP = wellControls->GetTargetBHP();
      real64 const normalizer = targetBHP > std::numeric_limits<real64>::epsilon()
                              ? 1.0 / targetBHP
                              : 1.0;

      // control equation is a normalized difference between current pressure and target pressure
      real64 const controlEqn = ( currentBHP - targetBHP ) * normalizer;
      real64 const dControlEqn_dPres = normalizer;

      globalIndex const elemOffset  = wellElemDofNumber[iwelemControl];
      globalIndex const eqnRowIndex = elemOffset + RowOffset::CONTROL;
      globalIndex const dofColIndex = elemOffset + ColOffset::DPRES;

      rhs->add( &eqnRowIndex,
                &controlEqn,
                1 );

      matrix->add( &eqnRowIndex,
                   &dofColIndex,
                   &dControlEqn_dPres,
                   1, 1 );

    }
    // rate control
    else
    {

      // get a reference to the primary variables on well element
      arrayView1d<real64 const> const & connRate  =
        subRegion->getReference<array1d<real64>>( viewKeyStruct::connRateString );

      arrayView1d<real64 const> const & dConnRate =
        subRegion->getReference<array1d<real64>>( viewKeyStruct::deltaConnRateString );

      // get rates and compute normalizer
      real64 const currentConnRate = connRate[iwelemControl] + dConnRate[iwelemControl];
      real64 const & targetConnRate = wellControls->GetTargetRate();
      real64 const normalizer = fabs(targetConnRate) > std::numeric_limits<real64>::min()
                              ? 1.0 / ( 1e-2 * fabs(targetConnRate) ) // hard-coded value comes from AD-GPRS
                              : 1.0;

      // for a producer, the actual (target) rate is negative

      // control equation is a normalized difference between current rate and target rate
      real64 const controlEqn = ( currentConnRate - targetConnRate ) * normalizer;
      real64 const dControlEqn_dRate = normalizer;

      globalIndex const elemOffset  = wellElemDofNumber[iwelemControl];
      globalIndex const eqnRowIndex = elemOffset + RowOffset::CONTROL;
      globalIndex const dofColIndex = elemOffset + ColOffset::DRATE;

      rhs->add( &eqnRowIndex,
                &controlEqn,
                1 );

      matrix->add( &eqnRowIndex,
                   &dofColIndex,
                   &dControlEqn_dRate,
                   1, 1 );
    }
  });
}


void SinglePhaseWell::ImplicitStepComplete( real64 const & GEOSX_UNUSED_ARG( time ),
                                            real64 const & GEOSX_UNUSED_ARG( dt ),
                                            DomainPartition * const domain )
{
  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion const * const subRegion )
  {

    // get a reference to the primary variables on well elements
    arrayView1d<real64> const & wellElemPressure  =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    arrayView1d<real64 const> const & dWellElemPressure =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    arrayView1d<real64> const & connRate  =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::connRateString );

    arrayView1d<real64 const> const & dConnRate =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::deltaConnRateString );

    for (localIndex iwelem = 0; iwelem < subRegion->size(); ++iwelem)
    {
      wellElemPressure[iwelem] += dWellElemPressure[iwelem];
      connRate[iwelem]         += dConnRate[iwelem];
    }

    // TODO: improve well data output
    /*
    int mpiSize = CommunicationTools::Comm_size(MPI_COMM_GEOSX) ;
    if (mpiSize == 1)
    {
      RecordWellData( subRegion );
    }
    */
  });
}

void SinglePhaseWell::RecordWellData( WellElementSubRegion const * const subRegion )
{
  // Note: this function is for debug and will go away

  WellControls const * const wellControls = GetWellControls( subRegion );
  PerforationData const * const perforationData = subRegion->GetPerforationData();

  arrayView1d<real64 const> const & wellElemPressure =
    subRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

  arrayView1d<real64 const> const & connRate =
    subRegion->getReference<array1d<real64>>( viewKeyStruct::connRateString );


  arrayView1d<real64 const> const & perfRate =
    perforationData->getReference<array1d<real64>>( viewKeyStruct::perforationRateString );

  // here, we will save the well info
  // for now, brute force: output to terminal

  std::cout << "Well : " << subRegion->getName() << std::endl;
  if (wellControls->GetType() == WellControls::Type::PRODUCER)
    std::cout << "Type : PRODUCER" << std::endl;
  else
    std::cout << "Type : INJECTOR" << std::endl;
  if (wellControls->GetControl() == WellControls::Control::BHP)
    std::cout << "Control : BHP" << std::endl;
  else
    std::cout << "Control : RATE" << std::endl;

  std::cout << "Below, positive perforation rate means flow from reservoir to well" << std::endl;
  std::cout << "Negative perforation rate means flow from well to reservoir" << std::endl;

  // output perforation rates
  for (localIndex iperf = 0; iperf < perforationData->size(); ++iperf)
  {
    std::cout << "Mass rate at perforation #" << iperf << ": " << perfRate[iperf] << std::endl;
  }

  // output the reference pressure
  if (subRegion->IsLocallyOwned())
  {
    localIndex const iwelemControl = wellControls->GetReferenceWellElementIndex();
    real64 const pressure = wellElemPressure[iwelemControl];
    real64 const targetPressure = wellControls->GetTargetBHP();

    if (wellControls->GetControl() == WellControls::Control::BHP)
    {
      std::cout << "Current reference pressure = " << pressure
                << ", targetPressure = "           << targetPressure
                << std::endl;
    }
    else
    {
      if (wellControls->GetType() == WellControls::Type::PRODUCER)
      {
        std::cout << "Current reference pressure = " << pressure
                  << ", min pressure = " << targetPressure
                  << std::endl;
      }
      else
      {
        std::cout << "Current reference pressure = " << pressure
                  << ", max pressure = " << targetPressure
                  << std::endl;
      }
    }
  }

  std::cout << "Below, negative connection rate means production" << std::endl;
  std::cout << "Positive connection rate means injection" << std::endl;

  for (localIndex iwelem = 0; iwelem < subRegion->size(); ++iwelem)
  {
    if (iwelem > 0 || wellControls->GetControl() == WellControls::Control::BHP)
      std::cout << "Mass rate at connection #"
                << iwelem
                << ": "
                << connRate[iwelem]
                << std::endl;
    else
      std::cout << "Mass rate at connection #" << iwelem << ": " << connRate[iwelem]
                << ", target rate : " << wellControls->GetTargetRate() << std::endl;
  }
}

REGISTER_CATALOG_ENTRY(SolverBase, SinglePhaseWell, string const &, Group * const)
}// namespace geosx
