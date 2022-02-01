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
 * @file SinglePhaseWellKernels.cpp
 */

#include "SinglePhaseWellKernels.hpp"

// TODO: move keys to WellControls
#include "SinglePhaseWell.hpp"

namespace geosx
{

namespace SinglePhaseWellKernels
{

/******************************** ControlEquationHelper ********************************/

GEOSX_HOST_DEVICE
void
ControlEquationHelper::
  switchControl( WellControls::Type const & wellType,
                 WellControls::Control const & currentControl,
                 real64 const & targetBHP,
                 real64 const & targetRate,
                 real64 const & currentBHP,
                 real64 const & currentVolRate,
                 WellControls::Control & newControl )
{
  // if isViable is true at the end of the following checks, no need to switch
  bool controlIsViable = false;

  // The limiting flow rates are treated as upper limits, while the pressure limits
  // are treated as lower limits in production wells and upper limits in injectors.
  // The well changes its mode of control whenever the existing control mode would
  // violate one of these limits.
  // BHP control
  if( currentControl == WellControls::Control::BHP )
  {
    // the control is viable if the reference rate is below the max rate
    controlIsViable = ( LvArray::math::abs( currentVolRate ) <= LvArray::math::abs( targetRate ) + EPS );
  }
  else // rate control
  {
    // the control is viable if the reference pressure is below/above the max/min pressure
    if( wellType == WellControls::Type::PRODUCER )
    {
      // targetBHP specifies a min pressure here
      controlIsViable = ( currentBHP >= targetBHP - EPS );
    }
    else
    {
      // targetBHP specifies a max pressure here
      controlIsViable = ( currentBHP <= targetBHP + EPS );
    }
  }

  if( controlIsViable )
  {
    newControl = currentControl;
  }
  else
  {
    // Note: if BHP control is not viable, we switch to TOTALVOLRATE
    //       if TOTALVOLRATE are not viable, we switch to BHP
    newControl = ( currentControl == WellControls::Control::BHP )
               ? WellControls::Control::TOTALVOLRATE
               : WellControls::Control::BHP;
  }
}

GEOSX_HOST_DEVICE
void
ControlEquationHelper::
  compute( globalIndex const rankOffset,
           WellControls::Control const currentControl,
           real64 const & targetBHP,
           real64 const & targetRate,
           real64 const & currentBHP,
           real64 const & dCurrentBHP_dPres,
           real64 const & currentVolRate,
           real64 const & dCurrentVolRate_dPres,
           real64 const & dCurrentVolRate_dRate,
           globalIndex const dofNumber,
           CRSMatrixView< real64, globalIndex const > const & localMatrix,
           arrayView1d< real64 > const & localRhs )
{
  localIndex const eqnRowIndex = dofNumber + ROFFSET::CONTROL - rankOffset;
  globalIndex const presDofColIndex = dofNumber + COFFSET::DPRES;
  globalIndex const rateDofColIndex = dofNumber + COFFSET::DRATE;

  real64 controlEqn = 0;
  real64 dControlEqn_dRate = 0;
  real64 dControlEqn_dPres = 0;

  // Note: We assume in the computation of currentBHP that the reference elevation
  //       is in the top well element. This is enforced by a check in the solver.
  //       If we wanted to allow the reference elevation to be outside the top
  //       well element, it would make more sense to check the BHP constraint in
  //       the well element that contains the reference elevation.

  // BHP control
  if( currentControl == WellControls::Control::BHP )
  {
    // control equation is a difference between current BHP and target BHP
    controlEqn = currentBHP - targetBHP;
    dControlEqn_dPres = dCurrentBHP_dPres;
  }
  // Total volumetric rate control
  else if( currentControl == WellControls::Control::TOTALVOLRATE )
  {
    // control equation is the difference between volumetric current rate and target rate
    controlEqn = currentVolRate - targetRate;
    dControlEqn_dRate = dCurrentVolRate_dRate;
    dControlEqn_dPres = dCurrentVolRate_dPres;
  }
  else
  {
    GEOSX_ERROR( "This constraint is not supported in SinglePhaseWell" );
  }

  localRhs[eqnRowIndex] += controlEqn;
  localMatrix.addToRow< serialAtomic >( eqnRowIndex,
                                        &presDofColIndex,
                                        &dControlEqn_dPres,
                                        1 );
  localMatrix.addToRow< serialAtomic >( eqnRowIndex,
                                        &rateDofColIndex,
                                        &dControlEqn_dRate,
                                        1 );
}

/******************************** FluxKernel ********************************/

void
FluxKernel::
  launch( localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & connRate,
          arrayView1d< real64 const > const & dConnRate,
          real64 const & dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  // loop over the well elements to compute the fluxes between elements
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
  {

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
    real64 const flux = dt * currentConnRate;
    real64 const dFlux_dRate = dt;

    // 2) Assemble the flux into residual and Jacobian
    if( iwelemNext < 0 )
    {
      // flux terms
      real64 const oneSidedLocalFlux = -flux;
      real64 const oneSidedLocalFluxJacobian_dRate = -dFlux_dRate;

      // jacobian indices
      globalIndex const oneSidedEqnRowIndex = wellElemDofNumber[iwelem] + ROFFSET::MASSBAL - rankOffset;
      globalIndex const oneSidedDofColIndex_dRate = wellElemDofNumber[iwelem] + COFFSET::DRATE;

      if( oneSidedEqnRowIndex >= 0 && oneSidedEqnRowIndex < localMatrix.numRows() )
      {
        localMatrix.addToRow< parallelDeviceAtomic >( oneSidedEqnRowIndex,
                                                      &oneSidedDofColIndex_dRate,
                                                      &oneSidedLocalFluxJacobian_dRate,
                                                      1 );
        atomicAdd( parallelDeviceAtomic{}, &localRhs[oneSidedEqnRowIndex], oneSidedLocalFlux );
      }
    }
    else
    {
      // local working variables and arrays
      globalIndex eqnRowIndices[2]{};

      real64 localFlux[2]{};
      real64 localFluxJacobian_dRate[2]{};

      // flux terms
      localFlux[TAG::NEXT]    =   flux;
      localFlux[TAG::CURRENT] = -flux;

      localFluxJacobian_dRate[TAG::NEXT]    =   dFlux_dRate;
      localFluxJacobian_dRate[TAG::CURRENT] = -dFlux_dRate;

      // indices
      eqnRowIndices[TAG::CURRENT] = wellElemDofNumber[iwelem] + ROFFSET::MASSBAL - rankOffset;
      eqnRowIndices[TAG::NEXT]    = wellElemDofNumber[iwelemNext] + ROFFSET::MASSBAL - rankOffset;
      globalIndex const dofColIndex_dRate = wellElemDofNumber[iwelem] + COFFSET::DRATE;

      for( localIndex i = 0; i < 2; ++i )
      {
        if( eqnRowIndices[i] >= 0 && eqnRowIndices[i] < localMatrix.numRows() )
        {
          localMatrix.addToRow< parallelDeviceAtomic >( eqnRowIndices[i],
                                                        &dofColIndex_dRate,
                                                        &localFluxJacobian_dRate[i],
                                                        1 );
          atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndices[i]], localFlux[i] );
        }
      }
    }
  } );
}

/******************************** PressureRelationKernel ********************************/

localIndex
PressureRelationKernel::
  launch( localIndex const size,
          globalIndex const rankOffset,
          bool const isLocallyOwned,
          localIndex const iwelemControl,
          WellControls const & wellControls,
          real64 const & timeAtEndOfStep,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & wellElemPressure,
          arrayView1d< real64 const > const & dWellElemPressure,
          arrayView2d< real64 const > const & wellElemDensity,
          arrayView2d< real64 const > const & dWellElemDensity_dPres,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  // static well control data
  WellControls::Type const wellType = wellControls.getType();
  WellControls::Control const currentControl = wellControls.getControl();
  real64 const targetBHP = wellControls.getTargetBHP( timeAtEndOfStep );
  real64 const targetRate = wellControls.getTargetTotalRate( timeAtEndOfStep );

  // dynamic well control data
  real64 const & currentBHP =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentBHPString() );
  real64 const & dCurrentBHP_dPres =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::dCurrentBHP_dPresString() );

  real64 const & currentVolRate =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentVolRateString() );
  real64 const & dCurrentVolRate_dPres =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::dCurrentVolRate_dPresString() );
  real64 const & dCurrentVolRate_dRate =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::dCurrentVolRate_dRateString() );

  RAJA::ReduceMax< parallelDeviceReduce, localIndex > switchControl( 0 );

  // loop over the well elements to compute the pressure relations between well elements
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
  {

    localIndex const iwelemNext = nextWellElemIndex[iwelem];

    if( iwelemNext < 0 && isLocallyOwned ) // if iwelemNext < 0, form control equation
    {
      WellControls::Control newControl = currentControl;
      ControlEquationHelper::switchControl( wellType,
                                            currentControl,
                                            targetBHP,
                                            targetRate,
                                            currentBHP,
                                            currentVolRate,
                                            newControl );
      if( currentControl != newControl )
      {
        switchControl.max( 1 );
      }

      ControlEquationHelper::compute( rankOffset,
                                      newControl,
                                      targetBHP,
                                      targetRate,
                                      currentBHP,
                                      dCurrentBHP_dPres,
                                      currentVolRate,
                                      dCurrentVolRate_dPres,
                                      dCurrentVolRate_dRate,
                                      wellElemDofNumber[iwelemControl],
                                      localMatrix,
                                      localRhs );
    }
    else if( iwelemNext >= 0 )  // if iwelemNext >= 0, form momentum equation
    {

      // local working variables and arrays
      globalIndex dofColIndices[2]{};
      real64 localPresRelJacobian[2]{};

      // compute avg density
      real64 const avgDensity = 0.5 * ( wellElemDensity[iwelem][0] + wellElemDensity[iwelemNext][0] );
      real64 const dAvgDensity_dPresNext    = 0.5 * dWellElemDensity_dPres[iwelemNext][0];
      real64 const dAvgDensity_dPresCurrent = 0.5 * dWellElemDensity_dPres[iwelem][0];

      // compute depth diff times acceleration
      real64 const gravD = wellElemGravCoef[iwelemNext] - wellElemGravCoef[iwelem];

      // compute the current pressure in the two well elements
      real64 const pressureCurrent = wellElemPressure[iwelem] + dWellElemPressure[iwelem];
      real64 const pressureNext    = wellElemPressure[iwelemNext] + dWellElemPressure[iwelemNext];

      // compute momentum flux and derivatives
      real64 const localPresRel = pressureNext - pressureCurrent - avgDensity * gravD;
      localPresRelJacobian[TAG::NEXT]    =  1 - dAvgDensity_dPresNext * gravD;
      localPresRelJacobian[TAG::CURRENT] = -1 - dAvgDensity_dPresCurrent * gravD;

      // TODO: add friction and acceleration terms

      // jacobian indices
      globalIndex const eqnRowIndex = wellElemDofNumber[iwelem] + ROFFSET::CONTROL - rankOffset;
      dofColIndices[TAG::NEXT]      = wellElemDofNumber[iwelemNext] + COFFSET::DPRES;
      dofColIndices[TAG::CURRENT]   = wellElemDofNumber[iwelem] + COFFSET::DPRES;

      if( eqnRowIndex >= 0 && eqnRowIndex < localMatrix.numRows() )
      {
        localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndex,
                                                                          &dofColIndices[0],
                                                                          &localPresRelJacobian[0],
                                                                          2 );
        atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndex], localPresRel );
      }
    }
  } );
  return switchControl.get();
}

/******************************** PerforationKernel ********************************/

GEOSX_HOST_DEVICE
void
PerforationKernel::
  compute( real64 const & resPressure,
           real64 const & dResPressure,
           real64 const & resDensity,
           real64 const & dResDensity_dPres,
           real64 const & resViscosity,
           real64 const & dResViscosity_dPres,
           real64 const & wellElemGravCoef,
           real64 const & wellElemPressure,
           real64 const & dWellElemPressure,
           real64 const & wellElemDensity,
           real64 const & dWellElemDensity_dPres,
           real64 const & wellElemViscosity,
           real64 const & dWellElemViscosity_dPres,
           real64 const & perfGravCoef,
           real64 const & trans,
           real64 & perfRate,
           arraySlice1d< real64 > const & dPerfRate_dPres )
{
  // local working variables and arrays
  real64 pressure[2]{};
  real64 dPressure_dP[2]{};
  real64 multiplier[2]{};

  perfRate = 0.0;
  dPerfRate_dPres[0] = 0.0;
  dPerfRate_dPres[1] = 0.0;

  // 1) Reservoir side


  // get reservoir variables
  pressure[TAG::RES] = resPressure + dResPressure;
  dPressure_dP[TAG::RES] = 1;

  // TODO: add a buoyancy term for the reservoir side here

  // multiplier for reservoir side in the flux
  multiplier[TAG::RES] = 1;

  // 2) Well side


  // get well variables
  pressure[TAG::WELL] = wellElemPressure + dWellElemPressure;
  dPressure_dP[TAG::WELL] = 1.0;

  real64 const gravD = ( perfGravCoef - wellElemGravCoef );
  pressure[TAG::WELL]     += wellElemDensity * gravD;
  dPressure_dP[TAG::WELL] += dWellElemDensity_dPres * gravD;

  // multiplier for well side in the flux
  multiplier[TAG::WELL] = -1;

  // compute potential difference
  real64 potDif = 0.0;
  for( localIndex i = 0; i < 2; ++i )
  {
    potDif += multiplier[i] * trans * pressure[i];
    dPerfRate_dPres[i] = multiplier[i] * trans * dPressure_dP[i];
  }

  // choose upstream cell based on potential difference
  localIndex const k_up = (potDif >= 0) ? TAG::RES : TAG::WELL;

  // compute upstream density, viscosity, and mobility
  real64 densityUp       = 0.0;
  real64 dDensityUp_dP   = 0.0;
  real64 viscosityUp     = 0.0;
  real64 dViscosityUp_dP = 0.0;

  // upwinding the variables
  if( k_up == TAG::RES ) // use reservoir vars
  {
    densityUp     = resDensity;
    dDensityUp_dP = dResDensity_dPres;

    viscosityUp     = resViscosity;
    dViscosityUp_dP = dResViscosity_dPres;
  }
  else // use well vars
  {
    densityUp = wellElemDensity;
    dDensityUp_dP = dWellElemDensity_dPres;

    viscosityUp = wellElemViscosity;
    dViscosityUp_dP = dWellElemViscosity_dPres;
  }

  // compute mobility
  real64 const mobilityUp     = densityUp / viscosityUp;
  real64 const dMobilityUp_dP = dDensityUp_dP / viscosityUp
                                - mobilityUp / viscosityUp * dViscosityUp_dP;

  perfRate = mobilityUp * potDif;
  for( localIndex ke = 0; ke < 2; ++ke )
  {
    dPerfRate_dPres[ke] *= mobilityUp;
  }
  dPerfRate_dPres[k_up] += dMobilityUp_dP * potDif;
}

void
PerforationKernel::
  launch( localIndex const size,
          ElementViewConst< arrayView1d< real64 const > > const & resPressure,
          ElementViewConst< arrayView1d< real64 const > > const & dResPressure,
          ElementViewConst< arrayView2d< real64 const > > const & resDensity,
          ElementViewConst< arrayView2d< real64 const > > const & dResDensity_dPres,
          ElementViewConst< arrayView2d< real64 const > > const & resViscosity,
          ElementViewConst< arrayView2d< real64 const > > const & dResViscosity_dPres,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< real64 const > const & wellElemPressure,
          arrayView1d< real64 const > const & dWellElemPressure,
          arrayView2d< real64 const > const & wellElemDensity,
          arrayView2d< real64 const > const & dWellElemDensity_dPres,
          arrayView2d< real64 const > const & wellElemViscosity,
          arrayView2d< real64 const > const & dWellElemViscosity_dPres,
          arrayView1d< real64 const > const & perfGravCoef,
          arrayView1d< localIndex const > const & perfWellElemIndex,
          arrayView1d< real64 const > const & perfTransmissibility,
          arrayView1d< localIndex const > const & resElementRegion,
          arrayView1d< localIndex const > const & resElementSubRegion,
          arrayView1d< localIndex const > const & resElementIndex,
          arrayView1d< real64 > const & perfRate,
          arrayView2d< real64 > const & dPerfRate_dPres )
{

  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const iperf )
  {

    // get the reservoir (sub)region and element indices
    localIndex const er  = resElementRegion[iperf];
    localIndex const esr = resElementSubRegion[iperf];
    localIndex const ei  = resElementIndex[iperf];

    // get the local index of the well element
    localIndex const iwelem = perfWellElemIndex[iperf];

    compute( resPressure[er][esr][ei],
             dResPressure[er][esr][ei],
             resDensity[er][esr][ei][0],
             dResDensity_dPres[er][esr][ei][0],
             resViscosity[er][esr][ei][0],
             dResViscosity_dPres[er][esr][ei][0],
             wellElemGravCoef[iwelem],
             wellElemPressure[iwelem],
             dWellElemPressure[iwelem],
             wellElemDensity[iwelem][0],
             dWellElemDensity_dPres[iwelem][0],
             wellElemViscosity[iwelem][0],
             dWellElemViscosity_dPres[iwelem][0],
             perfGravCoef[iperf],
             perfTransmissibility[iperf],
             perfRate[iperf],
             dPerfRate_dPres[iperf] );


  } );
}


/******************************** AccumulationKernel ********************************/

void
AccumulationKernel::
  launch( localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< integer const > const & wellElemGhostRank,
          arrayView1d< real64 const > const & wellElemVolume,
          arrayView2d< real64 const > const & wellElemDensity,
          arrayView2d< real64 const > const & dWellElemDensity_dPres,
          arrayView1d< real64 const > const & wellElemDensityOld,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
  {

    if( wellElemGhostRank[iwelem] >= 0 )
    {
      return;
    }

    localIndex const eqnRowIndex = wellElemDofNumber[iwelem] + ROFFSET::MASSBAL - rankOffset;
    globalIndex const presDofColIndex = wellElemDofNumber[iwelem] + COFFSET::DPRES;

    real64 const localAccum = wellElemVolume[iwelem] * ( wellElemDensity[iwelem][0] - wellElemDensityOld[iwelem] );
    real64 const localAccumJacobian = wellElemVolume[iwelem] * dWellElemDensity_dPres[iwelem][0];

    // add contribution to global residual and jacobian (no need for atomics here)
    localMatrix.addToRow< serialAtomic >( eqnRowIndex, &presDofColIndex, &localAccumJacobian, 1 );
    localRhs[eqnRowIndex] += localAccum;
  } );
}


/******************************** PressureInitializationKernel ********************************/

void
PresInitializationKernel::
  launch( localIndex const perforationSize,
          localIndex const subRegionSize,
          localIndex const numPerforations,
          WellControls const & wellControls,
          real64 const & currentTime,
          ElementViewConst< arrayView1d< real64 const > > const & resPressure,
          ElementViewConst< arrayView2d< real64 const > > const & resDensity,
          arrayView1d< localIndex const > const & resElementRegion,
          arrayView1d< localIndex const > const & resElementSubRegion,
          arrayView1d< localIndex const > const & resElementIndex,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< real64 > const & wellElemPressure )
{
  real64 const targetBHP = wellControls.getTargetBHP( currentTime );
  real64 const refWellElemGravCoef = wellControls.getReferenceGravityCoef();
  WellControls::Control const currentControl = wellControls.getControl();
  WellControls::Type const wellType = wellControls.getType();

  // loop over all perforations to compute an average density
  RAJA::ReduceSum< parallelDeviceReduce, real64 > sumDensity( 0 );
  RAJA::ReduceMin< parallelDeviceReduce, real64 > minResPressure( 1e10 );
  RAJA::ReduceMax< parallelDeviceReduce, real64 > maxResPressure( 0 );
  forAll< parallelDevicePolicy<> >( perforationSize, [=] GEOSX_HOST_DEVICE ( localIndex const iperf )
  {

    // get the reservoir (sub)region and element indices
    localIndex const er = resElementRegion[iperf];
    localIndex const esr = resElementSubRegion[iperf];
    localIndex const ei = resElementIndex[iperf];

    sumDensity += resDensity[er][esr][ei][0];
    minResPressure.min( resPressure[er][esr][ei] );
    maxResPressure.max( resPressure[er][esr][ei] );
  } );
  real64 const pres = ( wellControls.getType() == WellControls::Type::PRODUCER )
                    ? MpiWrapper::min( minResPressure.get() )
                    : MpiWrapper::max( maxResPressure.get() );
  real64 const avgDensity = MpiWrapper::sum( sumDensity.get() ) / numPerforations;

  real64 pressureControl = 0.0;
  real64 const gravCoefControl = refWellElemGravCoef;
  // initialize the pressure in the element where the BHP is controlled)
  if( currentControl == WellControls::Control::BHP )
  {
    // if pressure constraint, initialize the pressure at the constraint
    pressureControl = targetBHP;
  }
  else // rate control
  {
    // initialize the pressure in the element where the BHP is controlled slightly
    // above/below the target pressure depending on well type.
    // note: the targetBHP is not used here because we sometimes set targetBHP to a very large (unrealistic) value
    //       to keep the well in rate control during the full simulation, and we don't want this large targetBHP to
    //       be used for initialization
    pressureControl = ( wellType == WellControls::Type::PRODUCER )
                    ? 0.5 * pres
                    : 2.0 * pres;
  }

  GEOSX_ERROR_IF( pressureControl <= 0, "Invalid well initialization: negative pressure was found" );

  // estimate the pressures in the well elements using this avgDensity
  forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
  {
    wellElemPressure[iwelem] = pressureControl
                               + avgDensity * ( wellElemGravCoef[iwelem] - gravCoefControl );
  } );
}

/******************************** RateInitializationKernel ********************************/

void
RateInitializationKernel::
  launch( localIndex const subRegionSize,
          WellControls const & wellControls,
          real64 const & currentTime,
          arrayView2d< real64 const > const & wellElemDens,
          arrayView1d< real64 > const & connRate )
{
  real64 const targetRate = wellControls.getTargetTotalRate( currentTime );
  WellControls::Control const control = wellControls.getControl();
  WellControls::Type const wellType = wellControls.getType();

  // Estimate the connection rates
  forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
  {
    if( control == WellControls::Control::BHP )
    {
      // if BHP constraint set rate below the absolute max rate
      // with the appropriate sign (negative for prod, positive for inj)
      if( wellType == WellControls::Type::PRODUCER )
      {
        connRate[iwelem] = LvArray::math::max( 0.1 * targetRate * wellElemDens[iwelem][0], -1e3 );
      }
      else
      {
        connRate[iwelem] = LvArray::math::min( 0.1 * targetRate * wellElemDens[iwelem][0], 1e3 );
      }
    }
    else
    {
      connRate[iwelem] = targetRate * wellElemDens[iwelem][0];
    }
  } );
}


} // end namespace SinglePhaseWellKernels

} // end namespace geosx
