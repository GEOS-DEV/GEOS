/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
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
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWellFields.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidLayouts.hpp"
namespace geos
{

namespace singlePhaseWellKernels
{

/******************************** ControlEquationHelper ********************************/

GEOS_HOST_DEVICE
inline
void
ControlEquationHelper::
  switchControl( bool const isProducer,
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
    if( isProducer )
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

template< integer IS_THERMAL >
GEOS_HOST_DEVICE
inline
void
ControlEquationHelper::
  compute( globalIndex const rankOffset,
           WellControls::Control const currentControl,
           real64 const & targetBHP,
           real64 const & targetRate,
           real64 const & currentBHP,
           arrayView1d< real64 const > const & dCurrentBHP,
           real64 const & currentVolRate,
           arrayView1d< real64 const > const & dCurrentTotalVolRate,
           globalIndex const dofNumber,
           CRSMatrixView< real64, globalIndex const > const & localMatrix,
           arrayView1d< real64 > const & localRhs )
{
  using ROFFSET_WJ = singlePhaseWellKernels::RowOffset_WellJac< IS_THERMAL >;
  using COFFSET_WJ = singlePhaseWellKernels::ColOffset_WellJac< IS_THERMAL >;
  using Deriv = constitutive::singlefluid::DerivativeOffsetC< IS_THERMAL >;

  localIndex const eqnRowIndex = dofNumber + ROFFSET_WJ::CONTROL - rankOffset;
  globalIndex dofColIndices[COFFSET_WJ::nDer]{};
  for( integer i = 0; i < COFFSET_WJ::nDer; ++i )
  {
    dofColIndices[ i ] = dofNumber + i;
  }

  real64 controlEqn = 0;
  real64 dControlEqn[2+IS_THERMAL]{};

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
    dControlEqn[COFFSET_WJ::dP] = dCurrentBHP[Deriv::dP];
    if constexpr ( IS_THERMAL )
      dControlEqn[COFFSET_WJ::dT] = dCurrentBHP[Deriv::dT];

  }
  // Total volumetric rate control
  else if( currentControl == WellControls::Control::TOTALVOLRATE )
  {
    // control equation is the difference between volumetric current rate and target rate
    controlEqn = currentVolRate - targetRate;
    dControlEqn[COFFSET_WJ::dP] = dCurrentTotalVolRate[COFFSET_WJ::dP];
    dControlEqn[COFFSET_WJ::dQ] = dCurrentTotalVolRate[COFFSET_WJ::dQ];
    if constexpr ( IS_THERMAL )
      dControlEqn[COFFSET_WJ::dT] = dCurrentTotalVolRate[Deriv::dT];
  }
  else
  {
    GEOS_ERROR( "This constraint is not supported in SinglePhaseWell" );
  }

  localRhs[eqnRowIndex] += controlEqn;
  localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowIndex,
                                                            dofColIndices,
                                                            dControlEqn,
                                                            COFFSET_WJ::nDer );
}

/******************************** FluxKernel ********************************/
#define INST_FluxKernel( IS_THERMAL ) \
  template \
  void  \
  FluxKernel::  \
    launch< IS_THERMAL >( localIndex const size,  \
                          globalIndex const rankOffset,  \
                          arrayView1d< globalIndex const > const & wellElemDofNumber,  \
                          arrayView1d< localIndex const > const & nextWellElemIndex,  \
                          arrayView1d< real64 const > const & connRate,  \
                          real64 const & dt,  \
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,  \
                          arrayView1d< real64 > const & localRhs )

INST_FluxKernel( 0 );
INST_FluxKernel( 1 );

template< integer IS_THERMAL >
void
FluxKernel::
  launch( localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & connRate,
          real64 const & dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  // loop over the well elements to compute the fluxes between elements
  forAll< parallelDevicePolicy<> >( size, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
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
    real64 const currentConnRate = connRate[iwelem];
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
        RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[oneSidedEqnRowIndex], oneSidedLocalFlux );
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
          RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndices[i]], localFlux[i] );
        }
      }
    }
  } );
}

/******************************** PressureRelationKernel ********************************/

#define INST_PressureRelationKernel( IS_THERMAL ) \
  template \
  localIndex \
  PressureRelationKernel:: \
    launch< IS_THERMAL >( localIndex const size, \
                          globalIndex const rankOffset, \
                          bool const isLocallyOwned, \
                          localIndex const iwelemControl, \
                          WellControls const & wellControls, \
                          real64 const & timeAtEndOfStep, \
                          arrayView1d< globalIndex const > const & wellElemDofNumber, \
                          arrayView1d< real64 const > const & wellElemGravCoef, \
                          arrayView1d< localIndex const > const & nextWellElemIndex, \
                          arrayView1d< real64 const > const & wellElemPressure, \
                          arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & wellElemDensity, \
                          arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const & dWellElemDensity, \
                          CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                          arrayView1d< real64 > const & localRhs )

INST_PressureRelationKernel( 0 );
INST_PressureRelationKernel( 1 );

template< integer IS_THERMAL >
localIndex
PressureRelationKernel::
  launch( localIndex const size,
          globalIndex const rankOffset,
          bool const isLocallyOwned,
          localIndex const iwelemControl,
          WellControls const & wellControls,
          real64 const & time,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & wellElemPressure,
          arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & wellElemDensity,
          arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const & dWellElemDensity,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  using Deriv = constitutive::singlefluid::DerivativeOffset;
  using COFFSET_WJ = singlePhaseWellKernels::ColOffset_WellJac< IS_THERMAL >;
  // static well control data
  bool const isProducer = wellControls.isProducer();
  WellControls::Control const currentControl = wellControls.getControl();
  real64 const targetBHP = wellControls.getTargetBHP( time );
  real64 const targetRate = wellControls.getTargetTotalRate( time );

  // dynamic well control data
  real64 const & currentBHP =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentBHPString() );
  arrayView1d< real64 const > const & dCurrentBHP =
    wellControls.getReference< array1d< real64 > >( SinglePhaseWell::viewKeyStruct::dCurrentBHPString() );

  real64 const & currentVolRate =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentVolRateString() );
  arrayView1d< real64 const > const & dCurrentVolRate =
    wellControls.getReference< array1d< real64 > >( SinglePhaseWell::viewKeyStruct::dCurrentVolRateString() );
  
  RAJA::ReduceMax< parallelDeviceReduce, localIndex > switchControl( 0 );

  // loop over the well elements to compute the pressure relations between well elements
  forAll< parallelDevicePolicy<> >( size, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
  {

    localIndex const iwelemNext = nextWellElemIndex[iwelem];

    if( iwelemNext < 0 && isLocallyOwned ) // if iwelemNext < 0, form control equation
    {
      WellControls::Control newControl = currentControl;
      ControlEquationHelper::switchControl( isProducer,
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

      ControlEquationHelper::compute< IS_THERMAL >( rankOffset,
                                                    newControl,
                                                    targetBHP,
                                                    targetRate,
                                                    currentBHP,
                                                    dCurrentBHP,
                                                    currentVolRate,
                                                    dCurrentVolRate,
                                                    wellElemDofNumber[iwelemControl],
                                                    localMatrix,
                                                    localRhs );
    }
    else if( iwelemNext >= 0 )  // if iwelemNext >= 0, form momentum equation
    {

      // local working variables and arrays
      globalIndex dofColIndices[2*(1+IS_THERMAL)]{};
      real64 localPresRelJacobian[2*(1+IS_THERMAL)]{};

      // compute avg density
      real64 const avgDensity = 0.5 * ( wellElemDensity[iwelem][0] + wellElemDensity[iwelemNext][0] );
      real64 const dAvgDensity_dPresNext    = 0.5 * dWellElemDensity[iwelemNext][0][Deriv::dP];
      real64 const dAvgDensity_dPresCurrent = 0.5 * dWellElemDensity[iwelem][0][Deriv::dP];

      // compute depth diff times acceleration
      real64 const gravD = wellElemGravCoef[iwelemNext] - wellElemGravCoef[iwelem];

      // compute the current pressure in the two well elements
      real64 const pressureCurrent = wellElemPressure[iwelem];
      real64 const pressureNext    = wellElemPressure[iwelemNext];

      // compute momentum flux and derivatives
      real64 const localPresRel = pressureNext - pressureCurrent - avgDensity * gravD;
      localPresRelJacobian[TAG::NEXT *(1+IS_THERMAL)]    =  1 - dAvgDensity_dPresNext * gravD;
      localPresRelJacobian[TAG::CURRENT *(1+IS_THERMAL)]  = -1 - dAvgDensity_dPresCurrent * gravD;

      if constexpr ( IS_THERMAL )
      {
        localPresRelJacobian[TAG::NEXT *(1+IS_THERMAL)+1]    =  -0.5 * dWellElemDensity[iwelemNext][0][Deriv::dT]* gravD;
        localPresRelJacobian[TAG::CURRENT *(1+IS_THERMAL)+1] =  -0.5 * dWellElemDensity[iwelem][0][Deriv::dT]* gravD;
      }

      // TODO: add friction and acceleration terms

      // jacobian indices
      globalIndex const eqnRowIndex = wellElemDofNumber[iwelem] + ROFFSET::CONTROL - rankOffset;
      dofColIndices[TAG::NEXT *(1+IS_THERMAL)]      = wellElemDofNumber[iwelemNext] + COFFSET_WJ::dP;
      dofColIndices[TAG::CURRENT *(1+IS_THERMAL)]   = wellElemDofNumber[iwelem] + COFFSET_WJ::dP;

      if constexpr ( IS_THERMAL )
      {
        dofColIndices[TAG::NEXT *(1+IS_THERMAL)+1]    = wellElemDofNumber[iwelemNext] + COFFSET_WJ::dT;
        dofColIndices[TAG::CURRENT *(1+IS_THERMAL)+1] = wellElemDofNumber[iwelem] + COFFSET_WJ::dT;
      }
      if( eqnRowIndex >= 0 && eqnRowIndex < localMatrix.numRows() )
      {
        localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndex,
                                                                          dofColIndices,
                                                                          localPresRelJacobian,
                                                                          2 * (1+IS_THERMAL) );
        RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndex], localPresRel );
      }
    }
  } );
  return switchControl.get();
}

/******************************** AccumulationKernel ********************************/

void
AccumulationKernel::
  launch( localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< integer const > const & wellElemGhostRank,
          arrayView1d< real64 const > const & wellElemVolume,
          arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & wellElemDensity,
          arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const & dWellElemDensity,
          arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & wellElemDensity_n,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  using Deriv = constitutive::singlefluid::DerivativeOffset;
  forAll< parallelDevicePolicy<> >( size, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
  {

    if( wellElemGhostRank[iwelem] >= 0 )
    {
      return;
    }

    localIndex const eqnRowIndex = wellElemDofNumber[iwelem] + ROFFSET::MASSBAL - rankOffset;
    globalIndex const presDofColIndex = wellElemDofNumber[iwelem] + COFFSET::DPRES;

    real64 const localAccum = wellElemVolume[iwelem] * ( wellElemDensity[iwelem][0] - wellElemDensity_n[iwelem][0] );
    real64 const localAccumJacobian = wellElemVolume[iwelem] * dWellElemDensity[iwelem][0][Deriv::dP];

    // add contribution to global residual and jacobian (no need for atomics here)
    localMatrix.addToRow< serialAtomic >( eqnRowIndex, &presDofColIndex, &localAccumJacobian, 1 );
    localRhs[eqnRowIndex] += localAccum;
  } );
}


/******************************** PressureInitializationKernel ********************************/

void
PresTempInitializationKernel::
  launch( integer const isThermal,
          localIndex const perforationSize,
          localIndex const subRegionSize,
          localIndex const numPerforations,
          WellControls const & wellControls,
          real64 const & currentTime,
          ElementViewConst< arrayView1d< real64 const > > const & resPressure,
          ElementViewConst< arrayView1d< real64 const > > const & resTemp,
          ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const & resDensity,
          arrayView1d< localIndex const > const & resElementRegion,
          arrayView1d< localIndex const > const & resElementSubRegion,
          arrayView1d< localIndex const > const & resElementIndex,
          arrayView1d< real64 const > const & perfGravCoef,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< real64 > const & wellElemPressure,
          arrayView1d< real64 > const & wellElemTemperature )
{
  real64 const targetBHP = wellControls.getTargetBHP( currentTime );
  real64 const refWellElemGravCoef = wellControls.getReferenceGravityCoef();
  real64 const initialPressureCoef = wellControls.getInitialPressureCoefficient();
  WellControls::Control const currentControl = wellControls.getControl();
  bool const isProducer = wellControls.isProducer();



  // Step 1: we loop over all the perforations on this rank to compute the following quantities:
  //   - Sum of densities over the perforated reservoir elements
  // In passing, we save the min gravCoef difference between the reference depth and the perforation depth
  // Note that we use gravCoef instead of depth for the (unlikely) case in which the gravityVector is not aligned with z

  RAJA::ReduceSum< parallelDeviceReduce, real64 > sumDensity( 0 );
  RAJA::ReduceSum< parallelDeviceReduce, real64 > sumTemp( 0 );
  RAJA::ReduceMin< parallelDeviceReduce, real64 > localMinGravCoefDiff( 1e9 );

  forAll< parallelDevicePolicy<> >( perforationSize, [=] GEOS_HOST_DEVICE ( localIndex const iperf )
  {
    // get the reservoir (sub)region and element indices
    localIndex const er = resElementRegion[iperf];
    localIndex const esr = resElementSubRegion[iperf];
    localIndex const ei = resElementIndex[iperf];

    // save the min gravCoef difference between the reference depth and the perforation depth (times g)
    localMinGravCoefDiff.min( LvArray::math::abs( refWellElemGravCoef - perfGravCoef[iperf] ) );

    // increment the fluid density
    sumDensity += resDensity[er][esr][ei][0];

    // increment the temperature
    sumTemp += resTemp[er][esr][ei];
  } );

  real64 const minGravCoefDiff = MpiWrapper::min( localMinGravCoefDiff.get() );



  // Step 2: we assign average quantities over the well (i.e., over all the ranks)

  real64 const avgDensity = MpiWrapper::sum( sumDensity.get() ) / numPerforations;

  // Step 3: we compute the approximate pressure at the reference depth
  // We make a distinction between pressure-controlled wells and rate-controlled wells

  real64 refPres = 0.0;
  real64 avgTemp = 0;


  // for a producer, we use the temperature and component fractions from the reservoir
  if( isProducer )
  {
    // use average temperature from reservoir
    avgTemp = MpiWrapper::sum( sumTemp.get() ) / numPerforations;
  }
  // for an injector, we use the injection stream values
  else
  {
    // use temperature from injection stream
    avgTemp = wellControls.getInjectionTemperature();
  }
  // if the well is controlled by pressure, initialize the reference pressure at the target pressure
  if( currentControl == WellControls::Control::BHP )
  {
    refPres = targetBHP;
  }
  // if the well is controlled by rate, initialize the reference pressure using the pressure at the closest perforation
  else
  {
    RAJA::ReduceMin< parallelDeviceReduce, real64 > localRefPres( 1e9 );
    real64 const alpha = ( isProducer ) ? 1 - initialPressureCoef : 1 + initialPressureCoef;

    forAll< parallelDevicePolicy<> >( perforationSize, [=] GEOS_HOST_DEVICE ( localIndex const iperf )
    {
      // get the reservoir (sub)region and element indices
      localIndex const er = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei = resElementIndex[iperf];

      // get the perforation pressure and save the estimated reference pressure
      real64 const gravCoefDiff = LvArray::math::abs( refWellElemGravCoef - perfGravCoef[iperf] );
      if( isZero( gravCoefDiff - minGravCoefDiff ) )
      {
        localRefPres.min( alpha * resPressure[er][esr][ei] + avgDensity * ( refWellElemGravCoef - perfGravCoef[iperf] ) );
      }
    } );
    refPres = MpiWrapper::min( localRefPres.get() );
  }



  // Step 4: we are ready to assign the primary variables on the well elements:
  //  - pressure: hydrostatic pressure using our crude approximation of the total mass density

  RAJA::ReduceMax< parallelDeviceReduce, integer > foundNegativePressure( 0 );
  RAJA::ReduceMax< parallelDeviceReduce, integer > foundNegativeTemp( 0 );

  forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
  {
    wellElemPressure[iwelem] = refPres + avgDensity * ( wellElemGravCoef[iwelem] - refWellElemGravCoef );
    wellElemTemperature[iwelem] = avgTemp;
    if( wellElemPressure[iwelem] <= 0 )
    {
      foundNegativePressure.max( 1 );
    }
    if( wellElemTemperature[iwelem] < 0 )
    {
      foundNegativeTemp.max( 1 );
    }
  } );

  GEOS_THROW_IF( foundNegativePressure.get() == 1,
                 wellControls.getDataContext() << ": Invalid well initialization, negative pressure was found.",
                 InputError );
  if( isThermal )   // tjb change  temp in isothermal cases shouldnt be an issue (also what if temp in fluid prop calcs like compo)
  {
    GEOS_THROW_IF( foundNegativeTemp.get() == 1,
                   wellControls.getDataContext() << "Invalid well initialization, negative temperature was found.",
                   InputError );
  }
}

/******************************** RateInitializationKernel ********************************/

void
RateInitializationKernel::
  launch( localIndex const subRegionSize,
          WellControls const & wellControls,
          real64 const & currentTime,
          arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & wellElemDens,
          arrayView1d< real64 > const & connRate )
{
  real64 const targetRate = wellControls.getTargetTotalRate( currentTime );
  WellControls::Control const control = wellControls.getControl();
  bool const isProducer = wellControls.isProducer();

  // Estimate the connection rates
  forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
  {
    if( control == WellControls::Control::BHP )
    {
      // if BHP constraint set rate below the absolute max rate
      // with the appropriate sign (negative for prod, positive for inj)
      if( isProducer )
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

} // end namespace singlePhaseWellKernels

} // end namespace geos
