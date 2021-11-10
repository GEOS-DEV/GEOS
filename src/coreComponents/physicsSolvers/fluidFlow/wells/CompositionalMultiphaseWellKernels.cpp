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
 * @file CompositionalMultiphaseWellKernels.cpp
 */

#include "CompositionalMultiphaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"

namespace geosx
{

namespace CompositionalMultiphaseWellKernels
{

/******************************** ControlEquationHelper ********************************/

GEOSX_HOST_DEVICE
void
ControlEquationHelper::
  switchControl( WellControls::Type const & wellType,
                 WellControls::Control const & currentControl,
                 localIndex const phasePhaseIndex,
                 real64 const & targetBHP,
                 real64 const & targetPhaseRate,
                 real64 const & targetTotalRate,
                 real64 const & currentBHP,
                 arrayView1d< real64 const > const & currentPhaseVolRate,
                 real64 const & currentTotalVolRate,
                 WellControls::Control & newControl )
{
  // if isViable is true at the end of the following checks, no need to switch
  bool controlIsViable = false;

  // The limiting flow rates are treated as upper limits, while the pressure limits
  // are treated as lower limits in production wells and upper limits in injectors.
  // The well changes its mode of control whenever the existing control mode would
  // violate one of these limits.

  // Currently, the available constraints are:
  //   - Producer: BHP, PHASEVOLRATE
  //   - Injector: BHP, TOTALVOLRATE

  // TODO: support GRAT, WRAT, LIQUID for producers and check if any of the active constraint is violated

  // BHP control
  if( currentControl == WellControls::Control::BHP )
  {
    // the control is viable if the reference oil rate is below the max rate for producers
    if( wellType == WellControls::Type::PRODUCER )
    {
      controlIsViable = ( LvArray::math::abs( currentPhaseVolRate[phasePhaseIndex] ) <= LvArray::math::abs( targetPhaseRate ) );
    }
    // the control is viable if the reference total rate is below the max rate for injectors
    else
    {
      controlIsViable = ( LvArray::math::abs( currentTotalVolRate ) <= LvArray::math::abs( targetTotalRate ) );
    }
  }
  else // rate control
  {
    // the control is viable if the reference pressure is below/above the max/min pressure
    if( wellType == WellControls::Type::PRODUCER )
    {
      // targetBHP specifies a min pressure here
      controlIsViable = ( currentBHP >= targetBHP );
    }
    else
    {
      // targetBHP specifies a max pressure here
      controlIsViable = ( currentBHP <= targetBHP );
    }
  }

  if( controlIsViable )
  {
    newControl = currentControl;
  }
  else
  {
    if( wellType == WellControls::Type::PRODUCER )
    {
      newControl = ( currentControl == WellControls::Control::BHP )
           ? WellControls::Control::PHASEVOLRATE
           : WellControls::Control::BHP;
    }
    else
    {
      newControl = ( currentControl == WellControls::Control::BHP )
                 ? WellControls::Control::TOTALVOLRATE
                 : WellControls::Control::BHP;
    }
  }
}

template< localIndex NC >
GEOSX_HOST_DEVICE
void
ControlEquationHelper::
  compute( globalIndex const rankOffset,
           WellControls::Control const currentControl,
           localIndex const targetPhaseIndex,
           real64 const & targetBHP,
           real64 const & targetPhaseRate,
           real64 const & targetTotalRate,
           real64 const & currentBHP,
           real64 const & dCurrentBHP_dPres,
           arrayView1d< real64 const > const & dCurrentBHP_dCompDens,
           arrayView1d< real64 const > const & currentPhaseVolRate,
           arrayView1d< real64 const > const & dCurrentPhaseVolRate_dPres,
           arrayView2d< real64 const > const & dCurrentPhaseVolRate_dCompDens,
           arrayView1d< real64 const > const & dCurrentPhaseVolRate_dRate,
           real64 const & currentTotalVolRate,
           real64 const & dCurrentTotalVolRate_dPres,
           arrayView1d< real64 const > const & dCurrentTotalVolRate_dCompDens,
           real64 const & dCurrentTotalVolRate_dRate,
           globalIndex const dofNumber,
           CRSMatrixView< real64, globalIndex const > const & localMatrix,
           arrayView1d< real64 > const & localRhs )
{
  localIndex const eqnRowIndex      = dofNumber + ROFFSET::CONTROL - rankOffset;
  globalIndex const presDofColIndex = dofNumber + COFFSET::DPRES;
  globalIndex const rateDofColIndex = dofNumber + COFFSET::DCOMP + NC;

  globalIndex compDofColIndices[NC]{};
  for( localIndex ic = 0; ic < NC; ++ic )
  {
    compDofColIndices[ ic ] = presDofColIndex + ic + 1;
  }

  real64 controlEqn = 0;
  real64 dControlEqn_dPres = 0;
  real64 dControlEqn_dRate = 0;
  real64 dControlEqn_dComp[NC]{};

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
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dControlEqn_dComp[ic] = dCurrentBHP_dCompDens[ic];
    }
  }
  // Oil volumetric rate control
  else if( currentControl == WellControls::Control::PHASEVOLRATE )
  {
    controlEqn = currentPhaseVolRate[targetPhaseIndex] - targetPhaseRate;
    dControlEqn_dPres = dCurrentPhaseVolRate_dPres[targetPhaseIndex];
    dControlEqn_dRate = dCurrentPhaseVolRate_dRate[targetPhaseIndex];
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dControlEqn_dComp[ic] = dCurrentPhaseVolRate_dCompDens[targetPhaseIndex][ic];
    }
  }
  // Total volumetric rate control
  else if( currentControl == WellControls::Control::TOTALVOLRATE )
  {
    controlEqn = currentTotalVolRate - targetTotalRate;
    dControlEqn_dPres = dCurrentTotalVolRate_dPres;
    dControlEqn_dRate = dCurrentTotalVolRate_dRate;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dControlEqn_dComp[ic] = dCurrentTotalVolRate_dCompDens[ic];
    }
  }
  else
  {
    GEOSX_ERROR( "This constraint is not supported in CompositionalMultiphaseWell" );
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
  localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowIndex,
                                                            compDofColIndices,
                                                            dControlEqn_dComp,
                                                            NC );
}

/******************************** FluxKernel ********************************/

template< localIndex NC >
GEOSX_HOST_DEVICE
void
FluxKernel::
  computeExit( real64 const & dt,
               real64 const ( &compFlux )[NC],
               real64 const ( &dCompFlux_dRate )[NC],
               real64 const ( &dCompFlux_dPresUp )[NC],
               real64 const ( &dCompFlux_dCompDensUp )[NC][NC],
               real64 ( & oneSidedFlux )[NC],
               real64 ( & oneSidedFluxJacobian_dRate )[NC][1],
               real64 ( & oneSidedFluxJacobian_dPresCompUp )[NC][NC + 1] )
{
  for( localIndex ic = 0; ic < NC; ++ic )
  {
    oneSidedFlux[ic] = -dt * compFlux[ic];

    // derivative with respect to rate
    oneSidedFluxJacobian_dRate[ic][0] = -dt * dCompFlux_dRate[ic];

    // derivative with respect to upstream pressure
    oneSidedFluxJacobian_dPresCompUp[ic][0] = -dt * dCompFlux_dPresUp[ic];

    // derivatives with respect to upstream component densities
    for( localIndex jdof = 0; jdof < NC; ++jdof )
    {
      oneSidedFluxJacobian_dPresCompUp[ic][jdof+1] = -dt * dCompFlux_dCompDensUp[ic][jdof];
    }
  }
}

template< localIndex NC >
GEOSX_HOST_DEVICE
void
FluxKernel::
  compute( real64 const & dt,
           real64 const ( &compFlux )[NC],
           real64 const ( &dCompFlux_dRate )[NC],
           real64 const ( &dCompFlux_dPresUp )[NC],
           real64 const ( &dCompFlux_dCompDensUp )[NC][NC],
           real64 ( & localFlux )[2*NC],
           real64 ( & localFluxJacobian_dRate )[2*NC][1],
           real64 ( & localFluxJacobian_dPresCompUp )[2*NC][NC + 1] )
{
  // flux terms
  for( localIndex ic = 0; ic < NC; ++ic )
  {
    localFlux[TAG::NEXT *NC+ic]    = dt * compFlux[ic];
    localFlux[TAG::CURRENT *NC+ic] = -dt * compFlux[ic];

    // derivative with respect to rate
    localFluxJacobian_dRate[TAG::NEXT *NC+ic][0]    = dt * dCompFlux_dRate[ic];
    localFluxJacobian_dRate[TAG::CURRENT *NC+ic][0] = -dt * dCompFlux_dRate[ic];

    // derivative with respect to upstream pressure
    localFluxJacobian_dPresCompUp[TAG::NEXT *NC+ic][0]    = dt * dCompFlux_dPresUp[ic];
    localFluxJacobian_dPresCompUp[TAG::CURRENT *NC+ic][0] = -dt * dCompFlux_dPresUp[ic];

    // derivatives with respect to upstream component densities
    for( localIndex jdof = 0; jdof < NC; ++jdof )
    {
      localFluxJacobian_dPresCompUp[TAG::NEXT *NC+ic][jdof+1]    =  dt * dCompFlux_dCompDensUp[ic][jdof];
      localFluxJacobian_dPresCompUp[TAG::CURRENT *NC+ic][jdof+1] = -dt * dCompFlux_dCompDensUp[ic][jdof];
    }
  }
}

template< localIndex NC >
void
FluxKernel::
  launch( localIndex const size,
          globalIndex const rankOffset,
          WellControls const & wellControls,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & connRate,
          arrayView1d< real64 const > const & dConnRate,
          arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompFrac,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dWellElemCompFrac_dCompDens,
          real64 const & dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{

  using namespace CompositionalMultiphaseUtilities;

  WellControls::Type const wellType = wellControls.getType();
  arrayView1d< real64 const > const & injection = wellControls.getInjectionStream();

  // loop over the well elements to compute the fluxes between elements
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
  {

    // create local work arrays
    real64 compFracUp[NC]{};
    real64 dCompFrac_dCompDensUp[NC][NC]{};

    real64 compFlux[NC]{};
    real64 dCompFlux_dRate[NC]{};
    real64 dCompFlux_dPresUp[NC]{};
    real64 dCompFlux_dCompDensUp[NC][NC]{};

    // Step 1) decide the upwind well element

    /*  currentConnRate < 0 flow from iwelem to iwelemNext
     *  currentConnRate > 0 flow from iwelemNext to iwelem
     *  With this convention, currentConnRate < 0 at the last connection for a producer
     *                        currentConnRate > 0 at the last connection for a injector
     */

    localIndex const iwelemNext = nextWellElemIndex[iwelem];
    real64 const currentConnRate = connRate[iwelem] + dConnRate[iwelem];
    localIndex iwelemUp = -1;

    if( iwelemNext < 0 && wellType == WellControls::Type::INJECTOR ) // exit connection, injector
    {
      // we still need to define iwelemUp for Jacobian assembly
      iwelemUp = iwelem;

      // just copy the injection stream into compFrac
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        compFracUp[ic] = injection[ic];
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dCompFrac_dCompDensUp[ic][jc] = 0.0;
        }
      }
    }
    else
    {
      // first set iwelemUp to the upstream cell
      if( ( iwelemNext < 0 && wellType == WellControls::Type::PRODUCER )  // exit connection, producer
          || currentConnRate < 0 ) // not an exit connection, iwelem is upstream
      {
        iwelemUp = iwelem;
      }
      else // not an exit connection, iwelemNext is upstream
      {
        iwelemUp = iwelemNext;
      }

      // copy the vars of iwelemUp into compFrac
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        compFracUp[ic] = wellElemCompFrac[iwelemUp][ic];
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dCompFrac_dCompDensUp[ic][jc] = dWellElemCompFrac_dCompDens[iwelemUp][ic][jc];
        }
      }
    }

    // Step 2) compute upstream transport coefficient

    for( localIndex ic = 0; ic < NC; ++ic )
    {
      compFlux[ic]          = compFracUp[ic] * currentConnRate;
      dCompFlux_dRate[ic]   = compFracUp[ic];
      dCompFlux_dPresUp[ic] = 0.0; // none of these quantities depend on pressure
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dCompFlux_dCompDensUp[ic][jc] = dCompFrac_dCompDensUp[ic][jc] * currentConnRate;
      }
    }

    globalIndex const offsetUp = wellElemDofNumber[iwelemUp];
    globalIndex const offsetCurrent = wellElemDofNumber[iwelem];

    if( iwelemNext < 0 )  // exit connection
    {
      // for this case, we only need NC mass conservation equations
      // so we do not use the arrays initialized before the loop
      real64 oneSidedFlux[NC]{};
      real64 oneSidedFluxJacobian_dRate[NC][1]{};
      real64 oneSidedFluxJacobian_dPresCompUp[NC][NC+1]{};

      computeExit< NC >( dt,
                         compFlux,
                         dCompFlux_dRate,
                         dCompFlux_dPresUp,
                         dCompFlux_dCompDensUp,
                         oneSidedFlux,
                         oneSidedFluxJacobian_dRate,
                         oneSidedFluxJacobian_dPresCompUp );


      globalIndex oneSidedEqnRowIndices[NC]{};
      globalIndex oneSidedDofColIndices_dPresCompUp[NC+1]{};
      globalIndex oneSidedDofColIndices_dRate = 0;

      // jacobian indices
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        // mass balance equations for all components
        oneSidedEqnRowIndices[ic] = offsetUp + ROFFSET::MASSBAL + ic - rankOffset;
      }

      // in the dof ordering used in this class, there are 1 pressure dofs
      // and NC compDens dofs before the rate dof in this block
      localIndex const dRateColOffset = COFFSET::DCOMP + NC;
      oneSidedDofColIndices_dRate = offsetCurrent + dRateColOffset;

      for( localIndex jdof = 0; jdof < NC+1; ++jdof )
      {
        // dofs are the **upstream** pressure and component densities
        oneSidedDofColIndices_dPresCompUp[jdof] = offsetUp + COFFSET::DPRES + jdof;
      }

      // Apply equation/variable change transformation(s)
      real64 work[NC+1];
      shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, 1, oneSidedFluxJacobian_dRate, work );
      shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NC + 1, oneSidedFluxJacobian_dPresCompUp, work );
      shiftElementsAheadByOneAndReplaceFirstElementWithSum( NC, oneSidedFlux );

      for( localIndex i = 0; i < NC; ++i )
      {
        if( oneSidedEqnRowIndices[i] >= 0 && oneSidedEqnRowIndices[i] < localMatrix.numRows() )
        {
          localMatrix.addToRow< parallelDeviceAtomic >( oneSidedEqnRowIndices[i],
                                                        &oneSidedDofColIndices_dRate,
                                                        oneSidedFluxJacobian_dRate[i],
                                                        1 );
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( oneSidedEqnRowIndices[i],
                                                                            oneSidedDofColIndices_dPresCompUp,
                                                                            oneSidedFluxJacobian_dPresCompUp[i],
                                                                            NC+1 );
          atomicAdd( parallelDeviceAtomic{}, &localRhs[oneSidedEqnRowIndices[i]], oneSidedFlux[i] );
        }
      }
    }
    else // not an exit connection
    {
      real64 localFlux[2*NC]{};
      real64 localFluxJacobian_dRate[2*NC][1]{};
      real64 localFluxJacobian_dPresCompUp[2*NC][NC+1]{};

      compute< NC >( dt,
                     compFlux,
                     dCompFlux_dRate,
                     dCompFlux_dPresUp,
                     dCompFlux_dCompDensUp,
                     localFlux,
                     localFluxJacobian_dRate,
                     localFluxJacobian_dPresCompUp );


      globalIndex eqnRowIndices[2*NC]{};
      globalIndex dofColIndices_dPresCompUp[NC+1]{};
      globalIndex dofColIndices_dRate = 0;

      globalIndex const offsetNext = wellElemDofNumber[iwelemNext];

      // jacobian indices
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        // mass balance equations for all components
        eqnRowIndices[TAG::NEXT *NC+ic]    = offsetNext + ROFFSET::MASSBAL + ic - rankOffset;
        eqnRowIndices[TAG::CURRENT *NC+ic] = offsetCurrent + ROFFSET::MASSBAL + ic - rankOffset;
      }

      // in the dof ordering used in this class, there are 1 pressure dofs
      // and NC compDens dofs before the rate dof in this block
      localIndex const dRateColOffset = COFFSET::DCOMP + NC;
      dofColIndices_dRate = offsetCurrent + dRateColOffset;

      for( localIndex jdof = 0; jdof < NC+1; ++jdof )
      {
        // dofs are the **upstream** pressure and component densities
        dofColIndices_dPresCompUp[jdof] = offsetUp + COFFSET::DPRES + jdof;
      }

      // Apply equation/variable change transformation(s)
      real64 work[NC+1];
      shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, 1, 2, localFluxJacobian_dRate, work );
      shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NC + 1, 2, localFluxJacobian_dPresCompUp, work );
      shiftBlockElementsAheadByOneAndReplaceFirstElementWithSum( NC, 2, localFlux );

      for( localIndex i = 0; i < 2*NC; ++i )
      {
        if( eqnRowIndices[i] >= 0 && eqnRowIndices[i] < localMatrix.numRows() )
        {
          localMatrix.addToRow< parallelDeviceAtomic >( eqnRowIndices[i],
                                                        &dofColIndices_dRate,
                                                        localFluxJacobian_dRate[i],
                                                        1 );
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndices[i],
                                                                            dofColIndices_dPresCompUp,
                                                                            localFluxJacobian_dPresCompUp[i],
                                                                            NC+1 );
          atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndices[i]], localFlux[i] );
        }
      }
    }
  } );
}

#define INST_FluxKernel( NC ) \
  template \
  void FluxKernel:: \
    launch< NC >( localIndex const size, \
                  globalIndex const rankOffset, \
                  WellControls const & wellControls, \
                  arrayView1d< globalIndex const > const & wellElemDofNumber, \
                  arrayView1d< localIndex const > const & nextWellElemIndex, \
                  arrayView1d< real64 const > const & connRate, \
                  arrayView1d< real64 const > const & dConnRate, \
                  arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompFrac, \
                  arrayView3d< real64 const, compflow::USD_COMP_DC > const & dWellElemCompFrac_dCompDens, \
                  real64 const & dt, \
                  CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                  arrayView1d< real64 > const & localRhs )

INST_FluxKernel( 1 );
INST_FluxKernel( 2 );
INST_FluxKernel( 3 );
INST_FluxKernel( 4 );
INST_FluxKernel( 5 );

/******************************** PressureRelationKernel ********************************/

template< localIndex NC >
GEOSX_HOST_DEVICE
void
PressureRelationKernel::
  compute( real64 const & gravCoef,
           real64 const & gravCoefNext,
           real64 const & pres,
           real64 const & presNext,
           real64 const & dPres,
           real64 const & dPresNext,
           real64 const & totalMassDens,
           real64 const & totalMassDensNext,
           real64 const & dTotalMassDens_dPres,
           real64 const & dTotalMassDens_dPresNext,
           arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dTotalMassDens_dCompDens,
           arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dTotalMassDens_dCompDensNext,
           real64 & localPresRel,
           real64 ( & localPresRelJacobian )[2*(NC+1)] )
{
  // local working variables and arrays
  real64 dAvgMassDens_dCompCurrent[NC]{};
  real64 dAvgMassDens_dCompNext[NC]{};

  // compute the average density at the interface between well elements
  real64 const avgMassDens = 0.5 * ( totalMassDensNext + totalMassDens );
  real64 const dAvgMassDens_dPresNext    = 0.5 * dTotalMassDens_dPresNext;
  real64 const dAvgMassDens_dPresCurrent = 0.5 * dTotalMassDens_dPres;
  for( localIndex ic = 0; ic < NC; ++ic )
  {
    dAvgMassDens_dCompNext[ic]    = 0.5 * dTotalMassDens_dCompDensNext[ic];
    dAvgMassDens_dCompCurrent[ic] = 0.5 * dTotalMassDens_dCompDens[ic];
  }

  // compute depth diff times acceleration
  real64 const gravD = gravCoefNext - gravCoef;

  // TODO: add friction and acceleration terms

  localPresRel = ( presNext + dPresNext  - pres - dPres - avgMassDens * gravD );
  localPresRelJacobian[TAG::NEXT *(NC+1)]    = ( 1 - dAvgMassDens_dPresNext * gravD );
  localPresRelJacobian[TAG::CURRENT *(NC+1)] = ( -1 - dAvgMassDens_dPresCurrent * gravD );

  for( localIndex ic = 0; ic < NC; ++ic )
  {
    localPresRelJacobian[TAG::NEXT *(NC+1) + ic+1]    = -dAvgMassDens_dCompNext[ic] * gravD;
    localPresRelJacobian[TAG::CURRENT *(NC+1) + ic+1] = -dAvgMassDens_dCompCurrent[ic] * gravD;
  }
}

template< localIndex NC >
void
PressureRelationKernel::
  launch( localIndex const size,
          globalIndex const rankOffset,
          bool const isLocallyOwned,
          localIndex const iwelemControl,
          localIndex const targetPhaseIndex,
          WellControls const & wellControls,
          real64 const & timeAtEndOfStep,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & wellElemPressure,
          arrayView1d< real64 const > const & dWellElemPressure,
          arrayView1d< real64 const > const & wellElemTotalMassDens,
          arrayView1d< real64 const > const & dWellElemTotalMassDens_dPres,
          arrayView2d< real64 const, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens_dCompDens,
          bool & controlHasSwitched,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{

  // static well control data
  WellControls::Type const wellType = wellControls.getType();
  WellControls::Control const currentControl = wellControls.getControl();
  real64 const targetBHP = wellControls.getTargetBHP( timeAtEndOfStep );
  real64 const targetTotalRate = wellControls.getTargetTotalRate( timeAtEndOfStep );
  real64 const targetPhaseRate = wellControls.getTargetPhaseRate( timeAtEndOfStep );

  // dynamic well control data
  real64 const & currentBHP =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() );
  real64 const & dCurrentBHP_dPres =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentBHP_dPresString() );
  arrayView1d< real64 const > const & dCurrentBHP_dCompDens =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentBHP_dCompDensString() );

  arrayView1d< real64 const > const & currentPhaseVolRate =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() );
  arrayView1d< real64 const > const & dCurrentPhaseVolRate_dPres =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentPhaseVolRate_dPresString() );
  arrayView2d< real64 const > const & dCurrentPhaseVolRate_dCompDens =
    wellControls.getReference< array2d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentPhaseVolRate_dCompDensString() );
  arrayView1d< real64 const > const & dCurrentPhaseVolRate_dRate =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentPhaseVolRate_dRateString() );

  real64 const & currentTotalVolRate =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() );
  real64 const & dCurrentTotalVolRate_dPres =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentTotalVolRate_dPresString() );
  arrayView1d< real64 const > const & dCurrentTotalVolRate_dCompDens =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentTotalVolRate_dCompDensString() );
  real64 const & dCurrentTotalVolRate_dRate =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentTotalVolRate_dRateString() );

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
                                            targetPhaseIndex,
                                            targetBHP,
                                            targetPhaseRate,
                                            targetTotalRate,
                                            currentBHP,
                                            currentPhaseVolRate,
                                            currentTotalVolRate,
                                            newControl );
      if( currentControl != newControl )
      {
        switchControl.max( 1 );
      }

      ControlEquationHelper::compute< NC >( rankOffset,
                                            newControl,
                                            targetPhaseIndex,
                                            targetBHP,
                                            targetPhaseRate,
                                            targetTotalRate,
                                            currentBHP,
                                            dCurrentBHP_dPres,
                                            dCurrentBHP_dCompDens,
                                            currentPhaseVolRate,
                                            dCurrentPhaseVolRate_dPres,
                                            dCurrentPhaseVolRate_dCompDens,
                                            dCurrentPhaseVolRate_dRate,
                                            currentTotalVolRate,
                                            dCurrentTotalVolRate_dPres,
                                            dCurrentTotalVolRate_dCompDens,
                                            dCurrentTotalVolRate_dRate,
                                            wellElemDofNumber[iwelemControl],
                                            localMatrix,
                                            localRhs );

      // TODO: for consistency, we should assemble here, not in compute...

    }
    else if( iwelemNext >= 0 ) // if iwelemNext >= 0, form momentum equation
    {

      real64 localPresRel = 0;
      real64 localPresRelJacobian[2*(NC+1)]{};

      compute< NC >( wellElemGravCoef[iwelem],
                     wellElemGravCoef[iwelemNext],
                     wellElemPressure[iwelem],
                     wellElemPressure[iwelemNext],
                     dWellElemPressure[iwelem],
                     dWellElemPressure[iwelemNext],
                     wellElemTotalMassDens[iwelem],
                     wellElemTotalMassDens[iwelemNext],
                     dWellElemTotalMassDens_dPres[iwelem],
                     dWellElemTotalMassDens_dPres[iwelemNext],
                     dWellElemTotalMassDens_dCompDens[iwelem],
                     dWellElemTotalMassDens_dCompDens[iwelemNext],
                     localPresRel,
                     localPresRelJacobian );


      // local working variables and arrays
      globalIndex dofColIndices[2*(NC+1)];

      globalIndex const eqnRowIndex = wellElemDofNumber[iwelem] + ROFFSET::CONTROL - rankOffset;
      dofColIndices[TAG::NEXT *(NC+1)]    = wellElemDofNumber[iwelemNext] + COFFSET::DPRES;
      dofColIndices[TAG::CURRENT *(NC+1)] = wellElemDofNumber[iwelem] + COFFSET::DPRES;

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dofColIndices[TAG::NEXT *(NC+1) + ic+1]    = wellElemDofNumber[iwelemNext] + COFFSET::DCOMP + ic;
        dofColIndices[TAG::CURRENT *(NC+1) + ic+1] = wellElemDofNumber[iwelem] + COFFSET::DCOMP + ic;
      }

      if( eqnRowIndex >= 0 && eqnRowIndex < localMatrix.numRows() )
      {
        localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndex,
                                                                          dofColIndices,
                                                                          localPresRelJacobian,
                                                                          2 * (NC+1) );
        atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndex], localPresRel );
      }
    }
  } );
  controlHasSwitched = ( switchControl.get() == 1 );
}

#define INST_PressureRelationKernel( NC ) \
  template \
  void PressureRelationKernel:: \
    launch< NC >( localIndex const size, \
                  globalIndex const rankOffset, \
                  bool const isLocallyOwned, \
                  localIndex const iwelemControl, \
                  localIndex const targetPhaseIndex, \
                  WellControls const & wellControls, \
                  real64 const & timeAtEndOfStep, \
                  arrayView1d< globalIndex const > const & wellElemDofNumber, \
                  arrayView1d< real64 const > const & wellElemGravCoef, \
                  arrayView1d< localIndex const > const & nextWellElemIndex, \
                  arrayView1d< real64 const > const & wellElemPressure, \
                  arrayView1d< real64 const > const & dWellElemPressure, \
                  arrayView1d< real64 const > const & wellElemTotalMassDens, \
                  arrayView1d< real64 const > const & dWellElemTotalMassDens_dPres, \
                  arrayView2d< real64 const, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens_dCompDens, \
                  bool & controlHasSwitched, \
                  CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                  arrayView1d< real64 > const & localRhs )

INST_PressureRelationKernel( 1 );
INST_PressureRelationKernel( 2 );
INST_PressureRelationKernel( 3 );
INST_PressureRelationKernel( 4 );
INST_PressureRelationKernel( 5 );


/******************************** PerforationKernel ********************************/

template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
PerforationKernel::
  compute( real64 const & resPres,
           real64 const & dResPres,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & resPhaseVolFrac,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & dResPhaseVolFrac_dPres,
           arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dResPhaseVolFrac_dComp,
           arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dResCompFrac_dCompDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & resPhaseDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & dResPhaseDens_dPres,
           arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dResPhaseDens_dComp,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & resPhaseVisc,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & dResPhaseVisc_dPres,
           arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dResPhaseVisc_dComp,
           arraySlice2d< real64 const, multifluid::USD_PHASE_COMP - 2 > const & resPhaseCompFrac,
           arraySlice2d< real64 const, multifluid::USD_PHASE_COMP - 2 > const & dResPhaseCompFrac_dPres,
           arraySlice3d< real64 const, multifluid::USD_PHASE_COMP_DC - 2 > const & dResPhaseCompFrac_dComp,
           arraySlice1d< real64 const, relperm::USD_RELPERM - 2 > const & resPhaseRelPerm,
           arraySlice2d< real64 const, relperm::USD_RELPERM_DS - 2 > const & dResPhaseRelPerm_dPhaseVolFrac,
           real64 const & wellElemGravCoef,
           real64 const & wellElemPres,
           real64 const & dWellElemPres,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & wellElemCompDens,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & dWellElemCompDens,
           real64 const & wellElemTotalMassDens,
           real64 const & dWellElemTotalMassDens_dPres,
           arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dWellElemTotalMassDens_dCompDens,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & wellElemCompFrac,
           arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dWellElemCompFrac_dCompDens,
           real64 const & perfGravCoef,
           real64 const & trans,
           arraySlice1d< real64 > const & compPerfRate,
           arraySlice2d< real64 > const & dCompPerfRate_dPres,
           arraySlice3d< real64 > const & dCompPerfRate_dComp )
{

  // local working variables and arrays
  real64 pres[2]{};
  real64 dPres_dP[2]{};
  real64 dPres_dC[2][NC]{};
  real64 dFlux_dP[2]{};
  real64 dFlux_dC[2][NC]{};
  real64 dMult_dP[2]{};
  real64 dMult_dC[2][NC]{};
  real64 dPotDiff_dP[2]{};
  real64 dPotDiff_dC[2][NC]{};
  real64 multiplier[2]{};

  real64 dResTotalMob_dC[NC]{};
  real64 dDens_dC[NC]{};
  real64 dVisc_dC[NC]{};
  real64 dRelPerm_dC[NC]{};
  real64 dMob_dC[NC]{};
  real64 dCompFrac_dCompDens[NC]{};


  // Step 1: reset the perforation rates

  for( localIndex ic = 0; ic < NC; ++ic )
  {
    compPerfRate[ic] = 0.0;
    for( localIndex ke = 0; ke < 2; ++ke )
    {
      dCompPerfRate_dPres[ke][ic] = 0.0;
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dCompPerfRate_dComp[ke][ic][jc] = 0.0;
      }
    }
  }


  // Step 2: copy the variables from the reservoir and well element

  // a) get reservoir variables

  pres[TAG::RES] = resPres + dResPres;
  dPres_dP[TAG::RES] = 1.0;
  multiplier[TAG::RES] = 1.0;

  // Here in the absence of a buoyancy term we assume that the reservoir cell is perforated at its center
  // TODO: add a buoyancy term for the reservoir side here


  // b) get well variables

  pres[TAG::WELL] = wellElemPres + dWellElemPres;
  dPres_dP[TAG::WELL] = 1.0;
  multiplier[TAG::WELL] = -1.0;

  real64 const gravD = ( perfGravCoef - wellElemGravCoef );

  pres[TAG::WELL] += wellElemTotalMassDens * gravD;
  dPres_dP[TAG::WELL] += dWellElemTotalMassDens_dPres * gravD;
  for( localIndex ic = 0; ic < NC; ++ic )
  {
    dPres_dC[TAG::WELL][ic] += dWellElemTotalMassDens_dCompDens[ic] * gravD;
  }


  // Step 3: compute potential difference

  real64 potDiff = 0.0;
  for( localIndex i = 0; i < 2; ++i )
  {
    potDiff += multiplier[i] * trans * pres[i];
    dPotDiff_dP[i] += multiplier[i] * trans * dPres_dP[i];

    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dPotDiff_dC[i][ic] += multiplier[i] * trans * dPres_dC[i][ic];
    }
  }


  // Step 4: upwinding based on the flow direction

  real64 flux = 0.0;
  if( potDiff >= 0 )  // ** reservoir cell is upstream **
  {

    // loop over phases, compute and upwind phase flux
    // and sum contributions to each component's perforation rate
    for( localIndex ip = 0; ip < NP; ++ip )
    {

      // skip the rest of the calculation if the phase is absent
      bool const phaseExists = (resPhaseVolFrac[ip] > 0);
      if( !phaseExists )
      {
        continue;
      }

      // here, we have to recompute the reservoir phase mobility (not including density)

      // density
      real64 const resDens = resPhaseDens[ip];
      real64 const dResDens_dP  = dResPhaseDens_dPres[ip];
      applyChainRule( NC, dResCompFrac_dCompDens,
                      dResPhaseDens_dComp[ip],
                      dDens_dC );

      // viscosity
      real64 const resVisc = resPhaseVisc[ip];
      real64 const dResVisc_dP  = dResPhaseVisc_dPres[ip];
      applyChainRule( NC, dResCompFrac_dCompDens,
                      dResPhaseVisc_dComp[ip],
                      dVisc_dC );

      // relative permeability
      real64 const resRelPerm = resPhaseRelPerm[ip];
      real64 dResRelPerm_dP = 0.0;
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dRelPerm_dC[jc] = 0;
      }

      for( localIndex jp = 0; jp < NP; ++jp )
      {
        real64 const dResRelPerm_dS = dResPhaseRelPerm_dPhaseVolFrac[ip][jp];
        dResRelPerm_dP += dResRelPerm_dS * dResPhaseVolFrac_dPres[jp];

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dRelPerm_dC[jc] += dResRelPerm_dS * dResPhaseVolFrac_dComp[jp][jc];
        }
      }

      // compute the reservoir phase mobility, including phase density
      real64 const resPhaseMob = resDens * resRelPerm / resVisc;
      real64 const dResPhaseMob_dPres = dResRelPerm_dP * resDens / resVisc
                                        + resPhaseMob * (dResDens_dP / resDens - dResVisc_dP / resVisc);
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dMob_dC[jc] = dRelPerm_dC[jc] * resDens / resVisc
                      + resPhaseMob * (dDens_dC[jc] / resDens - dVisc_dC[jc] / resVisc);
      }

      // compute the phase flux and derivatives using upstream cell mobility
      flux = resPhaseMob * potDiff;
      dFlux_dP[TAG::RES]  = dResPhaseMob_dPres * potDiff + resPhaseMob * dPotDiff_dP[TAG::RES];
      dFlux_dP[TAG::WELL] = resPhaseMob *  dPotDiff_dP[TAG::WELL];

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dFlux_dC[TAG::RES][ic] = dMob_dC[ic] * potDiff + resPhaseMob * dPotDiff_dC[TAG::RES][ic];
        dFlux_dC[TAG::WELL][ic] = resPhaseMob * dPotDiff_dC[TAG::WELL][ic];
      }

      // increment component fluxes
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        compPerfRate[ic] += flux * resPhaseCompFrac[ip][ic];

        dCompPerfRate_dPres[TAG::RES][ic]  += resPhaseCompFrac[ip][ic] * dFlux_dP[TAG::RES];
        dCompPerfRate_dPres[TAG::RES][ic]  += dResPhaseCompFrac_dPres[ip][ic] * flux;
        dCompPerfRate_dPres[TAG::WELL][ic] += resPhaseCompFrac[ip][ic] * dFlux_dP[TAG::WELL];

        applyChainRule( NC,
                        dResCompFrac_dCompDens,
                        dResPhaseCompFrac_dComp[ip][ic],
                        dCompFrac_dCompDens );

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dCompPerfRate_dComp[TAG::RES][ic][jc]  += dFlux_dC[TAG::RES][jc] * resPhaseCompFrac[ip][ic];
          dCompPerfRate_dComp[TAG::RES][ic][jc]  += flux * dCompFrac_dCompDens[jc];
          dCompPerfRate_dComp[TAG::WELL][ic][jc] += dFlux_dC[TAG::WELL][jc] * resPhaseCompFrac[ip][ic];
        }
      }
    }
  }
  else // ** well is upstream **
  {

    real64 resTotalMob     = 0.0;
    real64 dResTotalMob_dP = 0.0;

    // we re-compute here the total mass (when useMass == 1) or molar (when useMass == 0) density
    real64 wellElemTotalDens = 0;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      wellElemTotalDens += wellElemCompDens[ic] + dWellElemCompDens[ic];
    }

    // first, compute the reservoir total mobility (excluding phase density)
    for( localIndex ip = 0; ip < NP; ++ip )
    {

      // skip the rest of the calculation if the phase is absent
      bool const phaseExists = (resPhaseVolFrac[ip] > 0);
      if( !phaseExists )
      {
        continue;
      }

      // viscosity
      real64 const resVisc = resPhaseVisc[ip];
      real64 const dResVisc_dP  = dResPhaseVisc_dPres[ip];
      applyChainRule( NC, dResCompFrac_dCompDens,
                      dResPhaseVisc_dComp[ip],
                      dVisc_dC );

      // relative permeability
      real64 const resRelPerm = resPhaseRelPerm[ip];
      real64 dResRelPerm_dP = 0.0;
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dRelPerm_dC[jc] = 0;
      }

      for( localIndex jp = 0; jp < NP; ++jp )
      {
        real64 const dResRelPerm_dS = dResPhaseRelPerm_dPhaseVolFrac[ip][jp];
        dResRelPerm_dP += dResRelPerm_dS * dResPhaseVolFrac_dPres[jp];

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dRelPerm_dC[jc] += dResRelPerm_dS * dResPhaseVolFrac_dComp[jp][jc];
        }
      }

      // increment total mobility
      resTotalMob     += resRelPerm / resVisc;
      dResTotalMob_dP += ( dResRelPerm_dP * resVisc - resRelPerm * dResVisc_dP )
                         / ( resVisc * resVisc );
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dResTotalMob_dC[ic] += ( dRelPerm_dC[ic] * resVisc - resRelPerm * dVisc_dC[ic] )
                               / ( resVisc * resVisc );
      }
    }

    // compute a potdiff multiplier = wellElemTotalDens * resTotalMob
    // wellElemTotalDens is a mass density if useMass == 1 and a molar density otherwise
    real64 const mult   = wellElemTotalDens * resTotalMob;
    dMult_dP[TAG::RES]  = wellElemTotalDens * dResTotalMob_dP;
    dMult_dP[TAG::WELL] = 0.0; // because totalDens does not depend on pressure

    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dMult_dC[TAG::RES][ic]  = wellElemTotalDens * dResTotalMob_dC[ic];
      dMult_dC[TAG::WELL][ic] = resTotalMob;
    }

    // compute the volumetric flux and derivatives using upstream cell mobility
    flux = mult * potDiff;
    dFlux_dP[TAG::RES]  = dMult_dP[TAG::RES] * potDiff + mult * dPotDiff_dP[TAG::RES];
    dFlux_dP[TAG::WELL] = dMult_dP[TAG::WELL] * potDiff + mult * dPotDiff_dP[TAG::WELL];

    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dFlux_dC[TAG::RES][ic]  = dMult_dC[TAG::RES][ic] * potDiff + mult * dPotDiff_dC[TAG::RES][ic];
      dFlux_dC[TAG::WELL][ic] = dMult_dC[TAG::WELL][ic] * potDiff + mult * dPotDiff_dC[TAG::WELL][ic];
    }

    // compute component fluxes
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      compPerfRate[ic] += wellElemCompFrac[ic] * flux;
      dCompPerfRate_dPres[TAG::RES][ic]  = wellElemCompFrac[ic] * dFlux_dP[TAG::RES];
      dCompPerfRate_dPres[TAG::WELL][ic] = wellElemCompFrac[ic] * dFlux_dP[TAG::WELL];

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dCompPerfRate_dComp[TAG::RES][ic][jc]  += wellElemCompFrac[ic] * dFlux_dC[TAG::RES][jc];
        dCompPerfRate_dComp[TAG::WELL][ic][jc] += wellElemCompFrac[ic] * dFlux_dC[TAG::WELL][jc];
        dCompPerfRate_dComp[TAG::WELL][ic][jc] += dWellElemCompFrac_dCompDens[ic][jc] * flux;
      }
    }
  }
}

template< localIndex NC, localIndex NP >
void
PerforationKernel::
  launch( localIndex const size,
          ElementViewConst< arrayView1d< real64 const > > const & resPres,
          ElementViewConst< arrayView1d< real64 const > > const & dResPres,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & resPhaseVolFrac,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dResPhaseVolFrac_dPres,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dResPhaseVolFrac_dComp,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dResCompFrac_dCompDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & resPhaseDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dResPhaseDens_dPres,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dResPhaseDens_dComp,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & resPhaseVisc,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dResPhaseVisc_dPres,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dResPhaseVisc_dComp,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & resPhaseCompFrac,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dResPhaseCompFrac_dPres,
          ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dResPhaseCompFrac_dComp,
          ElementViewConst< arrayView3d< real64 const, relperm::USD_RELPERM > > const & resPhaseRelPerm,
          ElementViewConst< arrayView4d< real64 const, relperm::USD_RELPERM_DS > > const & dResPhaseRelPerm_dPhaseVolFrac,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< real64 const > const & wellElemPres,
          arrayView1d< real64 const > const & dWellElemPres,
          arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dWellElemCompDens,
          arrayView1d< real64 const > const & wellElemTotalMassDens,
          arrayView1d< real64 const > const & dWellElemTotalMassDens_dPres,
          arrayView2d< real64 const, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens_dCompDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompFrac,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dWellElemCompFrac_dCompDens,
          arrayView1d< real64 const > const & perfGravCoef,
          arrayView1d< localIndex const > const & perfWellElemIndex,
          arrayView1d< real64 const > const & perfTrans,
          arrayView1d< localIndex const > const & resElementRegion,
          arrayView1d< localIndex const > const & resElementSubRegion,
          arrayView1d< localIndex const > const & resElementIndex,
          arrayView2d< real64 > const & compPerfRate,
          arrayView3d< real64 > const & dCompPerfRate_dPres,
          arrayView4d< real64 > const & dCompPerfRate_dComp )
{

  // loop over the perforations to compute the perforation rates
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const iperf )
  {

    // get the index of the reservoir elem
    localIndex const er  = resElementRegion[iperf];
    localIndex const esr = resElementSubRegion[iperf];
    localIndex const ei  = resElementIndex[iperf];

    // get the index of the well elem
    localIndex const iwelem = perfWellElemIndex[iperf];

    compute< NC, NP >( resPres[er][esr][ei],
                       dResPres[er][esr][ei],
                       resPhaseVolFrac[er][esr][ei],
                       dResPhaseVolFrac_dPres[er][esr][ei],
                       dResPhaseVolFrac_dComp[er][esr][ei],
                       dResCompFrac_dCompDens[er][esr][ei],
                       resPhaseDens[er][esr][ei][0],
                       dResPhaseDens_dPres[er][esr][ei][0],
                       dResPhaseDens_dComp[er][esr][ei][0],
                       resPhaseVisc[er][esr][ei][0],
                       dResPhaseVisc_dPres[er][esr][ei][0],
                       dResPhaseVisc_dComp[er][esr][ei][0],
                       resPhaseCompFrac[er][esr][ei][0],
                       dResPhaseCompFrac_dPres[er][esr][ei][0],
                       dResPhaseCompFrac_dComp[er][esr][ei][0],
                       resPhaseRelPerm[er][esr][ei][0],
                       dResPhaseRelPerm_dPhaseVolFrac[er][esr][ei][0],
                       wellElemGravCoef[iwelem],
                       wellElemPres[iwelem],
                       dWellElemPres[iwelem],
                       wellElemCompDens[iwelem],
                       dWellElemCompDens[iwelem],
                       wellElemTotalMassDens[iwelem],
                       dWellElemTotalMassDens_dPres[iwelem],
                       dWellElemTotalMassDens_dCompDens[iwelem],
                       wellElemCompFrac[iwelem],
                       dWellElemCompFrac_dCompDens[iwelem],
                       perfGravCoef[iperf],
                       perfTrans[iperf],
                       compPerfRate[iperf],
                       dCompPerfRate_dPres[iperf],
                       dCompPerfRate_dComp[iperf] );

  } );
}

#define INST_PerforationKernel( NC, NP ) \
  template \
  void PerforationKernel:: \
    launch< NC, NP >( localIndex const size, \
                      ElementViewConst< arrayView1d< real64 const > > const & resPres, \
                      ElementViewConst< arrayView1d< real64 const > > const & dResPres, \
                      ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & resPhaseVolFrac, \
                      ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dResPhaseVolFrac_dPres, \
                      ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dResPhaseVolFrac_dComp, \
                      ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dResCompFrac_dCompDens, \
                      ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & resPhaseDens, \
                      ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dResPhaseDens_dPres, \
                      ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dResPhaseDens_dComp, \
                      ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & resPhaseVisc, \
                      ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dResPhaseVisc_dPres, \
                      ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dResPhaseVisc_dComp, \
                      ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & resPhaseCompFrac, \
                      ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dResPhaseCompFrac_dPres, \
                      ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dResPhaseCompFrac_dComp, \
                      ElementViewConst< arrayView3d< real64 const, relperm::USD_RELPERM > > const & resPhaseRelPerm, \
                      ElementViewConst< arrayView4d< real64 const, relperm::USD_RELPERM_DS > > const & dResPhaseRelPerm_dPhaseVolFrac, \
                      arrayView1d< real64 const > const & wellElemGravCoef, \
                      arrayView1d< real64 const > const & wellElemPres, \
                      arrayView1d< real64 const > const & dWellElemPres, \
                      arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompDens, \
                      arrayView2d< real64 const, compflow::USD_COMP > const & dWellElemCompDens, \
                      arrayView1d< real64 const > const & wellElemTotalMassDens, \
                      arrayView1d< real64 const > const & dWellElemTotalMassDens_dPres, \
                      arrayView2d< real64 const, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens_dCompDens, \
                      arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompFrac, \
                      arrayView3d< real64 const, compflow::USD_COMP_DC > const & dWellElemCompFrac_dCompDens, \
                      arrayView1d< real64 const > const & perfGravCoef, \
                      arrayView1d< localIndex const > const & perfWellElemIndex, \
                      arrayView1d< real64 const > const & perfTrans, \
                      arrayView1d< localIndex const > const & resElementRegion, \
                      arrayView1d< localIndex const > const & resElementSubRegion, \
                      arrayView1d< localIndex const > const & resElementIndex, \
                      arrayView2d< real64 > const & compPerfRate, \
                      arrayView3d< real64 > const & dCompPerfRate_dPres, \
                      arrayView4d< real64 > const & dCompPerfRate_dComp )

INST_PerforationKernel( 1, 2 );
INST_PerforationKernel( 2, 2 );
INST_PerforationKernel( 3, 2 );
INST_PerforationKernel( 4, 2 );
INST_PerforationKernel( 5, 2 );
INST_PerforationKernel( 1, 3 );
INST_PerforationKernel( 2, 3 );
INST_PerforationKernel( 3, 3 );
INST_PerforationKernel( 4, 3 );
INST_PerforationKernel( 5, 3 );

/******************************** AccumulationKernel ********************************/

template< localIndex NC >
GEOSX_HOST_DEVICE
void
AccumulationKernel::
  compute( localIndex const numPhases,
           real64 const & volume,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dPhaseVolFrac_dCompDens,
           arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dCompFrac_dCompDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & dPhaseDens_dPres,
           arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dPhaseDens_dComp,
           arraySlice2d< real64 const, multifluid::USD_PHASE_COMP - 2 > const & phaseCompFrac,
           arraySlice2d< real64 const, multifluid::USD_PHASE_COMP - 2 > const & dPhaseCompFrac_dPres,
           arraySlice3d< real64 const, multifluid::USD_PHASE_COMP_DC - 2 > const & dPhaseCompFrac_dComp,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFracOld,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseDensOld,
           arraySlice2d< real64 const, compflow::USD_PHASE_COMP - 1 > const & phaseCompFracOld,
           real64 ( & localAccum )[NC],
           real64 ( & localAccumJacobian )[NC][NC + 1] )
{
  // temporary work arrays
  real64 dPhaseAmount_dC[NC]{};
  real64 dPhaseCompFrac_dC[NC]{};

  // reset the local values
  for( localIndex i = 0; i < NC; ++i )
  {
    localAccum[i] = 0.0;
    for( localIndex j = 0; j < NC+1; ++j )
    {
      localAccumJacobian[i][j] = 0.0;
    }
  }

  // sum contributions to component accumulation from each phase
  for( localIndex ip = 0; ip < numPhases; ++ip )
  {
    real64 const phaseAmountNew = volume * phaseVolFrac[ip] * phaseDens[ip];
    real64 const phaseAmountOld = volume * phaseVolFracOld[ip] * phaseDensOld[ip];

    real64 const dPhaseAmount_dP = volume * (dPhaseVolFrac_dPres[ip] * phaseDens[ip]
                                             + phaseVolFrac[ip] * dPhaseDens_dPres[ip]);

    // assemble density dependence
    applyChainRule( NC, dCompFrac_dCompDens, dPhaseDens_dComp[ip], dPhaseAmount_dC );
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc] * phaseVolFrac[ip]
                            + phaseDens[ip] * dPhaseVolFrac_dCompDens[ip][jc];
      dPhaseAmount_dC[jc] *= volume;
    }

    // ic - index of component whose conservation equation is assembled
    // (i.e. row number in local matrix)
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      real64 const phaseCompAmountNew = phaseAmountNew * phaseCompFrac[ip][ic];
      real64 const phaseCompAmountOld = phaseAmountOld * phaseCompFracOld[ip][ic];

      real64 const dPhaseCompAmount_dP = dPhaseAmount_dP * phaseCompFrac[ip][ic]
                                         + phaseAmountNew * dPhaseCompFrac_dPres[ip][ic];

      localAccum[ic] += phaseCompAmountNew - phaseCompAmountOld;
      localAccumJacobian[ic][0] += dPhaseCompAmount_dP;

      // jc - index of component w.r.t. whose compositional var the derivative is being taken
      // (i.e. col number in local matrix)

      // assemble phase composition dependence
      applyChainRule( NC, dCompFrac_dCompDens, dPhaseCompFrac_dComp[ip][ic], dPhaseCompFrac_dC );
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        real64 const dPhaseCompAmount_dC = dPhaseCompFrac_dC[jc] * phaseAmountNew
                                           + phaseCompFrac[ip][ic] * dPhaseAmount_dC[jc];
        localAccumJacobian[ic][jc + 1] += dPhaseCompAmount_dC;
      }
    }
  }
}

template< localIndex NC >
void
AccumulationKernel::
  launch( localIndex const size,
          localIndex const numPhases,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< integer const > const & wellElemGhostRank,
          arrayView1d< real64 const > const & wellElemVolume,
          arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseVolFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & dWellElemPhaseVolFrac_dPres,
          arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dWellElemPhaseVolFrac_dCompDens,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dWellElemCompFrac_dCompDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & wellElemPhaseDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dWellElemPhaseDens_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dWellElemPhaseDens_dComp,
          arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & wellElemPhaseCompFrac,
          arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & dWellElemPhaseCompFrac_dPres,
          arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > const & dWellElemPhaseCompFrac_dComp,
          arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseVolFracOld,
          arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseDensOld,
          arrayView3d< real64 const, compflow::USD_PHASE_COMP > const & wellElemPhaseCompFracOld,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{

  using namespace CompositionalMultiphaseUtilities;

  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
  {

    if( wellElemGhostRank[iwelem] >= 0 )
    {
      return;
    }

    real64 localAccum[NC]{};
    real64 localAccumJacobian[NC][NC+1]{};

    compute< NC >( numPhases,
                   wellElemVolume[iwelem],
                   wellElemPhaseVolFrac[iwelem],
                   dWellElemPhaseVolFrac_dPres[iwelem],
                   dWellElemPhaseVolFrac_dCompDens[iwelem],
                   dWellElemCompFrac_dCompDens[iwelem],
                   wellElemPhaseDens[iwelem][0],
                   dWellElemPhaseDens_dPres[iwelem][0],
                   dWellElemPhaseDens_dComp[iwelem][0],
                   wellElemPhaseCompFrac[iwelem][0],
                   dWellElemPhaseCompFrac_dPres[iwelem][0],
                   dWellElemPhaseCompFrac_dComp[iwelem][0],
                   wellElemPhaseVolFracOld[iwelem],
                   wellElemPhaseDensOld[iwelem],
                   wellElemPhaseCompFracOld[iwelem],
                   localAccum,
                   localAccumJacobian );

    // set the equation row indices to be the mass balance equations for all components
    localIndex eqnRowIndices[NC]{};
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      eqnRowIndices[ic] = wellElemDofNumber[iwelem] + ROFFSET::MASSBAL + ic - rankOffset;
    }

    // set DOF col indices for this block
    globalIndex dofColIndices[NC+1]{};
    for( localIndex idof = 0; idof < NC+1; ++idof )
    {
      dofColIndices[idof] = wellElemDofNumber[iwelem] + COFFSET::DPRES + idof;
    }

    // Apply equation/variable change transformation(s)
    real64 work[NC+1];
    shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NC + 1, localAccumJacobian, work );
    shiftElementsAheadByOneAndReplaceFirstElementWithSum( NC, localAccum );

    // add contribution to residual and jacobian
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      localRhs[eqnRowIndices[ic]] += localAccum[ic];
      localMatrix.addToRow< serialAtomic >( eqnRowIndices[ic],
                                            dofColIndices,
                                            localAccumJacobian[ic],
                                            NC+1 );
    }
  } );
}

#define INST_AccumulationKernel( NC ) \
  template \
  void AccumulationKernel:: \
    launch< NC >( localIndex const size, \
                  localIndex const numPhases, \
                  globalIndex const rankOffset, \
                  arrayView1d< globalIndex const > const & wellElemDofNumber, \
                  arrayView1d< integer const > const & wellElemGhostRank, \
                  arrayView1d< real64 const > const & wellElemVolume, \
                  arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseVolFrac, \
                  arrayView2d< real64 const, compflow::USD_PHASE > const & dWellElemPhaseVolFrac_dPres, \
                  arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dWellElemPhaseVolFrac_dCompDens, \
                  arrayView3d< real64 const, compflow::USD_COMP_DC > const & dWellElemCompFrac_dCompDens, \
                  arrayView3d< real64 const, multifluid::USD_PHASE > const & wellElemPhaseDens, \
                  arrayView3d< real64 const, multifluid::USD_PHASE > const & dWellElemPhaseDens_dPres, \
                  arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dWellElemPhaseDens_dComp, \
                  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & wellElemPhaseCompFrac, \
                  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & dWellElemPhaseCompFrac_dPres, \
                  arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > const & dWellElemPhaseCompFrac_dComp, \
                  arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseVolFracOld, \
                  arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseDensOld, \
                  arrayView3d< real64 const, compflow::USD_PHASE_COMP > const & wellElemPhaseCompFracOld, \
                  CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                  arrayView1d< real64 > const & localRhs )

INST_AccumulationKernel( 1 );
INST_AccumulationKernel( 2 );
INST_AccumulationKernel( 3 );
INST_AccumulationKernel( 4 );
INST_AccumulationKernel( 5 );

/******************************** VolumeBalanceKernel ********************************/

template< localIndex NC >
GEOSX_HOST_DEVICE
void
VolumeBalanceKernel::
  compute( localIndex const numPhases,
           real64 const & volume,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dPhaseVolFrac_dComp,
           real64 & localVolBalance,
           real64 ( & localVolBalanceJacobian )[NC+1] )
{
  localVolBalance = 1.0;
  for( localIndex ic = 0; ic < NC+1; ++ic )
  {
    localVolBalanceJacobian[ic] = 0.0;
  }

  // sum contributions to component accumulation from each phase
  for( localIndex ip = 0; ip < numPhases; ++ip )
  {
    localVolBalance -= phaseVolFrac[ip];
    localVolBalanceJacobian[0] -= dPhaseVolFrac_dPres[ip];

    for( localIndex jc = 0; jc < NC; ++jc )
    {
      localVolBalanceJacobian[jc + 1] -= dPhaseVolFrac_dComp[ip][jc];
    }
  }

  // scale saturation-based volume balance by pore volume (for better scaling w.r.t. other equations)
  for( localIndex idof = 0; idof < NC+1; ++idof )
  {
    localVolBalanceJacobian[idof] *= volume;
  }
  localVolBalance *= volume;
}

template< localIndex NC >
void
VolumeBalanceKernel::
  launch( localIndex const size,
          localIndex const numPhases,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< integer const > const & wellElemGhostRank,
          arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseVolFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & dWellElemPhaseVolFrac_dPres,
          arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dWellElemPhaseVolFrac_dComp,
          arrayView1d< real64 const > const & wellElemVolume,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
  {

    if( wellElemGhostRank[iwelem] >= 0 )
    {
      return;
    }

    real64 localVolBalance = 1.0;
    real64 localVolBalanceJacobian[NC+1]{};

    compute< NC >( numPhases,
                   wellElemVolume[iwelem],
                   wellElemPhaseVolFrac[iwelem],
                   dWellElemPhaseVolFrac_dPres[iwelem],
                   dWellElemPhaseVolFrac_dComp[iwelem],
                   localVolBalance,
                   localVolBalanceJacobian );

    // get equation/dof indices
    localIndex const localVolBalanceEqnIndex = wellElemDofNumber[iwelem] - rankOffset + ROFFSET::MASSBAL + NC;
    globalIndex localVolBalanceDOF[NC+1]{};
    for( localIndex jdof = 0; jdof < NC+1; ++jdof )
    {
      localVolBalanceDOF[jdof] = wellElemDofNumber[iwelem] + COFFSET::DPRES + jdof;
    }

    localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localVolBalanceEqnIndex,
                                                              localVolBalanceDOF,
                                                              localVolBalanceJacobian,
                                                              NC+1 );
    localRhs[localVolBalanceEqnIndex] += localVolBalance;
  } );
}

#define INST_VolumeBalanceKernel( NC ) \
  template \
  void VolumeBalanceKernel:: \
    launch< NC >( localIndex const size, \
                  localIndex const numPhases, \
                  globalIndex const rankOffset, \
                  arrayView1d< globalIndex const > const & wellElemDofNumber, \
                  arrayView1d< integer const > const & wellElemGhostRank, \
                  arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseVolFrac, \
                  arrayView2d< real64 const, compflow::USD_PHASE > const & dWellElemPhaseVolFrac_dPres, \
                  arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dWellElemPhaseVolFrac_dComp, \
                  arrayView1d< real64 const > const & wellElemVolume, \
                  CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                  arrayView1d< real64 > const & localRhs )

INST_VolumeBalanceKernel( 1 );
INST_VolumeBalanceKernel( 2 );
INST_VolumeBalanceKernel( 3 );
INST_VolumeBalanceKernel( 4 );
INST_VolumeBalanceKernel( 5 );

/******************************** PresTempCompFracInitializationKernel ********************************/

void
PresTempCompFracInitializationKernel::
  launch( localIndex const perforationSize,
          localIndex const subRegionSize,
          localIndex const numComps,
          localIndex const numPhases,
          localIndex const numPerforations,
          WellControls const & wellControls,
          real64 const & currentTime,
          ElementViewConst< arrayView1d< real64 const > > const & resPres,
          ElementViewConst< arrayView1d< real64 const > > const & resTemp,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_COMP > > const & resCompDens,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & resPhaseVolFrac,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & resPhaseMassDens,
          arrayView1d< localIndex const > const & resElementRegion,
          arrayView1d< localIndex const > const & resElementSubRegion,
          arrayView1d< localIndex const > const & resElementIndex,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< real64 > const & wellElemPres,
          arrayView1d< real64 > const & wellElemTemp,
          arrayView2d< real64, compflow::USD_COMP > const & wellElemCompFrac )
{
  localIndex constexpr MAX_NUM_COMP = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;

  real64 const targetBHP = wellControls.getTargetBHP( currentTime );
  real64 const refWellElemGravCoef = wellControls.getReferenceGravityCoef();
  WellControls::Control const currentControl = wellControls.getControl();
  WellControls::Type const wellType = wellControls.getType();

  // loop over all perforations to compute an average total mass density and component fraction
  RAJA::ReduceSum< parallelDeviceReduce, real64 > sumTotalMassDens( 0 );
  RAJA::ReduceSum< parallelDeviceReduce, real64 > sumTemp( 0 );
  RAJA::ReduceMin< parallelDeviceReduce, real64 > minResPres( 1e10 );
  RAJA::ReduceMax< parallelDeviceReduce, real64 > maxResPres( 0 );
  forAll< parallelDevicePolicy<> >( perforationSize, [=] GEOSX_HOST_DEVICE ( localIndex const iperf )
  {
    // get the reservoir (sub)region and element indices
    localIndex const er = resElementRegion[iperf];
    localIndex const esr = resElementSubRegion[iperf];
    localIndex const ei = resElementIndex[iperf];

    minResPres.min( resPres[er][esr][ei] );
    maxResPres.max( resPres[er][esr][ei] );

    // increment the temperature
    sumTemp += resTemp[er][esr][ei];

    // increment the average total mass density
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      sumTotalMassDens += resPhaseVolFrac[er][esr][ei][ip] * resPhaseMassDens[er][esr][ei][0][ip];
    }
  } );

  // TODO: there must a better way to do what is below
  // I would like to define an array of RAJA::ReduceSum to be able to do sum[ic] += ...
  // and put back what is below in the previous kernel.
  stackArray1d< real64, MAX_NUM_COMP > sumCompFrac( numComps );
  for( localIndex ic = 0; ic < numComps; ++ic )
  {
    RAJA::ReduceSum< parallelDeviceReduce, real64 > sum( 0.0 );
    forAll< parallelDevicePolicy<> >( perforationSize, [=] GEOSX_HOST_DEVICE ( localIndex const iperf )
    {
      // get the reservoir (sub)region and element indices
      localIndex const er = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei = resElementIndex[iperf];

      real64 perfTotalDens = 0.0;
      for( localIndex jc = 0; jc < numComps; ++jc )
      {
        perfTotalDens += resCompDens[er][esr][ei][jc];
      }
      sum += resCompDens[er][esr][ei][ic] / perfTotalDens;
    } );
    sumCompFrac[ic] = sum.get();
  }

  real64 const pres = ( wellControls.getType() == WellControls::Type::PRODUCER )
                      ? MpiWrapper::min( minResPres.get() )
                      : MpiWrapper::max( maxResPres.get() );
  real64 const avgTotalMassDens = MpiWrapper::sum( sumTotalMassDens.get() ) / numPerforations;

  stackArray1d< real64, MAX_NUM_COMP > avgCompFrac( numComps );
  real64 avgTemp = 0;
  // compute average component fraction
  if( wellControls.getType() == WellControls::Type::PRODUCER )
  {
    // use average temperature from reservoir
    avgTemp = MpiWrapper::sum( sumTemp.get() ) / numPerforations;

    // use average comp frac from reservoir
    real64 compFracSum = 0;
    for( localIndex ic = 0; ic < numComps; ++ic )
    {
      avgCompFrac[ic] = MpiWrapper::sum( sumCompFrac[ic] ) / numPerforations;
      compFracSum += avgCompFrac[ic];
    }

    real64 const tol = 1e-13;
    GEOSX_THROW_IF( compFracSum < 1 - tol || compFracSum > 1 + tol,
                    "Invalid well initialization: sum of component fractions should be between 0 and 1",
                    InputError );
  }
  else // injector
  {
    // use temperature from XML file
    avgTemp = wellControls.getInjectionTemperature();

    // use comp frac from XML file
    for( localIndex ic = 0; ic < numComps; ++ic )
    {
      avgCompFrac[ic] = wellControls.getInjectionStream()[ic];
    }
  }

  // set the global component fractions to avgCompFrac / temperature to avgTemp
  forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
  {
    wellElemTemp[iwelem] = avgTemp;
    for( localIndex ic = 0; ic < numComps; ++ic )
    {
      wellElemCompFrac[iwelem][ic] = avgCompFrac[ic];
    }
  } );

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

  GEOSX_THROW_IF( pressureControl <= 0,
                  "Invalid well initialization: negative pressure was found",
                  InputError );

  // estimate the pressures in the well elements using this avgDens
  forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
  {
    wellElemPres[iwelem] = pressureControl + avgTotalMassDens * ( wellElemGravCoef[iwelem] - gravCoefControl );
  } );
}

/******************************** CompDensInitializationKernel ********************************/

void
CompDensInitializationKernel::
  launch( localIndex const subRegionSize,
          localIndex const numComponents,
          arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompFrac,
          arrayView2d< real64 const, multifluid::USD_FLUID > const & wellElemTotalDens,
          arrayView2d< real64, compflow::USD_COMP > const & wellElemCompDens )
{
  forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
  {
    for( localIndex ic = 0; ic < numComponents; ++ic )
    {
      wellElemCompDens[iwelem][ic] = wellElemCompFrac[iwelem][ic] * wellElemTotalDens[iwelem][0];
    }
  } );
}

/******************************** RateInitializationKernel ********************************/

void
RateInitializationKernel::
  launch( localIndex const subRegionSize,
          localIndex const targetPhaseIndex,
          WellControls const & wellControls,
          real64 const & currentTime,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens,
          arrayView2d< real64 const, multifluid::USD_FLUID > const & totalDens,
          arrayView1d< real64 > const & connRate )
{
  WellControls::Control const control = wellControls.getControl();
  WellControls::Type const wellType = wellControls.getType();
  real64 const targetTotalRate = wellControls.getTargetTotalRate( currentTime );
  real64 const targetPhaseRate = wellControls.getTargetPhaseRate( currentTime );

  // Estimate the connection rates
  forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
  {
    if( control == WellControls::Control::BHP )
    {
      // if BHP constraint set rate below the absolute max rate
      // with the appropriate sign (negative for prod, positive for inj)
      if( wellType == WellControls::Type::PRODUCER )
      {
        connRate[iwelem] = LvArray::math::max( 0.1 * targetPhaseRate * phaseDens[iwelem][0][targetPhaseIndex], -1e3 );
      }
      else
      {
        connRate[iwelem] = LvArray::math::min( 0.1 * targetTotalRate * totalDens[iwelem][0], 1e3 );
      }
    }
    else
    {
      if( wellType == WellControls::Type::PRODUCER )
      {
        connRate[iwelem] = targetPhaseRate * phaseDens[iwelem][0][targetPhaseIndex];
      }
      else
      {
        connRate[iwelem] = targetTotalRate * totalDens[iwelem][0];
      }
    }
  } );
}

/******************************** TotalMassDensityKernel ****************************/

template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
TotalMassDensityKernel::
  compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dPhaseVolFrac_dCompDens,
           arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dCompFrac_dCompDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseMassDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & dPhaseMassDens_dPres,
           arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dPhaseMassDens_dComp,
           real64 & totalMassDens,
           real64 & dTotalMassDens_dPres,
           arraySlice1d< real64, compflow::USD_FLUID_DC - 1 > const & dTotalMassDens_dCompDens )
{
  real64 dMassDens_dC[NC]{};

  totalMassDens = 0.0;
  dTotalMassDens_dPres = 0.0;
  for( localIndex ic = 0; ic < NC; ++ic )
  {
    dTotalMassDens_dCompDens[ic] = 0.0;
  }

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    totalMassDens += phaseVolFrac[ip] * phaseMassDens[ip];
    dTotalMassDens_dPres += dPhaseVolFrac_dPres[ip] * phaseMassDens[ip]
                            + phaseVolFrac[ip] * dPhaseMassDens_dPres[ip];

    applyChainRule( NC, dCompFrac_dCompDens, dPhaseMassDens_dComp[ip], dMassDens_dC );
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dTotalMassDens_dCompDens[ic] += dPhaseVolFrac_dCompDens[ip][ic] * phaseMassDens[ip]
                                      + phaseVolFrac[ip] * dMassDens_dC[ic];
    }
  }
}

template< localIndex NC, localIndex NP >
void
TotalMassDensityKernel::
  launch( localIndex const size,
          arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseVolFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & dWellElemPhaseVolFrac_dPres,
          arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dWellElemPhaseVolFrac_dCompDens,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dWellElemCompFrac_dCompDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & wellElemPhaseMassDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dWellElemPhaseMassDens_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dWellElemPhaseMassDens_dComp,
          arrayView1d< real64 > const & wellElemTotalMassDens,
          arrayView1d< real64 > const & dWellElemTotalMassDens_dPres,
          arrayView2d< real64, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens_dCompDens )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
  {
    compute< NC, NP >( wellElemPhaseVolFrac[iwelem],
                       dWellElemPhaseVolFrac_dPres[iwelem],
                       dWellElemPhaseVolFrac_dCompDens[iwelem],
                       dWellElemCompFrac_dCompDens[iwelem],
                       wellElemPhaseMassDens[iwelem][0],
                       dWellElemPhaseMassDens_dPres[iwelem][0],
                       dWellElemPhaseMassDens_dComp[iwelem][0],
                       wellElemTotalMassDens[iwelem],
                       dWellElemTotalMassDens_dPres[iwelem],
                       dWellElemTotalMassDens_dCompDens[iwelem] );
  } );
}

#define INST_TotalMassDensityKernel( NC, NP ) \
  template \
  void TotalMassDensityKernel:: \
    launch< NC, NP >( localIndex const size, \
                      arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseVolFrac, \
                      arrayView2d< real64 const, compflow::USD_PHASE > const & dWellElemPhaseVolFrac_dPres, \
                      arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dWellElemPhaseVolFrac_dCompDens, \
                      arrayView3d< real64 const, compflow::USD_COMP_DC > const & dWellElemCompFrac_dCompDens, \
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & wellElemPhaseMassDens, \
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & dWellElemPhaseMassDens_dPres, \
                      arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dWellElemPhaseMassDens_dComp, \
                      arrayView1d< real64 > const & wellElemTotalMassDens, \
                      arrayView1d< real64 > const & dWellElemTotalMassDens_dPres, \
                      arrayView2d< real64, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens_dCompDens )

INST_TotalMassDensityKernel( 1, 2 );
INST_TotalMassDensityKernel( 2, 2 );
INST_TotalMassDensityKernel( 3, 2 );
INST_TotalMassDensityKernel( 4, 2 );
INST_TotalMassDensityKernel( 5, 2 );
INST_TotalMassDensityKernel( 1, 3 );
INST_TotalMassDensityKernel( 2, 3 );
INST_TotalMassDensityKernel( 3, 3 );
INST_TotalMassDensityKernel( 4, 3 );
INST_TotalMassDensityKernel( 5, 3 );


} // end namespace CompositionalMultiphaseWellKernels

} // end namespace geosx
