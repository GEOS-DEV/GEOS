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
 * @file CompositionalMultiphaseWellKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELLKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELLKERNELS_HPP

#include "common/DataTypes.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"

namespace geosx
{

namespace CompositionalMultiphaseWellKernels
{

static constexpr real64 minDensForDivision = 1e-10;

/******************************** ControlEquationHelper ********************************/

struct ControlEquationHelper
{

  GEOSX_HOST_DEVICE
  static void
  Switch( WellControls::Type const & wellType,
          WellControls::Control const & currentControl,
          real64 const & targetBHP,
          real64 const & targetConnRate,
          real64 const & wellElemPressure,
          real64 const & dWellElemPressure,
          real64 const & connRate,
          real64 const & dConnRate,
          WellControls::Control & newControl )
  {
    // TODO: check all inactive constraints (possibly more than one) and switch the one which is most violated
    // TODO: for the rate, use surface conditions (flash for compositional, easier for BO)

    // if isViable is true at the end of the following checks, no need to switch
    bool controlIsViable = false;

    real64 const refRate = connRate + dConnRate;
    real64 const refPressure = wellElemPressure + dWellElemPressure;

    // BHP control
    if( currentControl == WellControls::Control::BHP )
    {
      // the control is viable if the reference rate is below the max rate
      controlIsViable = ( fabs( refRate ) <= fabs( targetConnRate ) );
    }
    else // rate control
    {
      // the control is viable if the reference pressure is below/above the max/min pressure
      if( wellType == WellControls::Type::PRODUCER )
      {
        // targetBHP specifies a min pressure here
        controlIsViable = ( refPressure >= targetBHP );
      }
      else
      {
        // targetBHP specifies a max pressure here
        controlIsViable = ( refPressure <= targetBHP );
      }
    }

    if( controlIsViable )
    {
      newControl = currentControl;
    }
    else
    {
      newControl = ( currentControl == WellControls::Control::BHP )
                 ? WellControls::Control::LIQUIDRATE
                 : WellControls::Control::BHP;
    }
  }

  GEOSX_HOST_DEVICE
  static void
  Compute( globalIndex const rankOffset,
           localIndex const numComponents,
           WellControls::Control const currentControl,
           real64 const & targetBHP,
           real64 const & targetConnRate,
           globalIndex const wellElemDofNumber,
           real64 const & wellElemPressure,
           real64 const & dWellElemPressure,
           real64 const & connRate,
           real64 const & dConnRate,
           CRSMatrixView< real64, globalIndex const > const & localMatrix,
           arrayView1d< real64 > const & localRhs )
  {
    globalIndex eqnRowIndex = 0;
    globalIndex dofColIndex = 0;
    real64 controlEqn = 0;
    real64 dControlEqn_dX = 0;

    // BHP control
    if( currentControl == WellControls::Control::BHP )
    {

      // get the pressure and compute normalizer
      real64 const currentBHP = wellElemPressure + dWellElemPressure;

      // control equation is a normalized difference
      // between current pressure and target pressure
      controlEqn = ( currentBHP - targetBHP );
      dControlEqn_dX = 1.0;
      dofColIndex = wellElemDofNumber + CompositionalMultiphaseWell::ColOffset::DPRES;
      eqnRowIndex = wellElemDofNumber + CompositionalMultiphaseWell::RowOffset::CONTROL - rankOffset;
    }
    else if( currentControl == WellControls::Control::LIQUIDRATE ) // liquid rate control
    {
      // get rates and compute normalizer
      real64 const currentConnRate = connRate + dConnRate;

      // control equation is a normalized difference
      // between current rate and target rate
      controlEqn = ( currentConnRate - targetConnRate );
      dControlEqn_dX = 1.0;
      dofColIndex = wellElemDofNumber + CompositionalMultiphaseWell::ColOffset::DCOMP + numComponents;
      eqnRowIndex = wellElemDofNumber + CompositionalMultiphaseWell::RowOffset::CONTROL - rankOffset;
    }
    else
    {
      GEOSX_ERROR_IF( ( currentControl != WellControls::Control::BHP ) && ( currentControl != WellControls::Control::LIQUIDRATE ),
                      "Phase rate constraints for CompositionalMultiphaseWell will be implemented later" );
    }

    localMatrix.addToRow< serialAtomic >( eqnRowIndex,
                                          &dofColIndex,
                                          &dControlEqn_dX,
                                          1 );
    localRhs[eqnRowIndex] += controlEqn;
  }

};

/******************************** FluxKernel ********************************/

struct FluxKernel
{

  template< typename POLICY >
  static void
  Launch( localIndex const size,
          globalIndex const rankOffset,
          localIndex const numComponents,
          localIndex const numDofPerResElement,
          WellControls const & wellControls,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & connRate,
          arrayView1d< real64 const > const & dConnRate,
          arrayView2d< real64 const > const & wellElemCompFrac,
          arrayView3d< real64 const > const & dWellElemCompFrac_dCompDens,
          real64 const & dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
  {
    localIndex const NC = numComponents;
    localIndex const resNDOF = numDofPerResElement;

    WellControls::Type const wellType = wellControls.GetType();
    arrayView1d< real64 const > const & injection = wellControls.GetInjectionStream();

    // loop over the well elements to compute the fluxes between elements
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      localIndex constexpr maxNumComp = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;
      localIndex constexpr maxNumDof  = maxNumComp + 1;

      // create local work arrays
      stackArray1d< real64, maxNumComp > compFracUp( NC );
      stackArray2d< real64, maxNumComp * maxNumComp > dCompFrac_dCompDensUp( NC, NC );

      stackArray1d< real64, maxNumComp > compFlux( NC );
      stackArray1d< real64, maxNumComp > dCompFlux_dRate( NC );
      stackArray1d< real64, maxNumComp > dCompFlux_dPresUp( NC );
      stackArray2d< real64, maxNumComp * maxNumComp > dCompFlux_dCompDensUp( NC, NC );

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
        compFlux[ic] = compFracUp[ic] * currentConnRate;
        dCompFlux_dRate[ic] = compFracUp[ic];
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
        stackArray1d< real64, maxNumComp > oneSidedFlux( NC );
        stackArray2d< real64, maxNumComp > oneSidedFluxJacobian_dRate( NC, 1 );
        stackArray2d< real64, maxNumComp * maxNumDof > oneSidedFluxJacobian_dPresCompUp( NC, resNDOF );

        stackArray1d< globalIndex, maxNumComp > oneSidedEqnRowIndices( NC );
        stackArray1d< globalIndex, maxNumDof > oneSidedDofColIndices_dPresCompUp( resNDOF );
        globalIndex oneSidedDofColIndices_dRate = 0;

        // flux terms
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          oneSidedFlux[ic] = -dt * compFlux[ic];

          // derivative with respect to rate
          oneSidedFluxJacobian_dRate( ic, 0 ) = -dt * dCompFlux_dRate[ic];

          // derivative with respect to upstream pressure
          oneSidedFluxJacobian_dPresCompUp[ic][0] = -dt * dCompFlux_dPresUp[ic];

          // derivatives with respect to upstream component densities
          for( localIndex jdof = 0; jdof < NC; ++jdof )
          {
            oneSidedFluxJacobian_dPresCompUp[ic][jdof + 1] = -dt * dCompFlux_dCompDensUp[ic][jdof];
          }

        }

        // jacobian indices
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          // mass balance equations for all components
          oneSidedEqnRowIndices[ic] = offsetUp + CompositionalMultiphaseWell::RowOffset::MASSBAL + ic - rankOffset;
        }

        // in the dof ordering used in this class, there are 1 pressure dofs
        // and NC compDens dofs before the rate dof in this block
        localIndex const dRateColOffset = CompositionalMultiphaseWell::ColOffset::DCOMP + NC;
        oneSidedDofColIndices_dRate = offsetCurrent + dRateColOffset;

        for( localIndex jdof = 0; jdof < resNDOF; ++jdof )
        {
          // dofs are the **upstream** pressure and component densities
          oneSidedDofColIndices_dPresCompUp[jdof] = offsetUp + CompositionalMultiphaseWell::ColOffset::DPRES + jdof;
        }

        for( localIndex i = 0; i < oneSidedFlux.size(); ++i )
        {
          if( oneSidedEqnRowIndices[i] >= 0 && oneSidedEqnRowIndices[i] < localMatrix.numRows() )
          {
            localMatrix.addToRow< parallelDeviceAtomic >( oneSidedEqnRowIndices[i],
                                                          &oneSidedDofColIndices_dRate,
                                                          oneSidedFluxJacobian_dRate.data() + i,
                                                          1 );
            localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( oneSidedEqnRowIndices[i],
                                                                              oneSidedDofColIndices_dPresCompUp.data(),
                                                                              oneSidedFluxJacobian_dPresCompUp.data() + i * resNDOF,
                                                                              resNDOF );
            atomicAdd( parallelDeviceAtomic{}, &localRhs[oneSidedEqnRowIndices[i]], oneSidedFlux[i] );
          }
        }
      }
      else // not an exit connection
      {
        stackArray1d< real64, 2 * maxNumComp > localFlux( 2 * NC );
        stackArray2d< real64, 2 * maxNumComp > localFluxJacobian_dRate( 2 * NC, 1 );
        stackArray2d< real64, 2 * maxNumComp * maxNumDof > localFluxJacobian_dPresCompUp( 2 * NC, resNDOF );

        stackArray1d< globalIndex, 2 * maxNumComp > eqnRowIndices( 2 * NC );
        stackArray1d< globalIndex, maxNumDof > dofColIndices_dPresCompUp( resNDOF );
        globalIndex dofColIndices_dRate = 0;

        globalIndex const offsetNext = wellElemDofNumber[iwelemNext];

        // flux terms
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          localFlux[WellSolverBase::ElemTag::NEXT * NC + ic] = dt * compFlux[ic];
          localFlux[WellSolverBase::ElemTag::CURRENT * NC + ic] = -dt * compFlux[ic];

          // derivative with respect to rate
          localFluxJacobian_dRate( WellSolverBase::ElemTag::NEXT * NC + ic, 0 ) = dt * dCompFlux_dRate[ic];
          localFluxJacobian_dRate( WellSolverBase::ElemTag::CURRENT * NC + ic, 0 ) = -dt * dCompFlux_dRate[ic];

          // derivative with respect to upstream pressure
          localFluxJacobian_dPresCompUp[WellSolverBase::ElemTag::NEXT * NC + ic][0] = dt * dCompFlux_dPresUp[ic];
          localFluxJacobian_dPresCompUp[WellSolverBase::ElemTag::CURRENT * NC + ic][0] = -dt * dCompFlux_dPresUp[ic];

          // derivatives with respect to upstream component densities
          for( localIndex jdof = 0; jdof < NC; ++jdof )
          {
            localFluxJacobian_dPresCompUp[WellSolverBase::ElemTag::NEXT * NC + ic][jdof + 1] =
              dt * dCompFlux_dCompDensUp[ic][jdof];
            localFluxJacobian_dPresCompUp[WellSolverBase::ElemTag::CURRENT * NC + ic][jdof + 1] =
              -dt * dCompFlux_dCompDensUp[ic][jdof];
          }
        }

        // jacobian indices
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          // mass balance equations for all components
          eqnRowIndices[WellSolverBase::ElemTag::NEXT * NC + ic] = offsetNext
                                                                   + CompositionalMultiphaseWell::RowOffset::MASSBAL + ic - rankOffset;
          eqnRowIndices[WellSolverBase::ElemTag::CURRENT * NC + ic] = offsetCurrent
                                                                      + CompositionalMultiphaseWell::RowOffset::MASSBAL + ic - rankOffset;
        }

        // in the dof ordering used in this class, there are 1 pressure dofs
        // and NC compDens dofs before the rate dof in this block
        localIndex const dRateColOffset = CompositionalMultiphaseWell::ColOffset::DCOMP + NC;
        dofColIndices_dRate = offsetCurrent + dRateColOffset;

        for( localIndex jdof = 0; jdof < resNDOF; ++jdof )
        {
          // dofs are the **upstream** pressure and component densities
          dofColIndices_dPresCompUp[jdof] = offsetUp + CompositionalMultiphaseWell::ColOffset::DPRES + jdof;
        }

        for( localIndex i = 0; i < localFlux.size(); ++i )
        {
          if( eqnRowIndices[i] >= 0 && eqnRowIndices[i] < localMatrix.numRows() )
          {
            localMatrix.addToRow< parallelDeviceAtomic >( eqnRowIndices[i],
                                                          &dofColIndices_dRate,
                                                          localFluxJacobian_dRate.data() + i,
                                                          1 );
            localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndices[i],
                                                                              dofColIndices_dPresCompUp.data(),
                                                                              localFluxJacobian_dPresCompUp.data() + i * resNDOF,
                                                                              resNDOF );
            atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndices[i]], localFlux[i] );
          }
        }
      }
    } );
  }

};

/******************************** PressureRelationKernel ********************************/

struct PressureRelationKernel
{

  template< typename POLICY, typename REDUCE_POLICY >
  static localIndex
  Launch( localIndex const size,
          globalIndex const rankOffset,
          bool const isLocallyOwned,
          localIndex const numComponents,
          localIndex const numDofPerResElement,
          WellControls const & wellControls,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & connRate,
          arrayView1d< real64 const > const & dConnRate,
          arrayView1d< real64 const > const & wellElemPressure,
          arrayView1d< real64 const > const & dWellElemPressure,
          arrayView1d< real64 const > const & wellElemTotalMassDens,
          arrayView1d< real64 const > const & dWellElemTotalMassDens_dPres,
          arrayView2d< real64 const > const & dWellElemTotalMassDens_dCompDens,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
  {
    localIndex const NC = numComponents;
    localIndex const resNDOF = numDofPerResElement;

    real64 const targetBHP = wellControls.GetTargetBHP();
    real64 const targetRate = wellControls.GetTargetRate();
    WellControls::Control const currentControl = wellControls.GetControl();
    WellControls::Type const wellType = wellControls.GetType();
    localIndex const iwelemControl = wellControls.GetReferenceWellElementIndex();

    RAJA::ReduceMax< REDUCE_POLICY, localIndex > switchControl( 0 );

    // loop over the well elements to compute the pressure relations between well elements
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      localIndex const iwelemNext = nextWellElemIndex[iwelem];

      if( iwelemNext < 0 && isLocallyOwned ) // if iwelemNext < 0, form control equation
      {

        WellControls::Control newControl = currentControl;
        ControlEquationHelper::Switch( wellType,
                                       currentControl,
                                       targetBHP,
                                       targetRate,
                                       wellElemPressure[iwelemControl],
                                       dWellElemPressure[iwelemControl],
                                       connRate[iwelemControl],
                                       dConnRate[iwelemControl],
                                       newControl );
        if( currentControl != newControl )
        {
          switchControl.max( 1 );
        }

        ControlEquationHelper::Compute( rankOffset,
                                        NC,
                                        newControl,
                                        targetBHP,
                                        targetRate,
                                        wellElemDofNumber[iwelemControl],
                                        wellElemPressure[iwelemControl],
                                        dWellElemPressure[iwelemControl],
                                        connRate[iwelemControl],
                                        dConnRate[iwelemControl],
                                        localMatrix,
                                        localRhs );
      }
      else if( iwelemNext >= 0 ) // if iwelemNext >= 0, form momentum equation
      {
        localIndex constexpr maxNumComp = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;
        localIndex constexpr maxNumDof  = maxNumComp + 1;

        // local working variables and arrays
        stackArray1d< globalIndex, 2 * maxNumDof > dofColIndices( 2 * resNDOF );
        stackArray1d< real64, 2 * maxNumDof > localPresRelJacobian( 2 * resNDOF );

        stackArray1d< real64, maxNumComp > dAvgMassDens_dCompCurrent( NC );
        stackArray1d< real64, maxNumComp > dAvgMassDens_dCompNext( NC );

        // compute the average density at the interface between well elements
        real64 const avgMassDens = 0.5 * ( wellElemTotalMassDens[iwelemNext] + wellElemTotalMassDens[iwelem] );
        real64 const dAvgMassDens_dPresNext = 0.5 * dWellElemTotalMassDens_dPres[iwelemNext];
        real64 const dAvgMassDens_dPresCurrent = 0.5 * dWellElemTotalMassDens_dPres[iwelem];
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dAvgMassDens_dCompNext[ic] = 0.5 * dWellElemTotalMassDens_dCompDens[iwelemNext][ic];
          dAvgMassDens_dCompCurrent[ic] = 0.5 * dWellElemTotalMassDens_dCompDens[iwelem][ic];
        }

        // compute depth diff times acceleration
        real64 const gravD = wellElemGravCoef[iwelemNext] - wellElemGravCoef[iwelem];

        // compute the current pressure in the two well elements
        real64 const pressureNext = wellElemPressure[iwelemNext] + dWellElemPressure[iwelemNext];
        real64 const pressureCurrent = wellElemPressure[iwelem] + dWellElemPressure[iwelem];

        // compute momentum flux and derivatives
        localIndex const localDofIndexPresNext = WellSolverBase::ElemTag::NEXT * resNDOF;
        localIndex const localDofIndexPresCurrent = WellSolverBase::ElemTag::CURRENT * resNDOF;

        globalIndex const offsetNext = wellElemDofNumber[iwelemNext];
        globalIndex const offsetCurrent = wellElemDofNumber[iwelem];

        globalIndex const eqnRowIndex = offsetCurrent + CompositionalMultiphaseWell::RowOffset::CONTROL - rankOffset;
        dofColIndices[localDofIndexPresNext] = offsetNext + CompositionalMultiphaseWell::ColOffset::DPRES;
        dofColIndices[localDofIndexPresCurrent] = offsetCurrent + CompositionalMultiphaseWell::ColOffset::DPRES;

        real64 const localPresRel = ( pressureNext - pressureCurrent - avgMassDens * gravD );

        localPresRelJacobian[localDofIndexPresNext] = ( 1 - dAvgMassDens_dPresNext * gravD );
        localPresRelJacobian[localDofIndexPresCurrent] = ( -1 - dAvgMassDens_dPresCurrent * gravD );

        for( localIndex ic = 0; ic < NC; ++ic )
        {
          localIndex const localDofIndexCompNext = localDofIndexPresNext + ic + 1;
          localIndex const localDofIndexCompCurrent = localDofIndexPresCurrent + ic + 1;

          dofColIndices[localDofIndexCompNext] = offsetNext + CompositionalMultiphaseWell::ColOffset::DCOMP + ic;
          dofColIndices[localDofIndexCompCurrent] = offsetCurrent + CompositionalMultiphaseWell::ColOffset::DCOMP + ic;

          localPresRelJacobian[localDofIndexCompNext] = -dAvgMassDens_dCompNext[ic] * gravD;
          localPresRelJacobian[localDofIndexCompCurrent] = -dAvgMassDens_dCompCurrent[ic] * gravD;
        }

        // TODO: add friction and acceleration terms

        if( eqnRowIndex >= 0 && eqnRowIndex < localMatrix.numRows() )
        {
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndex,
                                                                            dofColIndices.data(),
                                                                            localPresRelJacobian.data(),
                                                                            2 * resNDOF );
          atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndex], localPresRel );
        }
      }
    } );
    return switchControl.get();
  }

};

/******************************** PerforationKernel ********************************/

struct PerforationKernel
{

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  template< typename POLICY >
  static void
  Launch( localIndex const size,
          localIndex const numComponents,
          localIndex const numPhases,
          ElementViewConst< arrayView1d< real64 const > > const & resPres,
          ElementViewConst< arrayView1d< real64 const > > const & dResPres,
          ElementViewConst< arrayView2d< real64 const > > const & resPhaseMob,
          ElementViewConst< arrayView2d< real64 const > > const & dResPhaseMob_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & dResPhaseMob_dComp,
          ElementViewConst< arrayView2d< real64 const > > const & dResPhaseVolFrac_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & dResPhaseVolFrac_dComp,
          ElementViewConst< arrayView3d< real64 const > > const & dResCompFrac_dCompDens,
          ElementViewConst< arrayView3d< real64 const > > const & resPhaseVisc,
          ElementViewConst< arrayView3d< real64 const > > const & dResPhaseVisc_dPres,
          ElementViewConst< arrayView4d< real64 const > > const & dResPhaseVisc_dComp,
          ElementViewConst< arrayView4d< real64 const > > const & resPhaseCompFrac,
          ElementViewConst< arrayView4d< real64 const > > const & dResPhaseCompFrac_dPres,
          ElementViewConst< arrayView5d< real64 const > > const & dResPhaseCompFrac_dComp,
          ElementViewConst< arrayView3d< real64 const > > const & resPhaseRelPerm,
          ElementViewConst< arrayView4d< real64 const > > const & dResPhaseRelPerm_dPhaseVolFrac,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< real64 const > const & wellElemPres,
          arrayView1d< real64 const > const & dWellElemPres,
          arrayView2d< real64 const > const & wellElemCompDens,
          arrayView2d< real64 const > const & dWellElemCompDens,
          arrayView1d< real64 const > const & wellElemTotalMassDens,
          arrayView1d< real64 const > const & dWellElemTotalMassDens_dPres,
          arrayView2d< real64 const > const & dWellElemTotalMassDens_dCompDens,
          arrayView2d< real64 const > const & wellElemCompFrac,
          arrayView3d< real64 const > const & dWellElemCompFrac_dCompDens,
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
    localIndex const NC = numComponents;
    localIndex const NP = numPhases;

    // loop over the perforations to compute the perforation rates
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const iperf )
    {
      localIndex constexpr maxNumComp = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;

      // local working variables and arrays
      stackArray1d< real64, maxNumComp > dPhaseCompFrac_dCompDens( NC );

      real64 pres[ 2 ] = { 0.0 };
      real64 dPres_dP[ 2 ] = { 0.0 };
      stackArray2d< real64, 2 * maxNumComp > dPres_dC( 2, NC );

      real64 dFlux_dP[ 2 ] = { 0.0 };
      stackArray2d< real64, 2 * maxNumComp > dFlux_dC( 2, NC );

      real64 dMult_dP[ 2 ] = { 0.0 };
      stackArray2d< real64, 2 * maxNumComp > dMult_dC( 2, NC );

      stackArray1d< real64, maxNumComp > dResTotalMob_dC( NC );

      stackArray2d< real64, 2 * maxNumComp > phaseCompFrac( 2, NC );
      stackArray2d< real64, 2 * maxNumComp > dPhaseCompFrac_dP( 2, NC );
      stackArray3d< real64, 2 * maxNumComp * maxNumComp > dPhaseCompFrac_dC( 2, NC, NC );

      stackArray1d< real64, maxNumComp > dVisc_dC( NC );
      stackArray1d< real64, maxNumComp > dRelPerm_dC( NC );

      real64 dPotDiff_dP[ 2 ] = { 0.0 };
      stackArray2d< real64, 2 * maxNumComp > dPotDiff_dC( 2, NC );

      real64 multiplier[ 2 ] = { 0.0 };

      // reset the perforation rates
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        compPerfRate[iperf][ic] = 0.0;
        for( localIndex ke = 0; ke < 2; ++ke )
        {
          dCompPerfRate_dPres[iperf][ke][ic] = 0.0;
          for( localIndex jc = 0; jc < NC; ++jc )
          {
            dCompPerfRate_dComp[iperf][ke][ic][jc] = 0.0;
          }
        }
      }

      // 1) copy the variables from the reservoir and well element

      // a) get reservoir variables

      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];
      // get the index of the well elem
      localIndex const iwelem = perfWellElemIndex[iperf];

      pres[CompositionalMultiphaseWell::SubRegionTag::RES] = resPres[er][esr][ei] + dResPres[er][esr][ei];
      dPres_dP[CompositionalMultiphaseWell::SubRegionTag::RES] = 1.0;

      // Here in the absence of a buoyancy term we assume that the reservoir cell is perforated at its center
      // TODO: add a buoyancy term for the reservoir side here

      multiplier[CompositionalMultiphaseWell::SubRegionTag::RES] = 1.0;

      // b) get well variables

      pres[CompositionalMultiphaseWell::SubRegionTag::WELL] = wellElemPres[iwelem] + dWellElemPres[iwelem];
      dPres_dP[CompositionalMultiphaseWell::SubRegionTag::WELL] = 1.0;

      multiplier[CompositionalMultiphaseWell::SubRegionTag::WELL] = -1.0;

      real64 const gravD = ( perfGravCoef[iperf] - wellElemGravCoef[iwelem] );

      pres[CompositionalMultiphaseWell::SubRegionTag::WELL] += wellElemTotalMassDens[iwelem] * gravD;
      dPres_dP[CompositionalMultiphaseWell::SubRegionTag::WELL] += dWellElemTotalMassDens_dPres[iwelem] * gravD;
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dPres_dC[CompositionalMultiphaseWell::SubRegionTag::WELL][ic] += dWellElemTotalMassDens_dCompDens[iwelem][ic] * gravD;
      }

      // get transmissibility at the interface
      real64 const trans = perfTrans[iperf];

      // 2) compute potential difference

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

      real64 flux = 0.0;

      // 3) upwinding

      if( potDiff >= 0 )  // ** reservoir cell is upstream **
      {

        // loop over phases, compute and upwind phase flux
        // and sum contributions to each component's perforation rate
        for( localIndex ip = 0; ip < NP; ++ip )
        {

          // compute the phase flux and derivatives using upstream cell mobility
          flux = resPhaseMob[er][esr][ei][ip] * potDiff;

          dFlux_dP[CompositionalMultiphaseWell::SubRegionTag::RES] =
            dResPhaseMob_dPres[er][esr][ei][ip] * potDiff
            + resPhaseMob[er][esr][ei][ip] * dPotDiff_dP[CompositionalMultiphaseWell::SubRegionTag::RES];

          dFlux_dP[CompositionalMultiphaseWell::SubRegionTag::WELL] =
            resPhaseMob[er][esr][ei][ip] * dPotDiff_dP[CompositionalMultiphaseWell::SubRegionTag::WELL];

          for( localIndex ic = 0; ic < NC; ++ic )
          {
            dFlux_dC[CompositionalMultiphaseWell::SubRegionTag::RES][ic] =
              dResPhaseMob_dComp[er][esr][ei][ip][ic] * potDiff
              + resPhaseMob[er][esr][ei][ip] * dPotDiff_dC[CompositionalMultiphaseWell::SubRegionTag::RES][ic];

            dFlux_dC[CompositionalMultiphaseWell::SubRegionTag::WELL][ic] =
              resPhaseMob[er][esr][ei][ip] * dPotDiff_dC[CompositionalMultiphaseWell::SubRegionTag::WELL][ic];
          }

          // increment component fluxes
          for( localIndex ic = 0; ic < NC; ++ic )
          {
            compPerfRate[iperf][ic] += flux * resPhaseCompFrac[er][esr][ei][0][ip][ic];

            dCompPerfRate_dPres[iperf][CompositionalMultiphaseWell::SubRegionTag::RES][ic] +=
              resPhaseCompFrac[er][esr][ei][0][ip][ic] * dFlux_dP[CompositionalMultiphaseWell::SubRegionTag::RES];

            dCompPerfRate_dPres[iperf][CompositionalMultiphaseWell::SubRegionTag::RES][ic] +=
              dResPhaseCompFrac_dPres[er][esr][ei][0][ip][ic] * flux;

            dCompPerfRate_dPres[iperf][CompositionalMultiphaseWell::SubRegionTag::WELL][ic] +=
              resPhaseCompFrac[er][esr][ei][0][ip][ic] * dFlux_dP[CompositionalMultiphaseWell::SubRegionTag::WELL];

            applyChainRule( NC,
                            dResCompFrac_dCompDens[er][esr][ei],
                            dResPhaseCompFrac_dComp[er][esr][ei][0][ip][ic],
                            dPhaseCompFrac_dCompDens );

            for( localIndex jc = 0; jc < NC; ++jc )
            {
              dCompPerfRate_dComp[iperf][CompositionalMultiphaseWell::SubRegionTag::RES][ic][jc] +=
                dFlux_dC[CompositionalMultiphaseWell::SubRegionTag::RES][jc]
                * resPhaseCompFrac[er][esr][ei][0][ip][ic];

              dCompPerfRate_dComp[iperf][CompositionalMultiphaseWell::SubRegionTag::RES][ic][jc] +=
                flux * dPhaseCompFrac_dCompDens[jc];

              dCompPerfRate_dComp[iperf][CompositionalMultiphaseWell::SubRegionTag::WELL][ic][jc] +=
                dFlux_dC[CompositionalMultiphaseWell::SubRegionTag::WELL][jc]
                * resPhaseCompFrac[er][esr][ei][0][ip][ic];
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
          wellElemTotalDens += wellElemCompDens[iwelem][ic] + dWellElemCompDens[iwelem][ic];
        }

        // first, compute the reservoir total mobility (excluding phase density)
        for( localIndex ip = 0; ip < NP; ++ip )
        {
          // viscosity
          real64 const resVisc = resPhaseVisc[er][esr][ei][0][ip];
          real64 const dResVisc_dP  = dResPhaseVisc_dPres[er][esr][ei][0][ip];
          applyChainRule( NC, dResCompFrac_dCompDens[er][esr][ei],
                          dResPhaseVisc_dComp[er][esr][ei][0][ip],
                          dVisc_dC );

          // relative permeability
          real64 const resRelPerm = resPhaseRelPerm[er][esr][ei][0][ip];
          real64 dResRelPerm_dP = 0.0;
          for( localIndex jc = 0; jc < NC; ++jc )
          {
            dRelPerm_dC[jc] = 0;
          }

          for( localIndex jp = 0; jp < NP; ++jp )
          {
            real64 const dResRelPerm_dS = dResPhaseRelPerm_dPhaseVolFrac[er][esr][ei][0][ip][jp];
            dResRelPerm_dP += dResRelPerm_dS * dResPhaseVolFrac_dPres[er][esr][ei][jp];

            for( localIndex jc = 0; jc < NC; ++jc )
            {
              dRelPerm_dC[jc] += dResRelPerm_dS * dResPhaseVolFrac_dComp[er][esr][ei][jp][jc];
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
        real64 const mult = wellElemTotalDens * resTotalMob;
        dMult_dP[CompositionalMultiphaseWell::SubRegionTag::RES] = wellElemTotalDens * dResTotalMob_dP;
        dMult_dP[CompositionalMultiphaseWell::SubRegionTag::WELL] = 0; // because totalDens does not depend on pressure

        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dMult_dC[CompositionalMultiphaseWell::SubRegionTag::RES][ic] =
            wellElemTotalDens * dResTotalMob_dC[ic];

          dMult_dC[CompositionalMultiphaseWell::SubRegionTag::WELL][ic] = resTotalMob;
        }

        // compute the volumetric flux and derivatives using upstream cell mobility
        flux = mult * potDiff;

        dFlux_dP[CompositionalMultiphaseWell::SubRegionTag::RES] =
          dMult_dP[CompositionalMultiphaseWell::SubRegionTag::RES] * potDiff
          + mult * dPotDiff_dP[CompositionalMultiphaseWell::SubRegionTag::RES];

        dFlux_dP[CompositionalMultiphaseWell::SubRegionTag::WELL] =
          dMult_dP[CompositionalMultiphaseWell::SubRegionTag::WELL] * potDiff
          + mult * dPotDiff_dP[CompositionalMultiphaseWell::SubRegionTag::WELL];

        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dFlux_dC[CompositionalMultiphaseWell::SubRegionTag::RES][ic] =
            dMult_dC[CompositionalMultiphaseWell::SubRegionTag::RES][ic] * potDiff
            + mult * dPotDiff_dC[CompositionalMultiphaseWell::SubRegionTag::RES][ic];

          dFlux_dC[CompositionalMultiphaseWell::SubRegionTag::WELL][ic] =
            dMult_dC[CompositionalMultiphaseWell::SubRegionTag::WELL][ic] * potDiff
            + mult * dPotDiff_dC[CompositionalMultiphaseWell::SubRegionTag::WELL][ic];
        }

        // compute component fluxes
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          compPerfRate[iperf][ic] += wellElemCompFrac[iwelem][ic] * flux;

          dCompPerfRate_dPres[iperf][CompositionalMultiphaseWell::SubRegionTag::RES][ic] =
            wellElemCompFrac[iwelem][ic] * dFlux_dP[CompositionalMultiphaseWell::SubRegionTag::RES];

          dCompPerfRate_dPres[iperf][CompositionalMultiphaseWell::SubRegionTag::WELL][ic] =
            wellElemCompFrac[iwelem][ic] * dFlux_dP[CompositionalMultiphaseWell::SubRegionTag::WELL];

          for( localIndex jc = 0; jc < NC; ++jc )
          {
            dCompPerfRate_dComp[iperf][CompositionalMultiphaseWell::SubRegionTag::RES][ic][jc] +=
              wellElemCompFrac[iwelem][ic] * dFlux_dC[CompositionalMultiphaseWell::SubRegionTag::RES][jc];

            dCompPerfRate_dComp[iperf][CompositionalMultiphaseWell::SubRegionTag::WELL][ic][jc] +=
              wellElemCompFrac[iwelem][ic] * dFlux_dC[CompositionalMultiphaseWell::SubRegionTag::WELL][jc];

            dCompPerfRate_dComp[iperf][CompositionalMultiphaseWell::SubRegionTag::WELL][ic][jc] +=
              dWellElemCompFrac_dCompDens[iwelem][ic][jc] * flux;
          }
        }
      }
    } );
  }

};

/******************************** VolumeBalanceKernel ********************************/

struct VolumeBalanceKernel
{

  template< typename POLICY >
  static void
  Launch( localIndex const size,
          localIndex const numComponents,
          localIndex const numPhases,
          localIndex const numDofPerWellElement,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< integer const > const & wellElemGhostRank,
          arrayView2d< real64 const > const & wellElemPhaseVolFrac,
          arrayView2d< real64 const > const & dWellElemPhaseVolFrac_dPres,
          arrayView3d< real64 const > const & dWellElemPhaseVolFrac_dComp,
          arrayView1d< real64 const > const & wellElemVolume,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
  {
    localIndex constexpr maxNumComp = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;
    localIndex constexpr maxNumDof  = maxNumComp + 1;

    localIndex const NC        = numComponents;
    localIndex const NP        = numPhases;
    localIndex const welemNDOF = numDofPerWellElement;

    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {

      if( wellElemGhostRank[iwelem] >= 0 )
      {
        return;
      }

      stackArray1d< globalIndex, maxNumDof > localVolBalanceDOF( welemNDOF );
      stackArray1d< real64, maxNumDof > localVolBalanceJacobian( welemNDOF );

      // get equation/dof indices
      globalIndex const offset = wellElemDofNumber[iwelem];
      localIndex const volBalRowOffset = CompositionalMultiphaseWell::RowOffset::MASSBAL + NC;
      globalIndex const localVolBalanceEqnIndex = LvArray::integerConversion< localIndex >( offset - rankOffset ) + volBalRowOffset;
      for( localIndex jdof = 0; jdof < welemNDOF; ++jdof )
      {
        localVolBalanceDOF[jdof] = offset + CompositionalMultiphaseWell::ColOffset::DPRES + jdof;
      }

      real64 localVolBalance = 1.0;

      // sum contributions to component accumulation from each phase
      for( localIndex ip = 0; ip < NP; ++ip )
      {
        localVolBalance -= wellElemPhaseVolFrac[iwelem][ip];
        localVolBalanceJacobian[0] -= dWellElemPhaseVolFrac_dPres[iwelem][ip];

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          localVolBalanceJacobian[jc + 1] -= dWellElemPhaseVolFrac_dComp[iwelem][ip][jc];
        }
      }

      // scale saturation-based volume balance by pore volume (for better scaling w.r.t. other equations)
      for( localIndex idof = 0; idof < welemNDOF; ++idof )
      {
        localVolBalanceJacobian[idof] *= wellElemVolume[iwelem];
      }
      localVolBalance *= wellElemVolume[iwelem];

      localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localVolBalanceEqnIndex,
                                                                localVolBalanceDOF.data(),
                                                                localVolBalanceJacobian.data(),
                                                                welemNDOF );
      localRhs[localVolBalanceEqnIndex] += localVolBalance;
    } );
  }

};

/******************************** PresCompFracInitializationKernel ********************************/

struct PresCompFracInitializationKernel
{

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  template< typename POLICY >
  static void
  Launch( localIndex const perforationSize,
          localIndex const subRegionSize,
          localIndex const numComponents,
          localIndex const numPhases,
          bool const isLocallyOwned,
          int const topRank,
          localIndex const numPerforations,
          WellControls const & wellControls,
          ElementViewConst< arrayView1d< real64 const > > const & resPressure,
          ElementViewConst< arrayView2d< real64 const > > const & resCompDens,
          ElementViewConst< arrayView2d< real64 const > > const & resPhaseVolFrac,
          ElementViewConst< arrayView3d< real64 const > > const & resPhaseMassDens,
          arrayView1d< localIndex const > const & resElementRegion,
          arrayView1d< localIndex const > const & resElementSubRegion,
          arrayView1d< localIndex const > const & resElementIndex,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< real64 > const & wellElemPressure,
          arrayView2d< real64 > const & wellElemCompFrac )
  {
    localIndex constexpr maxNumComp = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;
    localIndex const NC = numComponents;
    localIndex const NP = numPhases;

    real64 const targetBHP = wellControls.GetTargetBHP();
    WellControls::Control const currentControl = wellControls.GetControl();
    WellControls::Type const wellType = wellControls.GetType();
    localIndex const iwelemControl = wellControls.GetReferenceWellElementIndex();

    // loop over all perforations to compute an average total mass density and component fraction
    RAJA::ReduceSum< parallelDeviceReduce, real64 > sumTotalMassDens( 0 );
    RAJA::ReduceMin< parallelDeviceReduce, real64 > minResPressure( 1e10 );
    RAJA::ReduceMax< parallelDeviceReduce, real64 > maxResPressure( 0 );
    forAll< POLICY >( perforationSize, [=] GEOSX_HOST_DEVICE ( localIndex const iperf )
    {
      // get the reservoir (sub)region and element indices
      localIndex const er = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei = resElementIndex[iperf];

      minResPressure.min( resPressure[er][esr][ei] );
      maxResPressure.max( resPressure[er][esr][ei] );

      // increment the average total mass density
      for( localIndex ip = 0; ip < NP; ++ip )
      {
        sumTotalMassDens += resPhaseVolFrac[er][esr][ei][ip] * resPhaseMassDens[er][esr][ei][0][ip];
      }
    } );

    // TODO: there must a better way to do what is below
    // I would like to define an array of RAJA::ReduceSum to be able to do sum[ic] += ...
    // and put back what is below in the previous kernel.
    stackArray1d< real64, maxNumComp > sumCompFrac( NC );
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      RAJA::ReduceSum< parallelDeviceReduce, real64 > sum( 0.0 );
      forAll< POLICY >( perforationSize, [=] GEOSX_HOST_DEVICE ( localIndex const iperf )
      {
        // get the reservoir (sub)region and element indices
        localIndex const er = resElementRegion[iperf];
        localIndex const esr = resElementSubRegion[iperf];
        localIndex const ei = resElementIndex[iperf];

        real64 perfTotalDens = 0.0;
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          perfTotalDens += resCompDens[er][esr][ei][jc];
        }
        sum += resCompDens[er][esr][ei][ic] / perfTotalDens;
      } );
      sumCompFrac[ic] = sum.get();
    }

    real64 const pres = ( wellControls.GetType() == WellControls::Type::PRODUCER )
                        ? MpiWrapper::Min( minResPressure.get() )
                        : MpiWrapper::Max( maxResPressure.get() );
    real64 const avgTotalMassDens = MpiWrapper::Sum( sumTotalMassDens.get() ) / numPerforations;

    stackArray1d< real64, maxNumComp > avgCompFrac( NC );
    // compute average component fraction
    if( wellControls.GetType() == WellControls::Type::PRODUCER )
    {
      // use average comp frac from reservoir
      real64 compFracSum = 0;
      real64 const tol = 1e-13;
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        avgCompFrac[ic] = MpiWrapper::Sum( sumCompFrac[ic] ) / numPerforations;
        compFracSum += avgCompFrac[ic];
      }
      GEOSX_ERROR_IF( compFracSum < 1 - tol || compFracSum > 1 + tol,
                      "Invalid well initialization: sum of component fractions should be between 0 and 1" );
    }
    else // injector
    {
      // use average comp frac from XML file
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        avgCompFrac[ic] = wellControls.GetInjectionStream()[ic];
      }
    }

    // set the global component fractions to avgCompFrac
    forAll< POLICY >( subRegionSize, [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        wellElemCompFrac[iwelem][ic] = avgCompFrac[ic];
      }
    } );

    real64 pressureControl = 0.0;
    real64 gravCoefControl = 0.0;
    if( isLocallyOwned )
    {
      // initialize the reference pressure
      if( currentControl == WellControls::Control::BHP )
      {
        // if pressure constraint, set the ref pressure at the constraint
        pressureControl = targetBHP;
      }
      else // rate control
      {
        // if rate constraint, set the ref pressure slightly
        // above/below the target pressure depending on well type
        pressureControl = ( wellType == WellControls::Type::PRODUCER )
                        ? 0.5 * pres
                        : 2.0 * pres;
      }
      gravCoefControl = wellElemGravCoef[iwelemControl];
      wellElemPressure[iwelemControl] = pressureControl;
    }

    MpiWrapper::Broadcast( pressureControl, topRank );
    MpiWrapper::Broadcast( gravCoefControl, topRank );

    GEOSX_ERROR_IF( pressureControl <= 0, "Invalid well initialization: negative pressure was found" );

    // estimate the pressures in the well elements using this avgDens
    forAll< POLICY >( subRegionSize, [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      wellElemPressure[iwelem] = pressureControl
                                 + avgTotalMassDens * ( wellElemGravCoef[iwelem] - gravCoefControl );

    } );
  }

};

/******************************** CompDensInitializationKernel ********************************/

struct CompDensInitializationKernel
{

  template< typename POLICY >
  static void
  Launch( localIndex const subRegionSize,
          localIndex const numComponents,
          arrayView2d< real64 const > const & wellElemCompFrac,
          arrayView2d< real64 const > const & wellElemTotalDens,
          arrayView2d< real64 > const & wellElemCompDens )
  {
    forAll< POLICY >( subRegionSize, [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      for( localIndex ic = 0; ic < numComponents; ++ic )
      {
        wellElemCompDens[iwelem][ic] = wellElemCompFrac[iwelem][ic] * wellElemTotalDens[iwelem][0];
      }
    } );
  }

};

/******************************** TotalMassDensityKernel ****************************/

struct TotalMassDensityKernel
{

  template< typename POLICY >
  static void
  Launch( localIndex const subRegionSize,
          localIndex const numComponents,
          localIndex const numPhases,
          arrayView2d< real64 const > const & phaseVolFrac,
          arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens,
          arrayView3d< real64 const > const & dCompFrac_dCompDens,
          arrayView3d< real64 const > const & phaseMassDens,
          arrayView3d< real64 const > const & dPhaseMassDens_dPres,
          arrayView4d< real64 const > const & dPhaseMassDens_dComp,
          arrayView1d< real64 > const & totalMassDens,
          arrayView1d< real64 > const & dTotalMassDens_dPres,
          arrayView2d< real64 > const & dTotalMassDens_dCompDens )
  {
    localIndex constexpr maxNumComp = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;
    localIndex const NC = numComponents;
    localIndex const NP = numPhases;

    stackArray1d< real64, maxNumComp > dMassDens_dC( NC );

    forAll< POLICY >( subRegionSize, [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      totalMassDens[iwelem] = 0.0;
      dTotalMassDens_dPres[iwelem] = 0.0;
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dTotalMassDens_dCompDens[iwelem][ic] = 0.0;
      }

      for( localIndex ip = 0; ip < NP; ++ip )
      {
        totalMassDens[iwelem] += phaseVolFrac[iwelem][ip] * phaseMassDens[iwelem][0][ip];
        dTotalMassDens_dPres[iwelem] += dPhaseVolFrac_dPres[iwelem][ip] * phaseMassDens[iwelem][0][ip]
                                        + phaseVolFrac[iwelem][ip] * dPhaseMassDens_dPres[iwelem][0][ip];

        applyChainRule( NC, dCompFrac_dCompDens[iwelem], dPhaseMassDens_dComp[iwelem][0][ip], dMassDens_dC );
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dTotalMassDens_dCompDens[iwelem][ic] += dPhaseVolFrac_dCompDens[iwelem][ip][ic] * phaseMassDens[iwelem][0][ip]
                                                  + phaseVolFrac[iwelem][ip] * dMassDens_dC[ic];

        }
      }

    } );
  }

};

/******************************** ResidualNormKernel ********************************/

struct ResidualNormKernel
{

  template< typename POLICY, typename REDUCE_POLICY, typename LOCAL_VECTOR >
  static void
  Launch( LOCAL_VECTOR const localResidual,
          globalIndex const rankOffset,
          localIndex const numComponents,
          localIndex const numDofPerWellElement,
          WellControls const & wellControls,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< integer const > const & wellElemGhostRank,
          arrayView1d< real64 const > const & wellElemVolume,
          arrayView2d< real64 const > const & wellElemTotalDensity,
          real64 const dt,
          real64 * localResidualNorm )
  {
    localIndex const NC = numComponents;

    real64 const targetBHP = wellControls.GetTargetBHP();
    real64 const targetRate = wellControls.GetTargetRate();
    WellControls::Control const currentControl = wellControls.GetControl();
    localIndex const iwelemControl = wellControls.GetReferenceWellElementIndex();

    RAJA::ReduceSum< REDUCE_POLICY, real64 > sumScaled( 0.0 );

    forAll< POLICY >( wellElemDofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      if( wellElemGhostRank[iwelem] < 0 )
      {
        real64 normalizer = 0.0;
        for( localIndex idof = 0; idof < numDofPerWellElement; ++idof )
        {
          if( idof == CompositionalMultiphaseWell::RowOffset::CONTROL )
          {
            if( iwelem == iwelemControl )
            {
              normalizer = ( currentControl == WellControls::Control::BHP )
    ? targetBHP : targetRate;
            }
            else
            {
              normalizer = targetBHP;
            }
          }
          else if( idof >= CompositionalMultiphaseWell::RowOffset::MASSBAL
                   && idof < CompositionalMultiphaseWell::RowOffset::MASSBAL + NC )
          {
            normalizer = dt * targetRate;
          }
          else
          {
            normalizer = wellElemTotalDensity[iwelem][0] * wellElemVolume[iwelem];
          }
          localIndex const lid = wellElemDofNumber[iwelem] + idof - rankOffset;
          real64 const val = localResidual[lid] / normalizer;
          sumScaled += val * val;
        }
      }
    } );
    *localResidualNorm = *localResidualNorm + sumScaled.get();
  }

};

/******************************** SolutionScalingKernel ********************************/

struct SolutionScalingKernel
{
  template< typename POLICY, typename REDUCE_POLICY, typename LOCAL_VECTOR >
  static real64
  Launch( LOCAL_VECTOR const localSolution,
          globalIndex const rankOffset,
          localIndex const numComponents,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< integer const > const & wellElemGhostRank,
          arrayView2d< real64 const > const & wellElemCompDens,
          arrayView2d< real64 const > const & dWellElemCompDens,
          real64 const maxCompFracChange )
  {
    real64 constexpr eps = minDensForDivision;

    RAJA::ReduceMin< REDUCE_POLICY, real64 > minVal( 1.0 );

    forAll< POLICY >( wellElemDofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      if( wellElemGhostRank[iwelem] < 0 )
      {

        real64 prevTotalDens = 0;
        for( localIndex ic = 0; ic < numComponents; ++ic )
        {
          prevTotalDens += wellElemCompDens[iwelem][ic] + dWellElemCompDens[iwelem][ic];
        }

        for( localIndex ic = 0; ic < numComponents; ++ic )
        {
          localIndex const lid = wellElemDofNumber[iwelem] + ic + 1 - rankOffset;

          // compute scaling factor based on relative change in component densities
          real64 const absCompDensChange = fabs( localSolution[lid] );
          real64 const maxAbsCompDensChange = maxCompFracChange * prevTotalDens;

          // This actually checks the change in component fraction, using a lagged total density
          // Indeed we can rewrite the following check as:
          //    | prevCompDens / prevTotalDens - newCompDens / prevTotalDens | > maxCompFracChange
          // Note that the total density in the second term is lagged (i.e, we use prevTotalDens)
          // because I found it more robust than using directly newTotalDens (which can vary also
          // wildly when the compDens change is large)
          if( absCompDensChange > maxAbsCompDensChange && absCompDensChange > eps )
          {
            minVal.min( maxAbsCompDensChange / absCompDensChange );
          }
        }
      }
    } );
    return minVal.get();
  }

};


/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{
  template< typename POLICY, typename REDUCE_POLICY, typename LOCAL_VECTOR >
  static localIndex
  Launch( LOCAL_VECTOR const localSolution,
          globalIndex const rankOffset,
          localIndex const numComponents,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< integer const > const & wellElemGhostRank,
          arrayView1d< real64 const > const & wellElemPressure,
          arrayView1d< real64 const > const & dWellElemPressure,
          arrayView2d< real64 const > const & wellElemCompDens,
          arrayView2d< real64 const > const & dWellElemCompDens,
          integer const allowCompDensChopping,
          real64 const scalingFactor )
  {
    real64 constexpr eps = minDensForDivision;

    RAJA::ReduceMin< REDUCE_POLICY, localIndex > minVal( 1 );

    forAll< POLICY >( wellElemDofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      if( wellElemGhostRank[iwelem] < 0 )
      {
        // pressure
        localIndex lid = wellElemDofNumber[iwelem] + CompositionalMultiphaseWell::ColOffset::DPRES - rankOffset;
        real64 const newPres = wellElemPressure[iwelem] + dWellElemPressure[iwelem]
                               + scalingFactor * localSolution[lid];

        // the pressure must be positive
        if( newPres < 0.0 )
        {
          minVal.min( 0 );
        }

        // if component density is not allowed, the time step fails if a component density is negative
        // otherwise, we just check that the total density is positive, and negative component densities
        // will be chopped (i.e., set to zero) in ApplySystemSolution
        if( !allowCompDensChopping )
        {
          for( localIndex ic = 0; ic < numComponents; ++ic )
          {
            lid = wellElemDofNumber[iwelem] + ic + 1 - rankOffset;
            real64 const newDens = wellElemCompDens[iwelem][ic] + dWellElemCompDens[iwelem][ic]
                                   + scalingFactor * localSolution[lid];

            if( newDens < 0 )
            {
              minVal.min( 0 );
            }
          }
        }
        else
        {
          real64 totalDens = 0.0;
          for( localIndex ic = 0; ic < numComponents; ++ic )
          {
            lid = wellElemDofNumber[iwelem] + ic + 1 - rankOffset;
            real64 const newDens = wellElemCompDens[iwelem][ic] + dWellElemCompDens[iwelem][ic]
                                   + scalingFactor * localSolution[lid];
            totalDens += (newDens > 0.0) ? newDens : 0.0;
          }
          if( totalDens < eps )
          {
            minVal.min( 0 );
          }
        }
      }
    } );
    return minVal.get();
  }

};


} // end namespace CompositionalMultiphaseWellKernels

} // end namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELLKERNELS_HPP
