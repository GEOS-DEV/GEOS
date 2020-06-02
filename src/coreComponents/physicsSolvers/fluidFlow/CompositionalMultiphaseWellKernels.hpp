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
 * @file CompositionalMultiphaseWellKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEWELLKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEWELLKERNELS_HPP

#include "common/DataTypes.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseWell.hpp"
#include "wells/WellControls.hpp"

namespace geosx
{

namespace CompositionalMultiphaseWellKernels
{

/******************************** ControlEquationHelper ********************************/

struct ControlEquationHelper
{

  static void
  ComputeJacobianEntry( globalIndex const rankOffset,
                        localIndex const numComponents,
                        WellControls const & wellControls,
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

    // get well control and type
    WellControls::Control const currentControl = wellControls.GetControl();

    // BHP control
    if( currentControl == WellControls::Control::BHP )
    {

      // get the pressure and compute normalizer
      real64 const currentBHP = wellElemPressure + dWellElemPressure;
      real64 const targetBHP  = wellControls.GetTargetBHP();
      real64 const normalizer = targetBHP > std::numeric_limits< real64 >::epsilon()
                                ? 1.0 / targetBHP
                                : 1.0;

      // control equation is a normalized difference
      // between current pressure and target pressure
      controlEqn = ( currentBHP - targetBHP ) * normalizer;
      dControlEqn_dX = normalizer;
      dofColIndex = wellElemDofNumber + CompositionalMultiphaseWell::ColOffset::DPRES;
      eqnRowIndex = wellElemDofNumber + CompositionalMultiphaseWell::RowOffset::CONTROL - rankOffset;
    }
    else if( currentControl == WellControls::Control::LIQUIDRATE ) // liquid rate control
    {
      // get rates and compute normalizer
      real64 const currentConnRate = connRate + dConnRate;
      real64 const targetConnRate = wellControls.GetTargetRate();
      real64 const normalizer = targetConnRate > std::numeric_limits< real64 >::min()
                                ? 1.0 / ( 1e-2 * targetConnRate ) // hard-coded value comes from AD-GPRS
                                : 1.0;

      // control equation is a normalized difference
      // between current rate and target rate
      controlEqn = ( currentConnRate - targetConnRate ) * normalizer;
      dControlEqn_dX = normalizer;
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
    forAll< POLICY >( size, [=]( localIndex const iwelem )
    {
      localIndex constexpr maxNumComp = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;
      localIndex constexpr maxNumDof  = maxNumComp + 1;

      // create local work arrays
      stackArray1d< real64, maxNumComp > compFracUp( NC );
      stackArray1d< real64, maxNumComp > dCompFrac_dPresUp( NC );
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
        }
        // zero out the derivatives wrt composition
        dCompFrac_dCompDensUp = 0.0;
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

      compFlux = 0.0;
      dCompFlux_dRate = 0.0;
      dCompFlux_dPresUp = 0.0;
      dCompFlux_dCompDensUp = 0.0;

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
          oneSidedEqnRowIndices[ic] = offsetUp + CompositionalMultiphaseWell::RowOffset::MASSBAL + ic;
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
            localMatrix.addToRow< serialAtomic >( oneSidedEqnRowIndices[i],
                                                  &oneSidedDofColIndices_dRate,
                                                  oneSidedFluxJacobian_dRate.data() + i,
                                                  1 );
            localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( oneSidedEqnRowIndices[i],
                                                                      oneSidedDofColIndices_dPresCompUp.data(),
                                                                      oneSidedFluxJacobian_dPresCompUp.data() + i * resNDOF,
                                                                      resNDOF );
            atomicAdd( serialAtomic{}, &localRhs[oneSidedEqnRowIndices[i]], oneSidedFlux[i] );
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
            localMatrix.addToRow< serialAtomic >( eqnRowIndices[i],
                                                  &dofColIndices_dRate,
                                                  localFluxJacobian_dRate.data() + i,
                                                  1 );
            localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowIndices[i],
                                                                      dofColIndices_dPresCompUp.data(),
                                                                      localFluxJacobian_dPresCompUp.data() + i * resNDOF,
                                                                      resNDOF );
            atomicAdd( serialAtomic{}, &localRhs[eqnRowIndices[i]], localFlux[i] );
          }
        }
      }
    } );
  }

};

/******************************** PressureRelationKernel ********************************/

struct PressureRelationKernel
{

  template< typename POLICY >
  static void
  Launch( localIndex const size,
          globalIndex const rankOffset,
          localIndex const numComponents,
          localIndex const numDofPerResElement,
          WellControls const & wellControls,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & wellElemPressure,
          arrayView1d< real64 const > const & dWellElemPressure,
          arrayView1d< real64 const > const & wellElemMixtureDensity,
          arrayView1d< real64 const > const & dWellElemMixtureDensity_dPres,
          arrayView2d< real64 const > const & dWellElemMixtureDensity_dComp,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
  {
    localIndex const NC = numComponents;
    localIndex const resNDOF = numDofPerResElement;

    // compute a coefficient to normalize the momentum equation
    real64 const targetBHP = wellControls.GetTargetBHP();
    real64 const normalizer = targetBHP > std::numeric_limits< real64 >::epsilon()
                              ? 1.0 / targetBHP
                              : 1.0;

    // loop over the well elements to compute the pressure relations between well elements
    forAll< serialPolicy >( size, [=] ( localIndex const iwelem )
    {
      localIndex constexpr maxNumComp = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;
      localIndex constexpr maxNumDof  = maxNumComp + 1;

      localIndex const iwelemNext = nextWellElemIndex[iwelem];

      if( iwelemNext >= 0 )   // if iwelemNext < 0, form control equation, not momentum
      {

        // local working variables and arrays
        stackArray1d< globalIndex, 2 * maxNumDof > dofColIndices( 2 * resNDOF );
        stackArray1d< real64, 2 * maxNumDof > localPresRelJacobian( 2 * resNDOF );

        stackArray1d< real64, maxNumComp > dAvgDensity_dCompCurrent( NC );
        stackArray1d< real64, maxNumComp > dAvgDensity_dCompNext( NC );

        // reset local working variables and arrays
        dofColIndices = -1;
        localPresRelJacobian = 0.0;

        // compute the average density at the interface between well elements
        real64 const avgDensity = 0.5 * ( wellElemMixtureDensity[iwelemNext]
                                          + wellElemMixtureDensity[iwelem] );
        real64 const dAvgDensity_dPresNext = 0.5 * dWellElemMixtureDensity_dPres[iwelemNext];
        real64 const dAvgDensity_dPresCurrent = 0.5 * dWellElemMixtureDensity_dPres[iwelem];
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dAvgDensity_dCompNext[ic] = 0.5 * dWellElemMixtureDensity_dComp[iwelemNext][ic];
          dAvgDensity_dCompCurrent[ic] = 0.5 * dWellElemMixtureDensity_dComp[iwelem][ic];
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

        real64 const localPresRel = ( pressureNext - pressureCurrent - avgDensity * gravD ) * normalizer;

        localPresRelJacobian[localDofIndexPresNext] = ( 1 - dAvgDensity_dPresNext * gravD ) * normalizer;
        localPresRelJacobian[localDofIndexPresCurrent] = ( -1 - dAvgDensity_dPresCurrent * gravD ) * normalizer;

        for( localIndex ic = 0; ic < NC; ++ic )
        {
          localIndex const localDofIndexCompNext = localDofIndexPresNext + ic + 1;
          localIndex const localDofIndexCompCurrent = localDofIndexPresCurrent + ic + 1;

          dofColIndices[localDofIndexCompNext] = offsetNext + CompositionalMultiphaseWell::ColOffset::DCOMP + ic;
          dofColIndices[localDofIndexCompCurrent] = offsetCurrent + CompositionalMultiphaseWell::ColOffset::DCOMP + ic;

          localPresRelJacobian[localDofIndexCompNext] = -dAvgDensity_dCompNext[ic] * gravD * normalizer;
          localPresRelJacobian[localDofIndexCompCurrent] = -dAvgDensity_dCompCurrent[ic] * gravD * normalizer;
        }

        // TODO: add friction and acceleration terms

        if( eqnRowIndex >= 0 && eqnRowIndex < localMatrix.numRows() )
        {
          localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowIndex,
                                                                    dofColIndices.data(),
                                                                    localPresRelJacobian.data(),
                                                                    2 * resNDOF );
          atomicAdd( serialAtomic{}, &localRhs[eqnRowIndex], localPresRel );
        }
      }
    } );
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
  using ElementView = typename ElementRegionManager::ElementViewAccessor< VIEWTYPE >::ViewTypeConst;

  template< typename POLICY >
  static void
  Launch( localIndex const size,
          localIndex const numComponents,
          localIndex const numPhases,
          ElementView< arrayView1d< real64 > > const & resPressure,
          ElementView< arrayView1d< real64 > > const & dResPressure,
          ElementView< arrayView2d< real64 > > const & resPhaseMob,
          ElementView< arrayView2d< real64 > > const & dResPhaseMob_dPres,
          ElementView< arrayView3d< real64 > > const & dResPhaseMob_dComp,
          ElementView< arrayView2d< real64 > > const & dResPhaseVolFrac_dPres,
          ElementView< arrayView3d< real64 > > const & dResPhaseVolFrac_dComp,
          ElementView< arrayView3d< real64 > > const & dResCompFrac_dCompDens,
          ElementView< arrayView3d< real64 > > const & resPhaseVisc,
          ElementView< arrayView3d< real64 > > const & dResPhaseVisc_dPres,
          ElementView< arrayView4d< real64 > > const & dResPhaseVisc_dComp,
          ElementView< arrayView4d< real64 > > const & resPhaseCompFrac,
          ElementView< arrayView4d< real64 > > const & dResPhaseCompFrac_dPres,
          ElementView< arrayView5d< real64 > > const & dResPhaseCompFrac_dComp,
          ElementView< arrayView3d< real64 > > const & resPhaseRelPerm,
          ElementView< arrayView4d< real64 > > const & dResPhaseRelPerm_dPhaseVolFrac,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< real64 const > const & wellElemPressure,
          arrayView1d< real64 const > const & dWellElemPressure,
          arrayView1d< real64 const > const & wellElemMixtureDensity,
          arrayView1d< real64 const > const & dWellElemMixtureDensity_dPres,
          arrayView2d< real64 const > const & dWellElemMixtureDensity_dComp,
          arrayView2d< real64 const > const & wellElemCompFrac,
          arrayView3d< real64 const > const & dWellElemCompFrac_dCompDens,
          arrayView1d< real64 const > const & perfGravCoef,
          arrayView1d< localIndex const > const & perfWellElemIndex,
          arrayView1d< real64 const > const & perfTransmissibility,
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
    forAll< serialPolicy >( size, [=] ( localIndex const iperf )
    {
      localIndex constexpr maxNumComp = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;

      // local working variables and arrays
      stackArray1d< real64, maxNumComp > dPhaseCompFrac_dCompDens( NC );

      stackArray1d< real64, 2 > pressure( 2 );
      stackArray1d< real64, 2 > dPressure_dP( 2 );
      stackArray2d< real64, 2 * maxNumComp > dPressure_dC( 2, NC );

      stackArray1d< real64, 2 > dFlux_dP( 2 );
      stackArray2d< real64, 2 * maxNumComp > dFlux_dC( 2, NC );

      stackArray1d< real64, 2 > dMult_dP( 2 );
      stackArray2d< real64, 2 * maxNumComp > dMult_dC( 2, NC );

      stackArray1d< real64, maxNumComp > dMixtureDensity_dC( NC );
      stackArray1d< real64, maxNumComp > dResTotalMobility_dC( NC );

      stackArray2d< real64, 2 * maxNumComp > phaseCompFrac( 2, NC );
      stackArray2d< real64, 2 * maxNumComp > dPhaseCompFrac_dP( 2, NC );
      stackArray3d< real64, 2 * maxNumComp * maxNumComp > dPhaseCompFrac_dC( 2, NC, NC );

      stackArray1d< real64, maxNumComp > dVisc_dC( NC );
      stackArray1d< real64, maxNumComp > dRelPerm_dC( NC );

      stackArray1d< real64, 2 > dPotDiff_dP( 2 );
      stackArray2d< real64, 2 * maxNumComp > dPotDiff_dC( 2, NC );

      stackArray1d< real64, 2 > multiplier( 2 );

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

      // clear working arrays
      pressure = 0.0;
      dPressure_dP = 0.0;
      dPressure_dC = 0.0;

      // 1) copy the variables from the reservoir and well element

      // a) get reservoir variables

      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];
      // get the index of the well elem
      localIndex const iwelem = perfWellElemIndex[iperf];

      pressure[CompositionalMultiphaseWell::SubRegionTag::RES] = resPressure[er][esr][ei] + dResPressure[er][esr][ei];
      dPressure_dP[CompositionalMultiphaseWell::SubRegionTag::RES] = 1.0;

      // TODO: add a buoyancy term for the reservoir side here

      multiplier[CompositionalMultiphaseWell::SubRegionTag::RES] = 1.0;

      // b) get well variables

      pressure[CompositionalMultiphaseWell::SubRegionTag::WELL] = wellElemPressure[iwelem] + dWellElemPressure[iwelem];
      dPressure_dP[CompositionalMultiphaseWell::SubRegionTag::WELL] = 1.0;

      multiplier[CompositionalMultiphaseWell::SubRegionTag::WELL] = -1.0;

      real64 const gravD = ( perfGravCoef[iperf] - wellElemGravCoef[iwelem] );

      pressure[CompositionalMultiphaseWell::SubRegionTag::WELL] += wellElemMixtureDensity[iwelem] * gravD;
      dPressure_dP[CompositionalMultiphaseWell::SubRegionTag::WELL] += dWellElemMixtureDensity_dPres[iwelem] * gravD;
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dPressure_dC[CompositionalMultiphaseWell::SubRegionTag::WELL][ic] +=
          dWellElemMixtureDensity_dComp[iwelem][ic] * gravD;
      }

      // get transmissibility at the interface
      real64 const trans = perfTransmissibility[iperf];

      // 2) compute potential difference

      real64 potDiff = 0.0;
      dPotDiff_dP = 0.0;
      dPotDiff_dC = 0.0;

      for( localIndex i = 0; i < 2; ++i )
      {
        potDiff += multiplier[i] * trans * pressure[i]; // pressure = pres + dPres
        dPotDiff_dP[i] += multiplier[i] * trans * dPressure_dP[i];

        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dPotDiff_dC[i][ic] += multiplier[i] * trans * dPressure_dC[i][ic];
        }
      }

      real64 flux = 0.0;
      dFlux_dP = 0.0;
      dFlux_dC = 0.0;

      // 3) upwinding

      if( potDiff >= 0 )  // ** reservoir cell is upstream **
      {

        phaseCompFrac = 0.0;
        dPhaseCompFrac_dP = 0.0;
        dPhaseCompFrac_dC = 0.0;

        dPhaseCompFrac_dCompDens = 0;

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
              resPhaseMob[er][esr][ei][ip] * dPotDiff_dC[CompositionalMultiphaseWell::SubRegionTag::RES][ic];
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

        real64 resTotalMobility     = 0.0;
        real64 dResTotalMobility_dP = 0.0;
        dResTotalMobility_dC        = 0.0;

        // first, compute the reservoir total mobitity (excluding phase density)
        for( localIndex ip = 0; ip < NP; ++ip )
        {
          // viscosity
          real64 const resViscosity = resPhaseVisc[er][esr][ei][0][ip];
          real64 const dResVisc_dP  = dResPhaseVisc_dPres[er][esr][ei][0][ip];
          applyChainRule( NC, dResCompFrac_dCompDens[er][esr][ei],
                          dResPhaseVisc_dComp[er][esr][ei][0][ip],
                          dVisc_dC );

          // relative permeability
          real64 const resRelPerm = resPhaseRelPerm[er][esr][ei][0][ip];
          real64 dResRelPerm_dP = 0.0;
          dRelPerm_dC = 0.0;
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
          resTotalMobility     += resRelPerm / resViscosity;
          dResTotalMobility_dP += ( dResRelPerm_dP * resViscosity - resRelPerm * dResVisc_dP )
                                  / ( resViscosity * resViscosity );
          for( localIndex ic = 0; ic < NC; ++ic )
          {
            dResTotalMobility_dC[ic] += ( dRelPerm_dC[ic] * resViscosity - resRelPerm * dVisc_dC[ic] )
                                        / ( resViscosity * resViscosity );
          }
        }

        // compute a potdiff multiplier = wellElemMixtureDensity * resTotalMobility
        real64 const mult = wellElemMixtureDensity[iwelem] * resTotalMobility;
        dMult_dP[CompositionalMultiphaseWell::SubRegionTag::RES]  =
          wellElemMixtureDensity[iwelem] * dResTotalMobility_dP;
        dMult_dP[CompositionalMultiphaseWell::SubRegionTag::WELL] =
          dWellElemMixtureDensity_dPres[iwelem] * resTotalMobility;

        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dMult_dC[CompositionalMultiphaseWell::SubRegionTag::RES][ic]  =
            wellElemMixtureDensity[iwelem] * dResTotalMobility_dC[ic];

          dMult_dC[CompositionalMultiphaseWell::SubRegionTag::WELL][ic] =
            dWellElemMixtureDensity_dComp[iwelem][ic] * resTotalMobility;
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
            + mult * dPotDiff_dC[CompositionalMultiphaseWell::SubRegionTag::RES][ic];
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
            dCompPerfRate_dComp[iperf][CompositionalMultiphaseWell::SubRegionTag::RES][ic][jc]  +=
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

    forAll< serialPolicy >( size, [=] ( localIndex const iwelem )
    {

      if( wellElemGhostRank[iwelem] >= 0 )
      {
        return;
      }

      stackArray1d< globalIndex, maxNumDof > localVolBalanceDOF( welemNDOF );
      stackArray1d< real64, maxNumDof > localVolBalanceJacobian( welemNDOF );

      localVolBalanceDOF = 0;
      localVolBalanceJacobian = 0;

      // get equation/dof indices
      globalIndex const offset = wellElemDofNumber[iwelem];
      localIndex const volBalRowOffset = CompositionalMultiphaseWell::RowOffset::MASSBAL + NC;
      globalIndex const localVolBalanceEqnIndex = offset + volBalRowOffset;
      for( localIndex jdof = 0; jdof < welemNDOF; ++jdof )
      {
        localVolBalanceDOF[jdof] = offset + CompositionalMultiphaseWell::ColOffset::DPRES + jdof;
      }

      real64 localVolBalance = 1.0;
      localVolBalanceJacobian = 0.0;

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

/******************************** MixtureDensityKernel ********************************/

struct MixtureDensityKernel
{

  template< typename POLICY >
  static void
  Launch( localIndex const size,
          localIndex const numComponents,
          localIndex const numPhases,
          arrayView2d< real64 const > const & wellElemPhaseVolFrac,
          arrayView2d< real64 const > const & dWellElemPhaseVolFrac_dPres,
          arrayView3d< real64 const > const & dWellElemPhaseVolFrac_dComp,
          arrayView3d< real64 const > const & dWellElemCompFrac_dCompDens,
          arrayView3d< real64 const > const & wellElemPhaseDens,
          arrayView3d< real64 const > const & dWellElemPhaseDens_dPres,
          arrayView4d< real64 const > const & dWellElemPhaseDens_dComp,
          arrayView1d< real64 > const & wellElemMixtureDensity,
          arrayView1d< real64 > const & dWellElemMixtureDensity_dPres,
          arrayView2d< real64 > const & dWellElemMixtureDensity_dComp )
  {
    localIndex constexpr maxNumComp = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;
    localIndex const NC = numComponents;
    localIndex const NP = numPhases;

    forAll< serialPolicy >( size, [=] ( localIndex const iwelem )
    {
      stackArray1d< real64, maxNumComp > dDens_dC( NC );

      // reset to zero
      wellElemMixtureDensity[iwelem]        = 0.0;
      dWellElemMixtureDensity_dPres[iwelem] = 0.0;
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dWellElemMixtureDensity_dComp[iwelem][ic] = 0.0;
      }

      // increment mixture velocity
      for( localIndex ip = 0; ip < NP; ++ip )
      {
        wellElemMixtureDensity[iwelem] += wellElemPhaseVolFrac[iwelem][ip] * wellElemPhaseDens[iwelem][0][ip];
        dWellElemMixtureDensity_dPres[iwelem] += dWellElemPhaseVolFrac_dPres[iwelem][ip] * wellElemPhaseDens[iwelem][0][ip]
                                                 + wellElemPhaseVolFrac[iwelem][ip] * dWellElemPhaseDens_dPres[iwelem][0][ip];

        dDens_dC = 0.0;

        applyChainRule( NC, dWellElemCompFrac_dCompDens[iwelem], dWellElemPhaseDens_dComp[iwelem][0][ip], dDens_dC );
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dWellElemMixtureDensity_dComp[iwelem][ic] += dWellElemPhaseVolFrac_dComp[iwelem][ip][ic] * wellElemPhaseDens[iwelem][0][ip]
                                                       + wellElemPhaseVolFrac[iwelem][ip] * dDens_dC[ic];
        }
      }
    } );
  }

};

/******************************** ResidualNormKernel ********************************/

struct ResidualNormKernel
{

  template< typename POLICY, typename REDUCE_POLICY, typename LOCAL_VECTOR >
  static void Launch( LOCAL_VECTOR const localResidual,
                      globalIndex const rankOffset,
                      localIndex const numComponents,
                      localIndex const numDofPerWellElement,
                      arrayView1d< globalIndex const > const & wellElemDofNumber,
                      arrayView1d< integer const > const & wellElemGhostRank,
                      arrayView1d< real64 const > const & wellElemVolume,
                      arrayView2d< real64 const > const & wellElemTotalDensity,
                      real64 * localResidualNorm )
  {
    localIndex const NC = numComponents;

    RAJA::ReduceSum< REDUCE_POLICY, real64 > sumScaled( 0.0 );

    forAll< POLICY >( wellElemDofNumber.size(), [=] ( localIndex const iwelem )
    {
      if( wellElemGhostRank[iwelem] < 0 )
      {
        for( localIndex idof = 0; idof < numDofPerWellElement; ++idof )
        {
          real64 const normalizer = ( idof >= CompositionalMultiphaseWell::RowOffset::MASSBAL
                                      && idof < CompositionalMultiphaseWell::RowOffset::MASSBAL + NC )
                                    ? wellElemTotalDensity[iwelem][0] * wellElemVolume[iwelem]
                                    : 1;
          localIndex const lid = wellElemDofNumber[iwelem] + idof - rankOffset;
          real64 const val = localResidual[lid] / normalizer;
          sumScaled += val * val;
        }
      }
    } );
    *localResidualNorm = *localResidualNorm + sumScaled.get();
  }

};

/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{
  template< typename POLICY, typename REDUCE_POLICY, typename LOCAL_VECTOR >
  static localIndex Launch( LOCAL_VECTOR const localSolution,
                            globalIndex const rankOffset,
                            localIndex const numComponents,
                            arrayView1d< globalIndex const > const & wellElemDofNumber,
                            arrayView1d< integer const > const & wellElemGhostRank,
                            arrayView1d< real64 const > const & wellElemPressure,
                            arrayView1d< real64 const > const & dWellElemPressure,
                            arrayView2d< real64 const > const & wellElemCompDens,
                            arrayView2d< real64 const > const & dWellElemCompDens,
                            real64 const scalingFactor )
  {
    RAJA::ReduceMin< REDUCE_POLICY, localIndex > minVal( 1 );

    forAll< POLICY >( wellElemDofNumber.size(), [=] ( localIndex const iwelem )
    {
      if( wellElemGhostRank[iwelem] < 0 )
      {
        // pressure
        localIndex lid = wellElemDofNumber[iwelem] + CompositionalMultiphaseWell::ColOffset::DPRES - rankOffset;
        real64 const newPres = wellElemPressure[iwelem] + dWellElemPressure[iwelem]
                               + scalingFactor * localSolution[lid];
        if( newPres < 0.0 )
        {
          minVal.min( 0 );
        }

        // comp densities
        for( localIndex ic = 0; ic < numComponents; ++ic )
        {
          lid = wellElemDofNumber[iwelem] + ic + 1 - rankOffset;
          real64 const newDens = wellElemCompDens[iwelem][ic] + dWellElemCompDens[iwelem][ic]
                                 + scalingFactor * localSolution[lid];
          if( newDens < 0.0 )
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

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEWELLKERNELS_HPP
