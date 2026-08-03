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
 * @file CompositionalMultiphaseWellConstraintKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTRAINTKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTRAINTKERNELS_HPP

#include "codingUtilities/Utilities.hpp"
#include "common/DataLayouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"


#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/fluidFlow/wells/WellBHPConstraints.hpp"
#include "physicsSolvers/fluidFlow/wells/WellVolumeRateConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellPhaseVolumeRateConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellMassRateConstraint.hpp"
namespace geos
{

namespace wellConstraintKernels
{

/******************************** ControlEquationHelper ********************************/
//template< integer NC, integer IS_THERMAL, typname S, typename T  >
//struct ConstraintHelper< NC, IS_THERMAL   > {};

template< integer NC, integer IS_THERMAL, typename CONSTRAINT = BHPConstraint< BHPConstraintTypeId::MIN > >
struct ConstraintHelper
{
  template< BHPConstraintTypeId I >
  static void assembleConstraintEquation( real64 const & time_n,
                                          WellControls & wellControls,
                                          BHPConstraint< I > & constraint,
                                          WellElementSubRegion const & subRegion,
                                          string const & wellDofKey,
                                          localIndex const & rankOffset,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
  {
    // subRegion data
    localIndex const iwelemRef = subRegion.getTopWellElementIndex();
    arrayView1d< globalIndex const > const & wellElemDofNumber = subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< real64 const > const & pres = subRegion.getField< fields::well::pressure >();
    arrayView1d< real64 const > const & totalMassDens = subRegion.getField< fields::well::totalMassDensity >();
    arrayView2d< real64 const, constitutive::multifluid::USD_FLUID > const & dTotalMassDens = subRegion.getField< fields::well::dTotalMassDensity >();
    arrayView1d< real64 const > const wellElemGravCoef = subRegion.getField< fields::well::gravityCoefficient >();

    // setup row/column indices for constraint equation
    using COFFSET_WJ = compositionalMultiphaseWellKernels::ColOffset_WellJac< NC, IS_THERMAL >;
    using WJ_ROFFSET = compositionalMultiphaseWellKernels::RowOffset_WellJac< NC, IS_THERMAL >;
    using Deriv = constitutive::multifluid::DerivativeOffset;

    localIndex const eqnRowIndex      = wellElemDofNumber[iwelemRef] + WJ_ROFFSET::CONTROL - rankOffset;
    globalIndex dofColIndices[COFFSET_WJ::nDer]{};
    for( integer ic = 0; ic < COFFSET_WJ::nDer; ++ic )
    {
      dofColIndices[ ic ] =  wellElemDofNumber[iwelemRef] + ic;
    }

    // constraint data
    real64 const & targetBHP = constraint.getConstraintValue( time_n );
    real64 const & refGravCoef = constraint.getReferenceGravityCoef();

    // current constraint value
    real64 const & currentBHP =
      wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() );

    // residual
    real64 controlEqn = currentBHP - targetBHP;

    // setup Jacobian terms
    real64 dControlEqn[NC+2+IS_THERMAL]{};

    // bring everything back to host, capture the scalars by reference
    forAll< serialPolicy >( 1, [pres,
                                totalMassDens,
                                dTotalMassDens,
                                wellElemGravCoef,
                                &dControlEqn,
                                &iwelemRef,
                                localRhs,
                                controlEqn,
                                eqnRowIndex,
                                dofColIndices,
                                localMatrix,
                                &refGravCoef] ( localIndex const )
    {
      real64 const diffGravCoef = refGravCoef - wellElemGravCoef[iwelemRef];
      dControlEqn[COFFSET_WJ::dP] =   1 + dTotalMassDens[iwelemRef][Deriv::dP] * diffGravCoef;
      for( integer ic = 0; ic < NC; ++ic )
      {
        dControlEqn[COFFSET_WJ::dC+ic] = dTotalMassDens[iwelemRef][Deriv::dC+ic] * diffGravCoef;
      }
      if constexpr ( IS_THERMAL )
      {
        dControlEqn[COFFSET_WJ::dT] =  dTotalMassDens[iwelemRef][Deriv::dT] * diffGravCoef;
      }
      // add solver matrices
      localRhs[eqnRowIndex] += controlEqn;
      localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowIndex,
                                                                dofColIndices,
                                                                dControlEqn,
                                                                COFFSET_WJ::nDer );
    } );


  }


  template< template< typename U > class T, typename U=PhaseVolumeRateConstraint >
  static void assembleConstraintEquation( real64 const & time_n,
                                          WellControls & wellControls,
                                          T< PhaseVolumeRateConstraint > & constraint,
                                          WellElementSubRegion const & subRegion,
                                          string const & wellDofKey,
                                          localIndex const & rankOffset,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
  {
    // subRegion data

    localIndex const iwelemRef = subRegion.getTopWellElementIndex();
    arrayView1d< globalIndex const > const & wellElemDofNumber = subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< real64 const > const & connRate = subRegion.getField< fields::well::connectionRate >();
    arrayView2d< real64 const, compflow::USD_COMP > const & compFrac = subRegion.getField< fields::well::globalCompFraction >();
    arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens = subRegion.getField< fields::well::dGlobalCompFraction_dGlobalCompDensity >();

    // setup row/column indices for constraint equation
    using COFFSET_WJ = compositionalMultiphaseWellKernels::ColOffset_WellJac< NC, IS_THERMAL >;
    using WJ_ROFFSET = compositionalMultiphaseWellKernels::RowOffset_WellJac< NC, IS_THERMAL >;
    using Deriv = constitutive::multifluid::DerivativeOffset;

    localIndex const eqnRowIndex      = wellElemDofNumber[iwelemRef] + WJ_ROFFSET::CONTROL - rankOffset;
    globalIndex dofColIndices[COFFSET_WJ::nDer]{};
    for( integer ic = 0; ic < COFFSET_WJ::nDer; ++ic )
    {
      dofColIndices[ ic ] = wellElemDofNumber[iwelemRef] + ic;
    }

    // fluid data
    constitutive::MultiFluidBase & fluidSeparator =  wellControls.getMultiFluidSeparator();

    arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseFrac = fluidSeparator.phaseFraction();
    arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const & dPhaseFrac = fluidSeparator.dPhaseFraction();
    arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseDens = fluidSeparator.phaseDensity();
    arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const & dPhaseDens = fluidSeparator.dPhaseDensity();

    // constraint data
    integer ip = getPhaseIndexFromFluidModel( fluidSeparator, constraint.getPhaseName());
    real64 const & targetPhaseRate = constraint.getConstraintValue( time_n );

    // current constraint value
    arrayView1d< real64 > const & currentPhaseVolRate =
      wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() );
    integer const useSurfaceConditions = wellControls.useSurfaceConditions();

    // residual
    real64 controlEqn =  currentPhaseVolRate[ip] - targetPhaseRate;

    // setup Jacobian terms
    real64 dControlEqn[NC+2+IS_THERMAL]{};

    // bring everything back to host, capture the scalars by reference
    forAll< serialPolicy >( 1, [&ip,
                                connRate,
                                phaseDens,
                                dPhaseDens,
                                phaseFrac,
                                dPhaseFrac,
                                compFrac,
                                dCompFrac_dCompDens,
                                &dControlEqn,
                                &useSurfaceConditions,
                                localRhs,
                                controlEqn,
                                eqnRowIndex,
                                dofColIndices,
                                localMatrix,
                                &iwelemRef] ( localIndex const )
    {
      // skip the rest of this function if phase ip is absent
      bool const phaseExists = (phaseFrac[iwelemRef][0][ip] > 0);
      if( phaseExists )
      {
        stackArray1d< real64, NC > work( NC );
        real64 const currentTotalRate = connRate[iwelemRef];

        real64 const phaseDensInv =  1.0 / phaseDens[iwelemRef][0][ip];
        real64 const phaseFracTimesPhaseDensInv = phaseFrac[iwelemRef][0][ip] * phaseDensInv;
        real64 const dPhaseFracTimesPhaseDensInv_dPres = dPhaseFrac[iwelemRef][0][ip][Deriv::dP] * phaseDensInv
                                                         - dPhaseDens[iwelemRef][0][ip][Deriv::dP] * phaseFracTimesPhaseDensInv * phaseDensInv;

        // divide the total mass/molar rate by the (phase density * phase fraction) to get the phase volumetric rate
        dControlEqn[COFFSET_WJ::dP] = ( useSurfaceConditions ==  0 ) * currentTotalRate * dPhaseFracTimesPhaseDensInv_dPres;
        dControlEqn[COFFSET_WJ::dQ] = phaseFracTimesPhaseDensInv;
        if constexpr (IS_THERMAL )
        {
          real64 const dPhaseFracTimesPhaseDensInv_dTemp = dPhaseFrac[iwelemRef][0][ip][Deriv::dT] * phaseDensInv
                                                           - dPhaseDens[iwelemRef][0][ip][Deriv::dT] * phaseFracTimesPhaseDensInv * phaseDensInv;
          dControlEqn[COFFSET_WJ::dT] = ( useSurfaceConditions ==  0 ) * currentTotalRate * dPhaseFracTimesPhaseDensInv_dTemp;
        }

        for( integer ic = 0; ic < NC; ++ic )
        {
          dControlEqn[COFFSET_WJ::dC+ic] = -phaseFracTimesPhaseDensInv * dPhaseDens[iwelemRef][0][ip][Deriv::dC+ic] * phaseDensInv;
          dControlEqn[COFFSET_WJ::dC+ic] += dPhaseFrac[iwelemRef][0][ip][Deriv::dC+ic] * phaseDensInv;
          dControlEqn[COFFSET_WJ::dC+ic] *= currentTotalRate;
        }
        applyChainRuleInPlace( NC, dCompFrac_dCompDens[iwelemRef], &dControlEqn[COFFSET_WJ::dC], work.data() );
        // add solver matrices
        localRhs[eqnRowIndex] += controlEqn;
        localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowIndex,
                                                                  dofColIndices,
                                                                  dControlEqn,
                                                                  COFFSET_WJ::nDer );
      }
    } );


  }

  template< template< typename U > class T, typename U=LiquidRateConstraint >
  static void assembleConstraintEquation( real64 const & time_n,
                                          WellControls & wellControls,
                                          T< LiquidRateConstraint > & constraint,
                                          WellElementSubRegion const & subRegion,
                                          string const & wellDofKey,
                                          localIndex const & rankOffset,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
  {
    // subRegion data

    localIndex const iwelemRef = subRegion.getTopWellElementIndex();
    arrayView1d< globalIndex const > const & wellElemDofNumber = subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< real64 const > const & connRate = subRegion.getField< fields::well::connectionRate >();
    arrayView2d< real64 const, compflow::USD_COMP > const & compFrac = subRegion.getField< fields::well::globalCompFraction >();
    arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens = subRegion.getField< fields::well::dGlobalCompFraction_dGlobalCompDensity >();

    // setup row/column indices for constraint equation
    using COFFSET_WJ = compositionalMultiphaseWellKernels::ColOffset_WellJac< NC, IS_THERMAL >;
    using WJ_ROFFSET = compositionalMultiphaseWellKernels::RowOffset_WellJac< NC, IS_THERMAL >;
    using Deriv = constitutive::multifluid::DerivativeOffset;

    localIndex const eqnRowIndex      = wellElemDofNumber[iwelemRef] + WJ_ROFFSET::CONTROL - rankOffset;
    globalIndex dofColIndices[COFFSET_WJ::nDer]{};
    for( integer ic = 0; ic < COFFSET_WJ::nDer; ++ic )
    {
      dofColIndices[ ic ] = wellElemDofNumber[iwelemRef] + ic;
    }

    // fluid data
    constitutive::MultiFluidBase & fluidSeparator =  wellControls.getMultiFluidSeparator();

    arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseFrac = fluidSeparator.phaseFraction();
    arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const & dPhaseFrac = fluidSeparator.dPhaseFraction();
    arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseDens = fluidSeparator.phaseDensity();
    arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const & dPhaseDens = fluidSeparator.dPhaseDensity();

    // constraint data
    real64 const & targetPhaseRate = constraint.getConstraintValue( time_n );
    const array1d< int > & phaseIndices = constraint.getPhaseIndices();


    // current constraint value
    arrayView1d< real64 > const & currentPhaseVolRate =
      wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() );
    integer const useSurfaceConditions = wellControls.useSurfaceConditions();

    // residual
    real64 controlEqn =  -targetPhaseRate;
    for( integer ip= 0; ip < phaseIndices.size(); ++ip )
    {
      controlEqn += currentPhaseVolRate[phaseIndices[ip]];
      std::cout << "Phase " << phaseIndices[ip] << " currentPhaseVolRate: " << currentPhaseVolRate[phaseIndices[ip]] <<  " " << controlEqn << " " << targetPhaseRate << std::endl;
    }

    // setup Jacobian terms
    real64 dControlEqn[NC+2+IS_THERMAL]{};

    // bring everything back to host, capture the scalars by reference
    forAll< serialPolicy >( 1, [phaseIndices,
                                connRate,
                                phaseDens,
                                dPhaseDens,
                                phaseFrac,
                                dPhaseFrac,
                                compFrac,
                                dCompFrac_dCompDens,
                                &dControlEqn,
                                &useSurfaceConditions,
                                localRhs,
                                controlEqn,
                                eqnRowIndex,
                                dofColIndices,
                                localMatrix,
                                &iwelemRef] ( localIndex const )
    {

      stackArray1d< real64, NC > work( NC );
      for( integer i= 0; i < phaseIndices.size(); ++i )
      {
        integer ip = phaseIndices[i];
        bool const phaseExists = (phaseFrac[iwelemRef][0][ip] > 0);
        // skip the rest of this function if phase ip is absent
        if( phaseExists )
        {

          real64 const currentTotalRate = connRate[iwelemRef];

          real64 const phaseDensInv =  1.0 / phaseDens[iwelemRef][0][ip];
          real64 const phaseFracTimesPhaseDensInv = phaseFrac[iwelemRef][0][ip] * phaseDensInv;
          real64 const dPhaseFracTimesPhaseDensInv_dPres = dPhaseFrac[iwelemRef][0][ip][Deriv::dP] * phaseDensInv
                                                           - dPhaseDens[iwelemRef][0][ip][Deriv::dP] * phaseFracTimesPhaseDensInv * phaseDensInv;

          // divide the total mass/molar rate by the (phase density * phase fraction) to get the phase volumetric rate
          dControlEqn[COFFSET_WJ::dP] += ( useSurfaceConditions ==  0 ) * currentTotalRate * dPhaseFracTimesPhaseDensInv_dPres;
          dControlEqn[COFFSET_WJ::dQ] += phaseFracTimesPhaseDensInv;
          if constexpr (IS_THERMAL )
          {
            real64 const dPhaseFracTimesPhaseDensInv_dTemp = dPhaseFrac[iwelemRef][0][ip][Deriv::dT] * phaseDensInv
                                                             - dPhaseDens[iwelemRef][0][ip][Deriv::dT] * phaseFracTimesPhaseDensInv * phaseDensInv;
            dControlEqn[COFFSET_WJ::dT] += ( useSurfaceConditions ==  0 ) * currentTotalRate * dPhaseFracTimesPhaseDensInv_dTemp;
          }

          for( integer ic = 0; ic < NC; ++ic )
          {
            real64 temp = -phaseFracTimesPhaseDensInv * dPhaseDens[iwelemRef][0][ip][Deriv::dC+ic] * phaseDensInv;
            temp += dPhaseFrac[iwelemRef][0][ip][Deriv::dC+ic] * phaseDensInv;
            temp *= currentTotalRate;
            dControlEqn[COFFSET_WJ::dC+ic] += temp;
          }
        }
      }
      applyChainRuleInPlace( NC, dCompFrac_dCompDens[iwelemRef], &dControlEqn[COFFSET_WJ::dC], work.data() );
      // add solver matrices
      localRhs[eqnRowIndex] += controlEqn;
      localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowIndex,
                                                                dofColIndices,
                                                                dControlEqn,
                                                                COFFSET_WJ::nDer );
    } );


  }
  template< template< typename U > class T, typename U=VolumeRateConstraint >
  static void assembleConstraintEquation( real64 const & time_n,
                                          WellControls & wellControls,
                                          T< VolumeRateConstraint > & constraint,
                                          WellElementSubRegion const & subRegion,
                                          string const & wellDofKey,
                                          localIndex const & rankOffset,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
  {
    // subRegion data

    localIndex const iwelemRef = subRegion.getTopWellElementIndex();
    arrayView1d< globalIndex const > const & wellElemDofNumber = subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< real64 const > const & connRate = subRegion.getField< fields::well::connectionRate >();
    arrayView2d< real64 const, compflow::USD_COMP > const & compFrac = subRegion.getField< fields::well::globalCompFraction >();
    arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens = subRegion.getField< fields::well::dGlobalCompFraction_dGlobalCompDensity >();

    // setup row/column indices for constraint equation
    using COFFSET_WJ = compositionalMultiphaseWellKernels::ColOffset_WellJac< NC, IS_THERMAL >;
    using WJ_ROFFSET = compositionalMultiphaseWellKernels::RowOffset_WellJac< NC, IS_THERMAL >;
    using Deriv = constitutive::multifluid::DerivativeOffset;

    localIndex const eqnRowIndex      = wellElemDofNumber[iwelemRef] + WJ_ROFFSET::CONTROL - rankOffset;
    globalIndex dofColIndices[COFFSET_WJ::nDer]{};
    for( integer ic = 0; ic < COFFSET_WJ::nDer; ++ic )
    {
      dofColIndices[ ic ] = wellElemDofNumber[iwelemRef] + ic;
    }

    // fluid data
    constitutive::MultiFluidBase & fluidSeparator =  wellControls.getMultiFluidSeparator();

    arrayView2d< real64 const, constitutive::multifluid::USD_FLUID > const & totalDens = fluidSeparator.totalDensity();
    arrayView3d< real64 const, constitutive::multifluid::USD_FLUID_DC > const & dTotalDens = fluidSeparator.dTotalDensity();

    // constraint data
    real64 const & targetTotalVolRate = constraint.getConstraintValue( time_n );

    // current constraint value
    real64 const & currentTotalVolRate =
      wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() );
    integer const useSurfaceConditions = wellControls.useSurfaceConditions();

    // residual
    real64 controlEqn =  currentTotalVolRate - targetTotalVolRate;

    // setup Jacobian terms
    real64 dControlEqn[NC+2+IS_THERMAL]{};

    // bring everything back to host, capture the scalars by reference
    forAll< serialPolicy >( 1, [connRate,
                                totalDens,
                                dTotalDens,
                                compFrac,
                                dCompFrac_dCompDens,
                                &dControlEqn,
                                &useSurfaceConditions,
                                localRhs,
                                controlEqn,
                                eqnRowIndex,
                                dofColIndices,
                                localMatrix,
                                &iwelemRef] ( localIndex const )
    {
      stackArray1d< real64, NC > work( NC );

      real64 const currentTotalRate = connRate[iwelemRef];

      // compute the inverse of the total density and derivatives

      real64 const totalDensInv = 1.0 / totalDens[iwelemRef][0];

      stackArray1d< real64, NC > dTotalDensInv_dCompDens( NC );
      for( integer ic = 0; ic < NC; ++ic )
      {
        dTotalDensInv_dCompDens[ic] = -dTotalDens[iwelemRef][0][Deriv::dC+ic] * totalDensInv * totalDensInv;
      }
      applyChainRuleInPlace( NC, dCompFrac_dCompDens[iwelemRef], dTotalDensInv_dCompDens, work.data() );

      // Step 2.2: divide the total mass/molar rate by the total density to get the total volumetric rate

      // Compute derivatives  dP dT
      real64 const dTotalDensInv_dPres = -dTotalDens[iwelemRef][0][Deriv::dP] * totalDensInv * totalDensInv;
      dControlEqn[COFFSET_WJ::dP] = ( useSurfaceConditions ==  0 ) * currentTotalRate * dTotalDensInv_dPres;
      if constexpr ( IS_THERMAL )
      {
        dControlEqn[COFFSET_WJ::dT] = ( useSurfaceConditions ==  0 ) * currentTotalRate * -dTotalDens[iwelemRef][0][Deriv::dT] * totalDensInv * totalDensInv;
      }

      dControlEqn[COFFSET_WJ::dQ] = totalDensInv;
      for( integer ic = 0; ic < NC; ++ic )
      {
        dControlEqn[COFFSET_WJ::dC+ic] = currentTotalRate * dTotalDensInv_dCompDens[ic];
      }
      localRhs[eqnRowIndex] += controlEqn;
      localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowIndex,
                                                                dofColIndices,
                                                                dControlEqn,
                                                                COFFSET_WJ::nDer );
    } );

  }
  template< template< typename U > class T, typename U=MassRateConstraint >
  static void assembleConstraintEquation( real64 const & time_n,
                                          WellControls & wellControls,
                                          T< MassRateConstraint > & constraint,
                                          WellElementSubRegion const & subRegion,
                                          string const & wellDofKey,
                                          localIndex const & rankOffset,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
  {
    // subRegion data

    localIndex const iwelemRef = subRegion.getTopWellElementIndex();
    arrayView1d< globalIndex const > const & wellElemDofNumber = subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< real64 const > const & connRate = subRegion.getField< fields::well::connectionRate >();
    arrayView2d< real64 const, compflow::USD_COMP > const & compFrac = subRegion.getField< fields::well::globalCompFraction >();
    arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens = subRegion.getField< fields::well::dGlobalCompFraction_dGlobalCompDensity >();

    // setup row/column indices for constraint equation
    using COFFSET_WJ = compositionalMultiphaseWellKernels::ColOffset_WellJac< NC, IS_THERMAL >;
    using WJ_ROFFSET = compositionalMultiphaseWellKernels::RowOffset_WellJac< NC, IS_THERMAL >;
    using Deriv = constitutive::multifluid::DerivativeOffset;

    localIndex const eqnRowIndex      = wellElemDofNumber[iwelemRef] + WJ_ROFFSET::CONTROL - rankOffset;
    globalIndex dofColIndices[COFFSET_WJ::nDer]{};
    for( integer ic = 0; ic < COFFSET_WJ::nDer; ++ic )
    {
      dofColIndices[ ic ] = wellElemDofNumber[iwelemRef] + ic;
    }

    // fluid data
    constitutive::MultiFluidBase & fluidSeparator =  wellControls.getMultiFluidSeparator();
    arrayView2d< real64 const, constitutive::multifluid::USD_FLUID > const & totalDens = fluidSeparator.totalDensity();
    arrayView3d< real64 const, constitutive::multifluid::USD_FLUID_DC > const & dTotalDens = fluidSeparator.dTotalDensity();

    // constraint data
    real64 const & targetMassRate = constraint.getConstraintValue( time_n );

    // current constraint value
    real64 const & massDensity  =
      wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::massDensityString() );

    // fix to use stored massrate
    real64 const & currentTotalVolRate =
      wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() );
    integer const useSurfaceConditions = wellControls.useSurfaceConditions();

    // residual
    real64 controlEqn =  massDensity*currentTotalVolRate - targetMassRate;

    // setup Jacobian terms
    real64 dControlEqn[NC+2+IS_THERMAL]{};

    // bring everything back to host, capture the scalars by reference
    forAll< serialPolicy >( 1, [connRate,
                                massDensity,
                                totalDens,
                                dTotalDens,
                                compFrac,
                                dCompFrac_dCompDens,
                                &dControlEqn,
                                &useSurfaceConditions,
                                localRhs,
                                controlEqn,
                                eqnRowIndex,
                                dofColIndices,
                                localMatrix,
                                &iwelemRef] ( localIndex const )
    {
      stackArray1d< real64, NC > work( NC );

      real64 const currentTotalRate = connRate[iwelemRef];

      // compute the inverse of the total density and derivatives

      real64 const totalDensInv = 1.0 / totalDens[iwelemRef][0];

      stackArray1d< real64, NC > dTotalDensInv_dCompDens( NC );
      for( integer ic = 0; ic < NC; ++ic )
      {
        dTotalDensInv_dCompDens[ic] = -dTotalDens[iwelemRef][0][Deriv::dC+ic] * totalDensInv * totalDensInv;
      }
      applyChainRuleInPlace( NC, dCompFrac_dCompDens[iwelemRef], dTotalDensInv_dCompDens, work.data() );

      // Step 2.2: divide the total mass/molar rate by the total density to get the total volumetric rate

      // Compute derivatives  dP dT
      real64 const dTotalDensInv_dPres = -dTotalDens[iwelemRef][0][Deriv::dP] * totalDensInv * totalDensInv;
      dControlEqn[COFFSET_WJ::dP] = ( useSurfaceConditions ==  0 )*massDensity * currentTotalRate * dTotalDensInv_dPres;
      if constexpr ( IS_THERMAL )
      {
        dControlEqn[COFFSET_WJ::dT] = ( useSurfaceConditions ==  0 ) * massDensity* currentTotalRate * -dTotalDens[iwelemRef][0][Deriv::dT] * totalDensInv * totalDensInv;
      }

      dControlEqn[COFFSET_WJ::dQ] = massDensity*totalDensInv;
      for( integer ic = 0; ic < NC; ++ic )
      {
        dControlEqn[COFFSET_WJ::dC+ic] = massDensity* currentTotalRate * dTotalDensInv_dCompDens[ic];
      }

      // add solver matrices
      localRhs[eqnRowIndex] += controlEqn;
      localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowIndex,
                                                                dofColIndices,
                                                                dControlEqn,
                                                                COFFSET_WJ::nDer );

    } );

  }
};


} // end namespace wellConstraintKernels

} // end namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTRAINTKERNELS_HPP
