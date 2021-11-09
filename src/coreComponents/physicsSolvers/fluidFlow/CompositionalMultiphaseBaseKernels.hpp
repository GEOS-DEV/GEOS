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
 * @file CompositionalMultiphaseBaseKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEKERNELS_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/layouts.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "functions/TableFunction.hpp"
#include "mesh/ElementRegionManager.hpp"


namespace geosx
{

namespace CompositionalMultiphaseBaseKernels
{

using namespace constitutive;

static constexpr real64 minDensForDivision = 1e-10;

/******************************** ComponentFractionKernel ********************************/

/**
 * @brief Functions to compute component fractions from global component densities (mass or molar)
 */
struct ComponentFractionKernel
{
  template< localIndex NC >
  GEOSX_HOST_DEVICE
  static void
  compute( arraySlice1d< real64 const, compflow::USD_COMP - 1 > compDens,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > dCompDens,
           arraySlice1d< real64, compflow::USD_COMP - 1 > compFrac,
           arraySlice2d< real64, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens );

  template< localIndex NC >
  static void
  launch( localIndex const size,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens,
          arrayView2d< real64, compflow::USD_COMP > const & compFrac,
          arrayView3d< real64, compflow::USD_COMP_DC > const & dCompFrac_dCompDens );

  template< localIndex NC >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens,
          arrayView2d< real64, compflow::USD_COMP > const & compFrac,
          arrayView3d< real64, compflow::USD_COMP_DC > const & dCompFrac_dCompDens );
};

/******************************** PhaseVolumeFractionKernel ********************************/

/**
 * @brief Functions to compute phase volume fractions (saturations) and derivatives
 */
struct PhaseVolumeFractionKernel
{
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
  compute( arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & compDens,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & dCompDens,
           arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dCompFrac_dCompDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & dPhaseDens_dPres,
           arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dPhaseDens_dComp,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseFrac,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & dPhaseFrac_dPres,
           arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dPhaseFrac_dComp,
           arraySlice1d< real64, compflow::USD_PHASE - 1 > const & phaseVolFrac,
           arraySlice1d< real64, compflow::USD_PHASE - 1 > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64, compflow::USD_PHASE_DC - 1 > const & dPhaseVolFrac_dComp );

  template< localIndex NC, localIndex NP >
  static void
  launch( localIndex const size,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens_dComp,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseFrac,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseFrac_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseFrac_dComp,
          arrayView2d< real64, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView2d< real64, compflow::USD_PHASE > const & dPhaseVolFrac_dPres,
          arrayView3d< real64, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dComp );

  template< localIndex NC, localIndex NP >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens_dComp,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseFrac,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseFrac_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseFrac_dComp,
          arrayView2d< real64, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView2d< real64, compflow::USD_PHASE > const & dPhaseVolFrac_dPres,
          arrayView3d< real64, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dComp );
};


/******************************** FluidUpdateKernel ********************************/

struct FluidUpdateKernel
{
  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  launch( localIndex const size,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & temp,
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k], temp[k], compFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  launch( localIndex const size,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & dPres,
          arrayView1d< real64 const > const & temp,
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k] + dPres[k], temp[k], compFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & dPres,
          arrayView1d< real64 const > const & temp,
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k] + dPres[k], temp[k], compFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & temp,
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k], temp[k], compFrac[k] );
      }
    } );
  }
};

/******************************** RelativePermeabilityUpdateKernel ********************************/

struct RelativePermeabilityUpdateKernel
{
  template< typename POLICY, typename RELPERM_WRAPPER >
  static void
  launch( localIndex const size,
          RELPERM_WRAPPER const & relPermWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < relPermWrapper.numGauss(); ++q )
      {
        relPermWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename RELPERM_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          RELPERM_WRAPPER const & relPermWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < relPermWrapper.numGauss(); ++q )
      {
        relPermWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }
};

/******************************** CapillaryPressureUpdateKernel ********************************/

struct CapillaryPressureUpdateKernel
{
  template< typename POLICY, typename CAPPRES_WRAPPER >
  static void
  launch( localIndex const size,
          CAPPRES_WRAPPER const & capPresWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < capPresWrapper.numGauss(); ++q )
      {
        capPresWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename CAPPRES_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          CAPPRES_WRAPPER const & capPresWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < capPresWrapper.numGauss(); ++q )
      {
        capPresWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }
};

/******************************** AccumulationKernel ********************************/

/**
 * @brief Functions to assemble accumulation term contributions to residual and Jacobian
 */
struct AccumulationKernel
{
  template< localIndex NC >
  GEOSX_HOST_DEVICE
  static void
    compute( localIndex const numPhases,
             real64 const & poreVolOld,
             real64 const & poreVolNew,
             real64 const & dPoreVol_dP,
             arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dCompFrac_dCompDens,
             arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFracOld,
             arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac,
             arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & dPhaseVolFrac_dPres,
             arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dPhaseVolFrac_dCompDens,
             arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseDensOld,
             arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseDens,
             arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & dPhaseDens_dPres,
             arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dPhaseDens_dComp,
             arraySlice2d< real64 const, compflow::USD_PHASE_COMP-1 > const & phaseCompFracOld,
             arraySlice2d< real64 const, multifluid::USD_PHASE_COMP-2 > const & phaseCompFrac,
             arraySlice2d< real64 const, multifluid::USD_PHASE_COMP-2 > const & dPhaseCompFrac_dPres,
             arraySlice3d< real64 const, multifluid::USD_PHASE_COMP_DC-2 > const & dPhaseCompFrac_dComp,
             real64 ( &localAccum )[NC],
             real64 ( &localAccumJacobian )[NC][NC+1] );

  template< localIndex NC >
  static void
  launch( localIndex const numPhases,
          localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView1d< real64 const > const & volume,
          arrayView2d< real64 const > const & porosityOld,
          arrayView2d< real64 const > const & porosityNew,
          arrayView2d< real64 const > const & dPoro_dPres,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFracOld,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dCompDens,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseDensOld,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens_dComp,
          arrayView3d< real64 const, compflow::USD_PHASE_COMP > const & phaseCompFracOld,
          arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & phaseCompFrac,
          arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & dPhaseCompFrac_dPres,
          arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > const & dPhaseCompFrac_dComp,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );
};

/******************************** VolumeBalanceKernel ********************************/

/**
 * @brief Functions to assemble volume balance contributions to residual and Jacobian
 */
struct VolumeBalanceKernel
{
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
  compute( real64 const & volume,
           real64 const & porosityNew,
           real64 const & dPoro_dPres,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dPhaseVolFrac_dCompDens,
           real64 & localVolBalance,
           real64 * const localVolBalanceJacobian );

  template< localIndex NC, localIndex NP >
  static void
  launch( localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView1d< real64 const > const & volume,
          arrayView2d< real64 const > const & porosityNew,
          arrayView2d< real64 const > const & dPoro_dPres,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dCompDens,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );
};

/******************************** ResidualNormKernel ********************************/

struct ResidualNormKernel
{

  template< typename POLICY, typename REDUCE_POLICY >
  static void launch( arrayView1d< real64 const > const & localResidual,
                      globalIndex const rankOffset,
                      localIndex const numComponents,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< integer const > const & ghostRank,
                      arrayView1d< real64 const > const & refPoro,
                      arrayView1d< real64 const > const & volume,
                      arrayView1d< real64 const > const & totalDensOld,
                      real64 & localResidualNorm )
  {
    RAJA::ReduceSum< REDUCE_POLICY, real64 > localSum( 0.0 );

    forAll< POLICY >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 )
      {
        localIndex const localRow = dofNumber[ei] - rankOffset;
        real64 const normalizer = totalDensOld[ei] * refPoro[ei] * volume[ei];

        for( localIndex idof = 0; idof < numComponents + 1; ++idof )
        {
          real64 const val = localResidual[localRow + idof] / normalizer;
          localSum += val * val;
        }
      }
    } );
    localResidualNorm += localSum.get();
  }

};


/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{
  template< typename POLICY, typename REDUCE_POLICY >
  static localIndex
  launch( arrayView1d< real64 const > const & localSolution,
          globalIndex const rankOffset,
          localIndex const numComponents,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & ghostRank,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & dPres,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens,
          integer const allowCompDensChopping,
          real64 const scalingFactor )
  {
    real64 constexpr eps = minDensForDivision;

    RAJA::ReduceMin< REDUCE_POLICY, integer > check( 1 );

    forAll< POLICY >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 )
      {
        localIndex const localRow = dofNumber[ei] - rankOffset;
        {
          real64 const newPres = pres[ei] + dPres[ei] + scalingFactor * localSolution[localRow];
          check.min( newPres >= 0.0 );
        }

        // if component density chopping is not allowed, the time step fails if a component density is negative
        // otherwise, we just check that the total density is positive, and negative component densities
        // will be chopped (i.e., set to zero) in ApplySystemSolution)
        if( !allowCompDensChopping )
        {
          for( localIndex ic = 0; ic < numComponents; ++ic )
          {
            real64 const newDens = compDens[ei][ic] + dCompDens[ei][ic] + scalingFactor * localSolution[localRow + ic + 1];
            check.min( newDens >= 0.0 );
          }
        }
        else
        {
          real64 totalDens = 0.0;
          for( localIndex ic = 0; ic < numComponents; ++ic )
          {
            real64 const newDens = compDens[ei][ic] + dCompDens[ei][ic] + scalingFactor * localSolution[localRow + ic + 1];
            totalDens += (newDens > 0.0) ? newDens : 0.0;
          }
          check.min( totalDens >= eps );
        }
      }
    } );
    return check.get();
  }

};

/******************************** HydrostaticPressureKernel ********************************/

struct HydrostaticPressureKernel
{

  // TODO: this type of constants should be centralized somewhere or provided by fluid model
  static real64 constexpr MIN_FOR_PHASE_PRESENCE = 1e-12;

  enum class ReturnType : integer
  {
    FAILED_TO_CONVERGE = 0,
    DETECTED_MULTIPHASE_FLOW = 1,
    SUCCESS = 2
  };

  template< typename FLUID_WRAPPER >
  static ReturnType
  computeHydrostaticPressure( localIndex const numComps,
                              localIndex const numPhases,
                              integer const ipInit,
                              integer const maxNumEquilIterations,
                              real64 const & equilTolerance,
                              real64 const (&gravVector)[ 3 ],
                              FLUID_WRAPPER fluidWrapper,
                              arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
                              TableFunction::KernelWrapper tempTableWrapper,
                              real64 const & refElevation,
                              real64 const & refPres,
                              arraySlice1d< real64 const > const & refPhaseMassDens,
                              real64 const & newElevation,
                              real64 & newPres,
                              arraySlice1d< real64 > const & newPhaseMassDens )
  {
    // fluid properties at this elevation
    StackArray< real64, 2, constitutive::MultiFluidBase::MAX_NUM_COMPONENTS, compflow::LAYOUT_COMP > compFrac( 1, numComps );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, multifluid::LAYOUT_PHASE > phaseFrac( 1, 1, numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, multifluid::LAYOUT_PHASE > phaseDens( 1, 1, numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, multifluid::LAYOUT_PHASE > phaseMassDens( 1, 1, numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, multifluid::LAYOUT_PHASE > phaseVisc( 1, 1, numPhases );
    StackArray< real64, 4, constitutive::MultiFluidBase::MAX_NUM_PHASES *constitutive::MultiFluidBase::MAX_NUM_COMPONENTS,
                multifluid::LAYOUT_PHASE_COMP > phaseCompFrac( 1, 1, numPhases, numComps );
    real64 totalDens = 0.0;

    bool isSinglePhaseFlow = true;

    // Step 1: compute the hydrostatic pressure at the current elevation

    real64 const gravCoef = gravVector[2] * ( refElevation - newElevation );
    real64 const temp = tempTableWrapper.compute( &newElevation );
    for( localIndex ic = 0; ic < numComps; ++ic )
    {
      compFrac[0][ic] = compFracTableWrappers[ic].compute( &newElevation );
    }

    // Step 2: guess the pressure with the refPhaseMassDensity

    real64 pres0 = refPres - refPhaseMassDens[ipInit] * gravCoef;
    real64 pres1 = 0.0;

    // Step 3: compute the mass density at this elevation using the guess, and update pressure

    fluidWrapper.compute( pres0,
                          temp,
                          compFrac[0],
                          phaseFrac[0][0],
                          phaseDens[0][0],
                          phaseMassDens[0][0],
                          phaseVisc[0][0],
                          phaseCompFrac[0][0],
                          totalDens );
    pres1 = refPres - 0.5 * ( refPhaseMassDens[ipInit] + phaseMassDens[0][0][ipInit] ) * gravCoef;

    // Step 4: fixed-point iteration until convergence

    bool equilHasConverged = false;
    for( localIndex eqIter = 0; eqIter < maxNumEquilIterations; ++eqIter )
    {

      // check convergence
      equilHasConverged = ( LvArray::math::abs( pres0 - pres1 ) < equilTolerance );
      pres0 = pres1;

      // if converged, check number of phases and move on
      if( equilHasConverged )
      {
        // make sure that the fluid is single-phase, other we have to issue a warning (for now)
        // if only one phase is mobile, we are in good shape (unfortunately it is hard to access relperm from here)
        localIndex numberOfPhases = 0;
        for( localIndex ip = 0; ip < numPhases; ++ip )
        {
          if( phaseFrac[0][0][ip] > MIN_FOR_PHASE_PRESENCE )
          {
            numberOfPhases++;
          }
        }
        if( numberOfPhases > 1 )
        {
          isSinglePhaseFlow = false;
        }

        break;
      }

      // compute the mass density at this elevation using the previous pressure, and compute the new pressure
      fluidWrapper.compute( pres0,
                            temp,
                            compFrac[0],
                            phaseFrac[0][0],
                            phaseDens[0][0],
                            phaseMassDens[0][0],
                            phaseVisc[0][0],
                            phaseCompFrac[0][0],
                            totalDens );
      pres1 = refPres - 0.5 * ( refPhaseMassDens[ipInit] + phaseMassDens[0][0][ipInit] ) * gravCoef;
    }

    // Step 5: save the hydrostatic pressure and the corresponding density

    newPres = pres1;
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      newPhaseMassDens[ip] = phaseMassDens[0][0][ip];
    }

    if( !equilHasConverged )
    {
      return ReturnType::FAILED_TO_CONVERGE;
    }
    else if( !isSinglePhaseFlow )
    {
      return ReturnType::DETECTED_MULTIPHASE_FLOW;
    }
    else
    {
      return ReturnType::SUCCESS;
    }
  }

  template< typename FLUID_WRAPPER >
  static ReturnType
  launch( localIndex const size,
          localIndex const numComps,
          localIndex const numPhases,
          integer const ipInit,
          integer const maxNumEquilIterations,
          real64 const equilTolerance,
          real64 const (&gravVector)[ 3 ],
          real64 const & minElevation,
          real64 const & elevationIncrement,
          real64 const & datumElevation,
          real64 const & datumPres,
          FLUID_WRAPPER fluidWrapper,
          arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
          TableFunction::KernelWrapper tempTableWrapper,
          arrayView1d< arrayView1d< real64 > const > elevationValues,
          arrayView1d< real64 > pressureValues )
  {

    ReturnType returnVal = ReturnType::SUCCESS;

    // Step 1: compute the phase mass densities at datum

    // datum fluid properties
    array2d< real64, compflow::LAYOUT_COMP > datumCompFrac( 1, numComps );
    array3d< real64, multifluid::LAYOUT_PHASE > datumPhaseFrac( 1, 1, numPhases );
    array3d< real64, multifluid::LAYOUT_PHASE > datumPhaseDens( 1, 1, numPhases );
    array3d< real64, multifluid::LAYOUT_PHASE > datumPhaseMassDens( 1, 1, numPhases );
    array3d< real64, multifluid::LAYOUT_PHASE > datumPhaseVisc( 1, 1, numPhases );
    array4d< real64, multifluid::LAYOUT_PHASE_COMP > datumPhaseCompFrac( 1, 1, numPhases, numComps );
    real64 datumTotalDens = 0.0;

    real64 const datumTemp = tempTableWrapper.compute( &datumElevation );
    for( localIndex ic = 0; ic < numComps; ++ic )
    {
      datumCompFrac[0][ic] = compFracTableWrappers[ic].compute( &datumElevation );
    }
    fluidWrapper.compute( datumPres,
                          datumTemp,
                          datumCompFrac[0],
                          datumPhaseFrac[0][0],
                          datumPhaseDens[0][0],
                          datumPhaseMassDens[0][0],
                          datumPhaseVisc[0][0],
                          datumPhaseCompFrac[0][0],
                          datumTotalDens );

    // Step 2: find the closest elevation to datumElevation

    forAll< parallelHostPolicy >( size, [=] ( localIndex const i )
    {
      real64 const elevation = minElevation + i * elevationIncrement;
      elevationValues[0][i] = elevation;
    } );
    integer const iRef = LvArray::sortedArrayManipulation::find( elevationValues[0].begin(),
                                                                 elevationValues[0].size(),
                                                                 datumElevation );

    // Step 3: compute the mass density and pressure at the reference elevation

    array2d< real64 > phaseMassDens( pressureValues.size(), numPhases );
    // temporary array without permutation to compile on Lassen
    array1d< real64 > datumPhaseMassDensTmp( numPhases );
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      datumPhaseMassDensTmp[ip] = datumPhaseMassDens[0][0][ip];
    }

    ReturnType const refReturnVal =
      computeHydrostaticPressure( numComps,
                                  numPhases,
                                  ipInit,
                                  maxNumEquilIterations,
                                  equilTolerance,
                                  gravVector,
                                  fluidWrapper,
                                  compFracTableWrappers,
                                  tempTableWrapper,
                                  datumElevation,
                                  datumPres,
                                  datumPhaseMassDensTmp,
                                  elevationValues[0][iRef],
                                  pressureValues[iRef],
                                  phaseMassDens[iRef] );
    if( refReturnVal == ReturnType::FAILED_TO_CONVERGE )
    {
      return ReturnType::FAILED_TO_CONVERGE;
    }
    else if( refReturnVal == ReturnType::DETECTED_MULTIPHASE_FLOW )
    {
      returnVal = ReturnType::DETECTED_MULTIPHASE_FLOW;
    }

    // Step 4: for each elevation above the reference elevation, compute the pressure

    localIndex const numEntriesAboveRef = size - iRef - 1;
    forAll< serialPolicy >( numEntriesAboveRef, [=, &returnVal] ( localIndex const i )
    {
      ReturnType const returnValAboveRef =
        computeHydrostaticPressure( numComps,
                                    numPhases,
                                    ipInit,
                                    maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    compFracTableWrappers,
                                    tempTableWrapper,
                                    elevationValues[0][iRef+i],
                                    pressureValues[iRef+i],
                                    phaseMassDens[iRef+i],
                                    elevationValues[0][iRef+i+1],
                                    pressureValues[iRef+i+1],
                                    phaseMassDens[iRef+i+1] );
      if( returnValAboveRef == ReturnType::FAILED_TO_CONVERGE )
      {
        returnVal = ReturnType::FAILED_TO_CONVERGE;
      }
      else if( ( returnValAboveRef == ReturnType::DETECTED_MULTIPHASE_FLOW ) &&
               ( returnVal != ReturnType::FAILED_TO_CONVERGE ) )
      {
        returnVal = ReturnType::DETECTED_MULTIPHASE_FLOW;
      }

    } );

    // Step 5: for each elevation below the reference elevation, compute the pressure

    localIndex const numEntriesBelowRef = iRef;
    forAll< serialPolicy >( numEntriesBelowRef, [=, &returnVal] ( localIndex const i )
    {
      ReturnType const returnValBelowRef =
        computeHydrostaticPressure( numComps,
                                    numPhases,
                                    ipInit,
                                    maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    compFracTableWrappers,
                                    tempTableWrapper,
                                    elevationValues[0][iRef-i],
                                    pressureValues[iRef-i],
                                    phaseMassDens[iRef-i],
                                    elevationValues[0][iRef-i-1],
                                    pressureValues[iRef-i-1],
                                    phaseMassDens[iRef-i-1] );
      if( returnValBelowRef == ReturnType::FAILED_TO_CONVERGE )
      {
        returnVal = ReturnType::FAILED_TO_CONVERGE;
      }
      else if( ( returnValBelowRef == ReturnType::DETECTED_MULTIPHASE_FLOW ) &&
               ( returnVal != ReturnType::FAILED_TO_CONVERGE ) )
      {
        returnVal = ReturnType::DETECTED_MULTIPHASE_FLOW;
      }

    } );

    return returnVal;
  }

};


/******************************** Kernel launch machinery ********************************/

namespace internal
{

template< typename T, typename LAMBDA >
void KernelLaunchSelectorCompSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "KernelLaunchSelectorCompSwitch: type should be integral" );

  switch( value )
  {
    case 1:
    { lambda( std::integral_constant< T, 1 >() ); return; }
    case 2:
    { lambda( std::integral_constant< T, 2 >() ); return; }
    case 3:
    { lambda( std::integral_constant< T, 3 >() ); return; }
    case 4:
    { lambda( std::integral_constant< T, 4 >() ); return; }
    case 5:
    { lambda( std::integral_constant< T, 5 >() ); return; }
    default:
    { GEOSX_ERROR( "Unsupported number of components: " << value ); }
  }
}

} // namespace helpers

template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector1( localIndex const numComp, ARGS && ... args )
{
  internal::KernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
  {
    KERNELWRAPPER::template launch< NC() >( std::forward< ARGS >( args )... );
  } );
}

template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector2( localIndex const numComp, localIndex const numPhase, ARGS && ... args )
{
  // Ideally this would be inside the dispatch, but it breaks on Summit with GCC 9.1.0 and CUDA 11.0.3.
  if( numPhase == 2 )
  {
    internal::KernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
    {
      KERNELWRAPPER::template launch< NC(), 2 >( std::forward< ARGS >( args ) ... );
    } );
  }
  else if( numPhase == 3 )
  {
    internal::KernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
    {
      KERNELWRAPPER::template launch< NC(), 3 >( std::forward< ARGS >( args ) ... );
    } );
  }
  else
  {
    GEOSX_ERROR( "Unsupported number of phases: " << numPhase );
  }
}

} // namespace CompositionalMultiphaseBaseKernels

} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEKERNELS_HPP
