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
 * @file SinglePhaseWellKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELLKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELLKERNELS_HPP

#include "constitutive/fluid/SingleFluidExtrinsicData.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"

namespace geosx
{

namespace singlePhaseWellKernels
{

// tag to access well and reservoir elements in perforation rates computation
struct SubRegionTag
{
  static constexpr integer RES  = 0;
  static constexpr integer WELL = 1;
};

// tag to access the next and current well elements of a connection
struct ElemTag
{
  static constexpr integer CURRENT = 0;
  static constexpr integer NEXT    = 1;
};

// define the column offset of the derivatives
struct ColOffset
{
  static constexpr integer DPRES = 0;
  static constexpr integer DRATE = 1;
};

// define the row offset of the residual equations
struct RowOffset
{
  static constexpr integer CONTROL = 0;
  static constexpr integer MASSBAL = 1;
};


/******************************** ControlEquationHelper ********************************/

struct ControlEquationHelper
{

  using ROFFSET = singlePhaseWellKernels::RowOffset;
  using COFFSET = singlePhaseWellKernels::ColOffset;

  // add an epsilon to the checks to avoid control changes due to tiny pressure/rate updates
  static constexpr real64 EPS = 1e-15;

  GEOSX_HOST_DEVICE
  static void
  switchControl( WellControls::Type const & wellType,
                 WellControls::Control const & currentControl,
                 real64 const & targetBHP,
                 real64 const & targetRate,
                 real64 const & currentBHP,
                 real64 const & currentVolRate,
                 WellControls::Control & newControl );

  GEOSX_HOST_DEVICE
  static void
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
           arrayView1d< real64 > const & localRhs );

};


/******************************** FluxKernel ********************************/

struct FluxKernel
{

  using ROFFSET = singlePhaseWellKernels::RowOffset;
  using COFFSET = singlePhaseWellKernels::ColOffset;
  using TAG = singlePhaseWellKernels::ElemTag;

  static void
  launch( localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & connRate,
          arrayView1d< real64 const > const & dConnRate,
          real64 const & dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );

};


/******************************** PressureRelationKernel ********************************/

struct PressureRelationKernel
{

  using ROFFSET = singlePhaseWellKernels::RowOffset;
  using COFFSET = singlePhaseWellKernels::ColOffset;
  using TAG = singlePhaseWellKernels::ElemTag;

  static localIndex
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
          arrayView1d< real64 > const & localRhs );

};


/******************************** PerforationKernel ********************************/

struct PerforationKernel
{

  using TAG = singlePhaseWellKernels::SubRegionTag;

  using SinglePhaseFlowAccessors =
    StencilAccessors< extrinsicMeshData::flow::pressure,
                      extrinsicMeshData::flow::deltaPressure >;

  using SingleFluidAccessors =
    StencilMaterialAccessors< constitutive::SingleFluidBase,
                              extrinsicMeshData::singlefluid::density,
                              extrinsicMeshData::singlefluid::dDensity_dPressure,
                              extrinsicMeshData::singlefluid::viscosity,
                              extrinsicMeshData::singlefluid::dViscosity_dPressure >;

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  GEOSX_HOST_DEVICE
  static void
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
           arraySlice1d< real64 > const & dPerfRate_dPres );

  static void
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
          arrayView2d< real64 > const & dPerfRate_dPres );

};

/******************************** AccumulationKernel ********************************/

struct AccumulationKernel
{

  using ROFFSET = singlePhaseWellKernels::RowOffset;
  using COFFSET = singlePhaseWellKernels::ColOffset;

  static void
  launch( localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< integer const > const & wellElemGhostRank,
          arrayView1d< real64 const > const & wellElemVolume,
          arrayView2d< real64 const > const & wellElemDensity,
          arrayView2d< real64 const > const & dWellElemDensity_dPres,
          arrayView1d< real64 const > const & wellElemDensityOld,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );

};

/******************************** PressureInitializationKernel ********************************/

struct PresInitializationKernel
{

  using SinglePhaseFlowAccessors =
    StencilAccessors< extrinsicMeshData::flow::pressure >;

  using SingleFluidAccessors =
    StencilMaterialAccessors< constitutive::SingleFluidBase,
                              extrinsicMeshData::singlefluid::density >;

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  static void
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
          arrayView1d< real64 > const & wellElemPressure );

};

/******************************** RateInitializationKernel ********************************/

struct RateInitializationKernel
{

  static void
  launch( localIndex const subRegionSize,
          WellControls const & wellControls,
          real64 const & currentTime,
          arrayView2d< real64 const > const & wellElemDens,
          arrayView1d< real64 > const & connRate );

};


/******************************** ResidualNormKernel ********************************/

struct ResidualNormKernel
{

  template< typename POLICY, typename REDUCE_POLICY, typename LOCAL_VECTOR >
  static void
  launch( LOCAL_VECTOR const localResidual,
          globalIndex const rankOffset,
          bool const isLocallyOwned,
          localIndex const iwelemControl,
          WellControls const & wellControls,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< integer const > const & wellElemGhostRank,
          arrayView1d< real64 const > wellElemVolume,
          arrayView1d< real64 const > const & wellElemDensOld,
          real64 const & timeAtEndOfStep,
          real64 const dt,
          real64 * localResidualNorm )
  {
    WellControls::Control const currentControl = wellControls.getControl();
    real64 const targetBHP = wellControls.getTargetBHP( timeAtEndOfStep );
    real64 const targetRate = wellControls.getTargetTotalRate( timeAtEndOfStep );
    real64 const absTargetRate = LvArray::math::abs( targetRate );

    RAJA::ReduceSum< REDUCE_POLICY, real64 > sumScaled( 0.0 );

    forAll< POLICY >( wellElemDofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      if( wellElemGhostRank[iwelem] < 0 )
      {
        for( localIndex idof = 0; idof < 2; ++idof )
        {
          real64 normalizer = 0.0;
          if( idof == singlePhaseWellKernels::RowOffset::CONTROL )
          {
            // for the top well element, normalize using the current control
            if( isLocallyOwned && iwelem == iwelemControl )
            {
              if( currentControl == WellControls::Control::BHP )
              {
                normalizer = targetBHP;
              }
              else if( currentControl == WellControls::Control::TOTALVOLRATE )
              {
                normalizer = LvArray::math::max( absTargetRate, 1e-12 );
              }
            }
            // for the pressure difference equation, always normalize by the BHP
            else
            {
              normalizer = targetBHP;
            }
          }
          else // SinglePhaseWell::RowOffset::MASSBAL
          {
            normalizer = dt * absTargetRate * wellElemDensOld[iwelem];

            // to make sure that everything still works well if the rate is zero, we add this check
            normalizer = LvArray::math::max( normalizer, wellElemVolume[iwelem] * wellElemDensOld[iwelem] );
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

/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{
  template< typename POLICY, typename REDUCE_POLICY, typename LOCAL_VECTOR >
  static localIndex
  launch( LOCAL_VECTOR const localSolution,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & presDofNumber,
          arrayView1d< integer const > const & ghostRank,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & dPres,
          real64 const scalingFactor )
  {
    RAJA::ReduceMin< REDUCE_POLICY, localIndex > minVal( 1 );

    forAll< POLICY >( presDofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 && presDofNumber[ei] >= 0 )
      {
        localIndex const lid = presDofNumber[ei] - rankOffset;
        real64 const newPres = pres[ei] + dPres[ei] + scalingFactor * localSolution[lid];

        if( newPres < 0.0 )
        {
          minVal.min( 0 );
        }
      }
    } );
    return minVal.get();
  }
};

} // end namespace singlePhaseWellKernels

} // end namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELLKERNELS_HPP
