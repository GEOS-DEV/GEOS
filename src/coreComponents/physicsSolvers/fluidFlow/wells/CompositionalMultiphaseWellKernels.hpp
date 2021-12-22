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
 * @file CompositionalMultiphaseWellKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELLKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELLKERNELS_HPP

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/relativePermeability/layouts.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"

namespace geosx
{

namespace CompositionalMultiphaseWellKernels
{

using namespace constitutive;

static constexpr real64 minDensForDivision = 1e-10;

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
  static constexpr integer DCOMP = 1;
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

  using ROFFSET = CompositionalMultiphaseWellKernels::RowOffset;
  using COFFSET = CompositionalMultiphaseWellKernels::ColOffset;

  GEOSX_HOST_DEVICE
  static void
  switchControl( WellControls::Type const & wellType,
                 WellControls::Control const & currentControl,
                 localIndex const phasePhaseIndex,
                 real64 const & targetBHP,
                 real64 const & targetPhaseRate,
                 real64 const & targetTotalRate,
                 real64 const & currentBHP,
                 arrayView1d< real64 const > const & currentPhaseVolRate,
                 real64 const & currentTotalVolRate,
                 WellControls::Control & newControl );

  template< localIndex NC >
  GEOSX_HOST_DEVICE
  static void
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
           arrayView1d< real64 > const & localRhs );

};

/******************************** FluxKernel ********************************/

struct FluxKernel
{

  using TAG = CompositionalMultiphaseWellKernels::ElemTag;
  using ROFFSET = CompositionalMultiphaseWellKernels::RowOffset;
  using COFFSET = CompositionalMultiphaseWellKernels::ColOffset;

  template< localIndex NC >
  GEOSX_HOST_DEVICE
  static void
    computeExit( real64 const & dt,
                 real64 const ( &compFlux )[NC],
                 real64 const ( &dCompFlux_dRate )[NC],
                 real64 const ( &dCompFlux_dPresUp )[NC],
                 real64 const ( &dCompFlux_dCompDensUp )[NC][NC],
                 real64 ( &oneSidedFlux )[NC],
                 real64 ( &oneSidedFluxJacobian_dRate )[NC][1],
                 real64 ( &oneSidedFluxJacobian_dPresCompUp )[NC][NC + 1] );

  template< localIndex NC >
  GEOSX_HOST_DEVICE
  static void
    compute( real64 const & dt,
             real64 const ( &compFlux )[NC],
             real64 const ( &dCompFlux_dRate )[NC],
             real64 const ( &dCompFlux_dPresUp )[NC],
             real64 const ( &dCompFlux_dCompDensUp )[NC][NC],
             real64 ( &localFlux )[2*NC],
             real64 ( &localFluxJacobian_dRate )[2*NC][1],
             real64 ( &localFluxJacobian_dPresCompUp )[2*NC][NC + 1] );

  template< localIndex NC >
  static void
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
          arrayView1d< real64 > const & localRhs );

};

/******************************** PressureRelationKernel ********************************/

struct PressureRelationKernel
{

  using TAG = CompositionalMultiphaseWellKernels::ElemTag;
  using ROFFSET = CompositionalMultiphaseWellKernels::RowOffset;
  using COFFSET = CompositionalMultiphaseWellKernels::ColOffset;

  template< localIndex NC >
  GEOSX_HOST_DEVICE
  static void
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
             real64 ( &localPresRelJacobian )[2*(NC+1)] );

  template< localIndex NC >
  static void
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
          arrayView1d< real64 > const & localRhs );

};

/******************************** PerforationKernel ********************************/

struct PerforationKernel
{

  using TAG = CompositionalMultiphaseWellKernels::SubRegionTag;

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;


  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
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
           arraySlice3d< real64 > const & dCompPerfRate_dComp );

  template< localIndex NC, localIndex NP >
  static void
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
          arrayView4d< real64 > const & dCompPerfRate_dComp );

};

/******************************** AccumulationKernel ********************************/

struct AccumulationKernel
{

  using ROFFSET = CompositionalMultiphaseWellKernels::RowOffset;
  using COFFSET = CompositionalMultiphaseWellKernels::ColOffset;

  template< localIndex NC >
  GEOSX_HOST_DEVICE
  static void
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
             real64 ( &localAccum )[NC],
             real64 ( &localAccumJacobian )[NC][NC + 1] );

  template< localIndex NC >
  static void
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
          arrayView1d< real64 > const & localRhs );

};

/******************************** VolumeBalanceKernel ********************************/

struct VolumeBalanceKernel
{

  using ROFFSET = CompositionalMultiphaseWellKernels::RowOffset;
  using COFFSET = CompositionalMultiphaseWellKernels::ColOffset;

  template< localIndex NC >
  GEOSX_HOST_DEVICE
  static void
    compute( localIndex const numPhases,
             real64 const & volume,
             arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac,
             arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & dPhaseVolFrac_dPres,
             arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dPhaseVolFrac_dComp,
             real64 & localVolBalance,
             real64 ( &localVolBalanceJacobian )[NC+1] );

  template< localIndex NC >
  static void
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
          arrayView1d< real64 > const & localRhs );

};

/******************************** PresTempCompFracInitializationKernel ********************************/

struct PresTempCompFracInitializationKernel
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

  static void
  launch( localIndex const perforationSize,
          localIndex const subRegionSize,
          localIndex const numComponents,
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
          arrayView2d< real64, compflow::USD_COMP > const & wellElemCompFrac );

};

/******************************** CompDensInitializationKernel ********************************/

struct CompDensInitializationKernel
{

  static void
  launch( localIndex const subRegionSize,
          localIndex const numComponents,
          arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompFrac,
          arrayView2d< real64 const, multifluid::USD_FLUID > const & wellElemTotalDens,
          arrayView2d< real64, compflow::USD_COMP > const & wellElemCompDens );

};

/******************************** RateInitializationKernel ********************************/

struct RateInitializationKernel
{

  static void
  launch( localIndex const subRegionSize,
          localIndex const targetPhaseIndex,
          WellControls const & wellControls,
          real64 const & currentTime,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens,
          arrayView2d< real64 const, multifluid::USD_FLUID > const & totalDens,
          arrayView1d< real64 > const & connRate );

};


/******************************** TotalMassDensityKernel ****************************/

/**
 * @class TotalMassDensityKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_PHASE number of fluid phases
 * @brief Define the interface for the property kernel in charge of computing the total mass density
 */
template< integer NUM_COMP, integer NUM_PHASE >
class TotalMassDensityKernel : public CompositionalMultiphaseBaseKernels::PropertyKernelBase< NUM_COMP >
{
public:

  using Base = CompositionalMultiphaseBaseKernels::PropertyKernelBase< NUM_COMP >;

  using Base::numComp;

  /// Compile time value for the number of phases
  static constexpr integer numPhase = NUM_PHASE;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  TotalMassDensityKernel( ObjectManagerBase & subRegion,
                          MultiFluidBase const & fluid )
    : Base(),
    m_phaseVolFrac( subRegion.getExtrinsicData< extrinsicMeshData::well::phaseVolumeFraction >() ),
    m_dPhaseVolFrac_dPres( subRegion.getExtrinsicData< extrinsicMeshData::well::dPhaseVolumeFraction_dPressure >() ),
    m_dPhaseVolFrac_dCompDens( subRegion.getExtrinsicData< extrinsicMeshData::well::dPhaseVolumeFraction_dGlobalCompDensity >() ),
    m_dCompFrac_dCompDens( subRegion.getExtrinsicData< extrinsicMeshData::well::dGlobalCompFraction_dGlobalCompDensity >() ),
    m_phaseMassDens( fluid.phaseMassDensity() ),
    m_dPhaseMassDens_dPres( fluid.dPhaseMassDensity_dPressure() ),
    m_dPhaseMassDens_dComp( fluid.dPhaseMassDensity_dGlobalCompFraction() ),
    m_totalMassDens( subRegion.getExtrinsicData< extrinsicMeshData::well::totalMassDensity >() ),
    m_dTotalMassDens_dPres( subRegion.getExtrinsicData< extrinsicMeshData::well::dTotalMassDensity_dPressure >() ),
    m_dTotalMassDens_dCompDens( subRegion.getExtrinsicData< extrinsicMeshData::well::dTotalMassDensity_dGlobalCompDensity >() )
  {}

  /**
   * @brief Compute the total mass density in an element
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[in] totalMassDensityKernelOp the function used to customize the kernel
   */
  template< typename FUNC = CompositionalMultiphaseBaseKernels::NoOpFunc >
  GEOSX_HOST_DEVICE
  void compute( localIndex const ei,
                FUNC && totalMassDensityKernelOp = CompositionalMultiphaseBaseKernels::NoOpFunc{} ) const
  {
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac = m_phaseVolFrac[ei];
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > dPhaseVolFrac_dPres = m_dPhaseVolFrac_dPres[ei];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac_dCompDens = m_dPhaseVolFrac_dCompDens[ei];
    arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];
    arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseMassDens = m_phaseMassDens[ei][0];
    arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > dPhaseMassDens_dPres = m_dPhaseMassDens_dPres[ei][0];
    arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > dPhaseMassDens_dComp = m_dPhaseMassDens_dComp[ei][0];
    real64 & totalMassDens = m_totalMassDens[ei];
    real64 & dTotalMassDens_dPres = m_dTotalMassDens_dPres[ei];
    arraySlice1d< real64, compflow::USD_FLUID_DC - 1 > dTotalMassDens_dCompDens = m_dTotalMassDens_dCompDens[ei];

    real64 dMassDens_dC[numComp]{};

    totalMassDens = 0.0;
    dTotalMassDens_dPres = 0.0;
    for( integer ic = 0; ic < numComp; ++ic )
    {
      dTotalMassDens_dCompDens[ic] = 0.0;
    }

    for( integer ip = 0; ip < numPhase; ++ip )
    {
      totalMassDens += phaseVolFrac[ip] * phaseMassDens[ip];
      dTotalMassDens_dPres += dPhaseVolFrac_dPres[ip] * phaseMassDens[ip]
                              + phaseVolFrac[ip] * dPhaseMassDens_dPres[ip];

      applyChainRule( numComp, dCompFrac_dCompDens, dPhaseMassDens_dComp[ip], dMassDens_dC );
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dTotalMassDens_dCompDens[ic] += dPhaseVolFrac_dCompDens[ip][ic] * phaseMassDens[ip]
                                        + phaseVolFrac[ip] * dMassDens_dC[ic];
      }

      totalMassDensityKernelOp( ip, totalMassDens, dTotalMassDens_dPres, dTotalMassDens_dCompDens );
    }
  }

protected:

  // inputs

  /// Views on phase volume fractions
  arrayView2d< real64 const, compflow::USD_PHASE > m_phaseVolFrac;
  arrayView2d< real64 const, compflow::USD_PHASE > m_dPhaseVolFrac_dPres;
  arrayView3d< real64 const, compflow::USD_PHASE_DC > m_dPhaseVolFrac_dCompDens;
  arrayView3d< real64 const, compflow::USD_COMP_DC > m_dCompFrac_dCompDens;

  /// Views on phase mass densities
  arrayView3d< real64 const, multifluid::USD_PHASE > m_phaseMassDens;
  arrayView3d< real64 const, multifluid::USD_PHASE > m_dPhaseMassDens_dPres;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > m_dPhaseMassDens_dComp;

  // outputs

  /// Views on total mass densities
  arrayView1d< real64 > m_totalMassDens;
  arrayView1d< real64 > m_dTotalMassDens_dPres;
  arrayView2d< real64, compflow::USD_FLUID_DC > m_dTotalMassDens_dCompDens;

};

/**
 * @class TotalMassDensityKernelFactory
 */
class TotalMassDensityKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComp the number of fluid components
   * @param[in] numPhase the number of fluid phases
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComp,
                   integer const numPhase,
                   ObjectManagerBase & subRegion,
                   MultiFluidBase const & fluid )
  {
    if( numPhase == 2 )
    {
      CompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        TotalMassDensityKernel< NUM_COMP, 2 > kernel( subRegion, fluid );
        TotalMassDensityKernel< NUM_COMP, 2 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
    else if( numPhase == 3 )
    {
      CompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        TotalMassDensityKernel< NUM_COMP, 3 > kernel( subRegion, fluid );
        TotalMassDensityKernel< NUM_COMP, 3 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
  }
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
          localIndex const numComponents,
          localIndex const numDofPerWellElement,
          localIndex const targetPhaseIndex,
          WellControls const & wellControls,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< integer const > const & wellElemGhostRank,
          arrayView1d< real64 const > wellElemVolume,
          arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseDensOld,
          arrayView1d< real64 const > const & wellElemTotalDensOld,
          real64 const & timeAtEndOfStep,
          real64 const dt,
          real64 * localResidualNorm )
  {
    using ROFFSET = CompositionalMultiphaseWellKernels::RowOffset;

    WellControls::Type const wellType = wellControls.getType();
    WellControls::Control const currentControl = wellControls.getControl();
    real64 const targetBHP = wellControls.getTargetBHP( timeAtEndOfStep );
    real64 const targetTotalRate = wellControls.getTargetTotalRate( timeAtEndOfStep );
    real64 const targetPhaseRate = wellControls.getTargetPhaseRate( timeAtEndOfStep );
    real64 const absTargetTotalRate = LvArray::math::abs( targetTotalRate );
    real64 const absTargetPhaseRate = LvArray::math::abs( targetPhaseRate );

    RAJA::ReduceSum< REDUCE_POLICY, real64 > sumScaled( 0.0 );

    forAll< POLICY >( wellElemDofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      if( wellElemGhostRank[iwelem] < 0 )
      {
        real64 normalizer = 0.0;
        for( localIndex idof = 0; idof < numDofPerWellElement; ++idof )
        {

          // Step 1: compute a normalizer for the control or pressure equation

          // for the control equation, we distinguish two cases
          if( idof == ROFFSET::CONTROL )
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
                normalizer = LvArray::math::max( absTargetTotalRate, 1e-12 );
              }
              else if( currentControl == WellControls::Control::PHASEVOLRATE )
              {
                normalizer = LvArray::math::max( absTargetPhaseRate, 1e-12 );
              }
            }
            // for the pressure difference equation, always normalize by the BHP
            else
            {
              normalizer = targetBHP;
            }
          }
          // Step 2: compute a normalizer for the mass balance equations

          else if( idof >= ROFFSET::MASSBAL && idof < ROFFSET::MASSBAL + numComponents )
          {
            if( wellType == WellControls::Type::PRODUCER ) // only PHASEVOLRATE is supported for now
            {
              normalizer = dt * absTargetPhaseRate * wellElemPhaseDensOld[iwelem][targetPhaseIndex];
            }
            else // Type::INJECTOR, only TOTALVOLRATE is supported for now
            {
              normalizer = dt * absTargetTotalRate * wellElemTotalDensOld[iwelem];
            }

            // to make sure that everything still works well if the rate is zero, we add this check
            normalizer = LvArray::math::max( normalizer, wellElemVolume[iwelem] * wellElemTotalDensOld[iwelem] );
          }
          // Step 3: compute a normalizer for the volume balance equations

          else
          {
            if( wellType == WellControls::Type::PRODUCER ) // only PHASEVOLRATE is supported for now
            {
              normalizer = dt * absTargetPhaseRate;
            }
            else // Type::INJECTOR, only TOTALVOLRATE is supported for now
            {
              normalizer = dt * absTargetTotalRate;
            }

            // to make sure that everything still works well if the rate is zero, we add this check
            normalizer = LvArray::math::max( normalizer, wellElemVolume[iwelem] );
          }

          // Step 4: compute the contribution to the residual

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
  launch( LOCAL_VECTOR const localSolution,
          globalIndex const rankOffset,
          localIndex const numComponents,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< integer const > const & wellElemGhostRank,
          arrayView1d< real64 const > const & wellElemPres,
          arrayView1d< real64 const > const & dWellElemPres,
          arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dWellElemCompDens,
          real64 const maxRelativePresChange,
          real64 const maxCompFracChange )
  {
    real64 constexpr eps = minDensForDivision;

    RAJA::ReduceMin< REDUCE_POLICY, real64 > minVal( 1.0 );

    forAll< POLICY >( wellElemDofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      if( wellElemGhostRank[iwelem] < 0 )
      {

        // the scaling of the pressures is particularly useful for the beginning of the simulation
        // with active rate control, but is useless otherwise
        real64 const pres = wellElemPres[iwelem] + dWellElemPres[iwelem];
        real64 const absPresChange = LvArray::math::abs( localSolution[wellElemDofNumber[iwelem] - rankOffset] );
        if( pres < eps )
        {
          real64 const relativePresChange = LvArray::math::abs( absPresChange ) / pres;
          if( relativePresChange > maxRelativePresChange )
          {
            minVal.min( maxRelativePresChange / relativePresChange );
          }
        }

        real64 prevTotalDens = 0;
        for( localIndex ic = 0; ic < numComponents; ++ic )
        {
          prevTotalDens += wellElemCompDens[iwelem][ic] + dWellElemCompDens[iwelem][ic];
        }

        for( localIndex ic = 0; ic < numComponents; ++ic )
        {
          localIndex const lid = wellElemDofNumber[iwelem] + ic + 1 - rankOffset;

          // compute scaling factor based on relative change in component densities
          real64 const absCompDensChange = LvArray::math::abs( localSolution[lid] );
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
  launch( LOCAL_VECTOR const localSolution,
          globalIndex const rankOffset,
          localIndex const numComponents,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< integer const > const & wellElemGhostRank,
          arrayView1d< real64 const > const & wellElemPressure,
          arrayView1d< real64 const > const & dWellElemPressure,
          arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dWellElemCompDens,
          integer const allowCompDensChopping,
          real64 const scalingFactor )
  {
    using COFFSET = CompositionalMultiphaseWellKernels::ColOffset;

    real64 constexpr eps = minDensForDivision;

    RAJA::ReduceMin< REDUCE_POLICY, localIndex > minVal( 1 );

    forAll< POLICY >( wellElemDofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      if( wellElemGhostRank[iwelem] < 0 )
      {
        // pressure
        localIndex lid = wellElemDofNumber[iwelem] + COFFSET::DPRES - rankOffset;
        real64 const newPres = wellElemPressure[iwelem] + dWellElemPressure[iwelem]
                               + scalingFactor * localSolution[lid];

        // the pressure must be positive
        if( newPres < 0.0 )
        {
          minVal.min( 0 );
        }

        // if component density is not allowed, the time step fails if a component density is negative
        // otherwise, we just check that the total density is positive, and negative component densities
        // will be chopped (i.e., set to zero) in applySystemSolution
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
