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
 * @file ThermalSinglePhaseMixedMFDKernels.hpp
 *
 * Energy equations of the mixed mimetic single-phase solver. Unknowns: face heat flux q_f, cell temperature T_K
 * and cell specific enthalpy h_K. With sigma_f the orientation of the face f relative to the cell K and m_f the
 * face mass flux:
 *
 *   Fourier law      sum_j (M_e)_ij sigma_j q_j - T_K + T_f = 0      M_e the inner product weighted by K_e^{-1}
 *   energy balance   E_K - E_K^n + dt sum_f sigma_f ( m_f h_f^up + q_f ) = 0
 *   closure          h_K - h_eos( p_K, T_K ) = 0
 *
 * h_f^up is the enthalpy of the upstream cell, or the prescribed value on an inflow boundary face.
 * The heat operator uses the cell classification eta of the flow operator.
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMALMIXEDMFDKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMALMIXEDMFDKERNELS_HPP

#include "constitutive/thermalConductivity/SinglePhaseThermalConductivityBase.hpp"
#include "constitutive/thermalConductivity/ThermalConductivityFields.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/SinglePhaseMixedMFDKernels.hpp"

namespace geos
{

namespace thermalSinglePhaseMixedMFDKernels
{

/// components of the cell unknown: p_K, T_K, h_K
struct CellDof
{
  static constexpr integer pressure = 0;
  static constexpr integer temperature = 1;
  static constexpr integer enthalpy = 2;
  static constexpr integer num = 3;
};

/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @class ElementBasedAssemblyKernel
 * @brief Cell contributions: temperature derivatives of the isothermal equations, Fourier law on the
 *        non-condensed faces, energy balance.
 */
template< integer NUM_FACE, typename IP >
class ElementBasedAssemblyKernel : public singlePhaseMixedMFDKernels::ElementBasedAssemblyKernel< NUM_FACE, IP >
{
public:

  using Base = singlePhaseMixedMFDKernels::ElementBasedAssemblyKernel< NUM_FACE, IP >;
  using ThermalDerivOffset = constitutive::singlefluid::DerivativeOffsetC< 1 >;

  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using Base::m_rankOffset;
  using Base::m_dt;
  using Base::m_elemGhostRank;
  using Base::m_elemDofNumber;
  using Base::m_faceGhostRank;
  using Base::m_elemToFaces;
  using Base::m_elemCenter;
  using Base::m_elemVolume;
  using Base::m_faceToNodes;
  using Base::m_elemRegionList;
  using Base::m_elemSubRegionList;
  using Base::m_elemList;
  using Base::m_nodePosition;
  using Base::m_lengthTolerance;
  using Base::m_mfdFlag;
  using Base::m_myElemLocalToGlobal;
  using Base::m_faceStencilLabel;
  using Base::m_mob;
  using Base::m_dMob;
  using Base::m_dElemDens;
  using Base::m_localMatrix;
  using Base::m_localRhs;

  ElementBasedAssemblyKernel( globalIndex const rankOffset,
                              real64 const & lengthTolerance,
                              string const elemDofKey,
                              string const faceDofKey,
                              string const faceHeatDofKey,
                              NodeManager const & nodeManager,
                              FaceManager const & faceManager,
                              ElementRegionManager const & elemManager,
                              CellElementSubRegion const & subRegion,
                              constitutive::SingleFluidBase const & fluid,
                              constitutive::PermeabilityBase const & permeability,
                              constitutive::SinglePhaseThermalConductivityBase const & conductivity,
                              ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const & enthalpyAcc,
                              ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const & elemDofNumberAcc,
                              ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const & elemLocalToGlobalAcc,
                              SortedArrayView< localIndex const > const & regionFilter,
                              real64 const & dt,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
    : Base( rankOffset, lengthTolerance, elemDofKey, faceDofKey, nodeManager, faceManager,
            subRegion, fluid, permeability, dt, localMatrix, localRhs ),
    m_faceHeatDofNumber( faceManager.getReference< array1d< globalIndex > >( faceHeatDofKey ) ),
    m_faceArea( faceManager.faceArea() ),
    m_faceHeatFlux( faceManager.getField< fields::mixedMimetic::faceHeatFlux >() ),
    m_heatBc{ faceManager.getField< fields::mixedMimetic::heatBoundaryType >(),
             faceManager.getField< fields::mixedMimetic::heatBoundaryValue >(),
             faceManager.getField< fields::mixedMimetic::heatBoundaryCoefficient >() },
    m_isEnthalpyBcFace( faceManager.getField< fields::mixedMimetic::isEnthalpyBcFace >() ),
    m_bcEnthalpy( faceManager.getField< fields::mixedMimetic::bcEnthalpy >() ),
    m_elemTemp( subRegion.getField< fields::flow::temperature >() ),
    m_elemEnthalpy( subRegion.getField< fields::mixedMimetic::enthalpy >() ),
    m_effCond( conductivity.effectiveConductivity() ),
    m_regionFilter( regionFilter ),
    m_enthalpyAcc( enthalpyAcc.toNestedViewConst() ),
    m_elemDofNumberAcc( elemDofNumberAcc.toNestedViewConst() ),
    m_elemLocalToGlobalAcc( elemLocalToGlobalAcc.toNestedViewConst() )
  {
    GEOS_UNUSED_VAR( elemManager );
  }

  /**
   * @struct StackVariables
   * @brief Stack variables of the energy equations
   */
  struct StackVariables : public Base::StackVariables
  {
    GEOS_HOST_DEVICE
    StackVariables()
      : Base::StackVariables(),
                                       heatMatrix( NUM_FACE, NUM_FACE )
    {}

    /// M_e = chi M_mfd( K_e ) + ( 1 - chi ) M_tpfa( K_e )
    stackArray2d< real64, NUM_FACE *NUM_FACE > heatMatrix;

    /// d( Darcy residual )/dT_K
    real64 dFaceResidual_dTemp[NUM_FACE]{};

    /// 1 if q_f is prescribed (Neumann); 1 if q_f is condensed into a two-point relation
    integer isEssentialHeatFace[NUM_FACE]{};
    integer isCondensedHeatFace[NUM_FACE]{};

    /// sigma_i q_i, d( sigma_i q_i )/dq_i (0 if prescribed), Fourier residuals and their derivatives
    real64 localHeatFlux[NUM_FACE]{};
    real64 dLocalHeatFlux[NUM_FACE]{};
    real64 heatResidual[NUM_FACE]{};
    real64 dHeatResidual_dHeatFlux[NUM_FACE][NUM_FACE]{};

    /// dt sum_f sigma_f ( m_f h_f^up + q_f ) and its derivatives
    real64 energyFlux = 0.0;
    real64 dEnergyFlux_dFlux[NUM_FACE]{};
    real64 dEnergyFlux_dHeatFlux[NUM_FACE]{};
    real64 dEnergyFlux_dOwnH = 0.0;
    real64 dEnergyFlux_dNbrH[NUM_FACE]{};
    globalIndex nbrHDofColIndex[NUM_FACE]{};

    localIndex faceHeatEqnRowIndex[NUM_FACE]{};
    globalIndex faceHeatDofColIndices[NUM_FACE]{};
  };

  GEOS_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    Base::setup( ei, stack );
    for( integer i = 0; i < NUM_FACE; ++i )
    {
      localIndex const kf = m_elemToFaces[ei][i];
      stack.faceHeatEqnRowIndex[i] = m_faceHeatDofNumber[kf] - m_rankOffset;
      stack.faceHeatDofColIndices[i] = m_faceHeatDofNumber[kf];
      stack.nbrHDofColIndex[i] = -1;
    }
  }

  GEOS_HOST_DEVICE
  void compute( localIndex const ei,
                StackVariables & stack ) const
  {
    Base::compute( ei, stack );

    // d( Darcy residual )/dT_K = d( mu / rho )/dT ( M sigma m )_i + d( rho )/dT ( gamma_K - gamma_f )
    real64 const invMob = 1.0 / m_mob[ei];
    real64 const dInvMob_dTemp = -m_dMob[ei][ThermalDerivOffset::dT] * invMob * invMob;
    real64 const dDens_dTemp = m_dElemDens[ei][0][ThermalDerivOffset::dT];
    for( integer i = 0; i < NUM_FACE; ++i )
    {
      stack.dFaceResidual_dTemp[i] = dInvMob_dTemp * stack.mDotFlux[i] + dDens_dTemp * stack.gravCoefDif[i];
    }

    // status of q_f: prescribed, condensed or non-condensed
    bool anyLiveHeatFace = false;
    for( integer i = 0; i < NUM_FACE; ++i )
    {
      localIndex const kf = m_elemToFaces[ei][i];
      stack.isEssentialHeatFace[i] = m_heatBc[kf].isEssential() ? 1 : 0;
      stack.isCondensedHeatFace[i] = ( m_faceStencilLabel[kf] == 0 && stack.isEssentialHeatFace[i] == 0 ) ? 1 : 0;
      anyLiveHeatFace = anyLiveHeatFace || ( stack.isEssentialHeatFace[i] == 0 && stack.isCondensedHeatFace[i] == 0 );

      // prescribed flux: sigma q_f = g |f|, d( sigma q_f )/dq_f = 0
      real64 const heatFlux = stack.isEssentialHeatFace[i] == 1 ? m_heatBc[kf].essentialFlux( stack.orientation[i], m_faceArea[kf] )
                                                                : m_faceHeatFlux[kf];
      stack.localHeatFlux[i] = stack.orientation[i] * heatFlux;
      stack.dLocalHeatFlux[i] = stack.isEssentialHeatFace[i] == 1 ? 0.0 : stack.orientation[i];
    }

    // M_e, needed if one face of the cell is non-condensed
    if( anyLiveHeatFace )
    {
      real64 const cond[3] = { m_effCond[ei][0][0], m_effCond[ei][0][1], m_effCond[ei][0][2] };
      real64 const chi = static_cast< real64 >( m_mfdFlag[ei] );
      mimeticInnerProduct::AdaptiveInnerProduct< IP >::template computeM< NUM_FACE >( m_nodePosition,
                                                                                      m_faceToNodes,
                                                                                      m_elemToFaces[ei],
                                                                                      m_elemCenter[ei],
                                                                                      m_elemVolume[ei],
                                                                                      cond,
                                                                                      m_lengthTolerance,
                                                                                      chi,
                                                                                      stack.heatMatrix );
    }

    // Fourier law on the non-condensed faces: sum_j (M_e)_ij sigma_j q_j - T_K + T_f = 0.
    // T_f cancels between the two cells of an interior face; the boundary condition defines it on a boundary face
    for( integer i = 0; i < NUM_FACE; ++i )
    {
      if( stack.isEssentialHeatFace[i] == 1 || stack.isCondensedHeatFace[i] == 1 )
      {
        continue;
      }
      localIndex const kf = m_elemToFaces[ei][i];
      real64 res = -m_elemTemp[ei];
      for( integer j = 0; j < NUM_FACE; ++j )
      {
        res += stack.heatMatrix( i, j ) * stack.localHeatFlux[j];
        stack.dHeatResidual_dHeatFlux[i][j] = stack.heatMatrix( i, j ) * stack.dLocalHeatFlux[j];
      }
      mixedMimeticBoundary::FaceBoundaryCondition const bc = m_heatBc[kf];
      if( bc.type != mixedMimeticBoundary::BoundaryType::interior )
      {
        real64 dTrace_dFlux = 0.0;
        res += bc.trace( stack.orientation[i], m_faceArea[kf], stack.localHeatFlux[i], dTrace_dFlux );
        stack.dHeatResidual_dHeatFlux[i][i] += dTrace_dFlux * stack.orientation[i];
      }
      stack.heatResidual[i] = res;
    }

    // dt sum_f sigma_f ( m_f h_f^up + q_f ) over the non-condensed faces of each operator
    stack.energyFlux = 0.0;
    stack.dEnergyFlux_dOwnH = 0.0;
    for( integer i = 0; i < NUM_FACE; ++i )
    {
      localIndex const kf = m_elemToFaces[ei][i];
      stack.dEnergyFlux_dFlux[i] = 0.0;
      stack.dEnergyFlux_dHeatFlux[i] = 0.0;
      stack.dEnergyFlux_dNbrH[i] = 0.0;

      if( stack.isCondensedFace[i] == 0 )
      {
        real64 const outflow = stack.localFlux[i];
        real64 hUp = m_elemEnthalpy[ei];
        bool upwindIsOwn = true;
        if( outflow < 0.0 )
        {
          bool const onBoundary = ( m_elemRegionList[kf][0] < 0 || m_elemRegionList[kf][1] < 0 );
          if( onBoundary )
          {
            if( m_isEnthalpyBcFace[kf] == 1 )
            {
              hUp = m_bcEnthalpy[kf];
              upwindIsOwn = false;
            }
          }
          else
          {
            // upstream cell: the other cell of the face
            integer const k = ( m_elemLocalToGlobalAcc[m_elemRegionList[kf][0]][m_elemSubRegionList[kf][0]][m_elemList[kf][0]]
                                == m_myElemLocalToGlobal[ei] ) ? 1 : 0;
            localIndex const er = m_elemRegionList[kf][k];
            localIndex const esr = m_elemSubRegionList[kf][k];
            localIndex const en = m_elemList[kf][k];
            hUp = m_enthalpyAcc[er][esr][en];
            stack.nbrHDofColIndex[i] = m_elemDofNumberAcc[er][esr][en] + CellDof::enthalpy;
            stack.dEnergyFlux_dNbrH[i] = m_dt * outflow;
            upwindIsOwn = false;
          }
        }
        stack.energyFlux += m_dt * outflow * hUp;
        stack.dEnergyFlux_dFlux[i] += m_dt * stack.dLocalFlux[i] * hUp;
        if( upwindIsOwn )
        {
          stack.dEnergyFlux_dOwnH += m_dt * outflow;
        }
      }

      if( stack.isCondensedHeatFace[i] == 0 )
      {
        stack.energyFlux += m_dt * stack.localHeatFlux[i];
        stack.dEnergyFlux_dHeatFlux[i] = m_dt * stack.dLocalHeatFlux[i];
      }
    }
  }

  GEOS_HOST_DEVICE
  void complete( localIndex const ei,
                 StackVariables & stack ) const
  {
    Base::complete( ei, stack );

    globalIndex const tempDofCol = stack.elemDofColIndex + CellDof::temperature;
    globalIndex const ownHDofCol = stack.elemDofColIndex + CellDof::enthalpy;

    // d( Darcy residual )/dT_K on the non-condensed faces
    for( integer i = 0; i < NUM_FACE; ++i )
    {
      localIndex const kf = m_elemToFaces[ei][i];
      if( m_faceGhostRank[kf] >= 0 || stack.isCondensedFace[i] == 1 || stack.isEssentialFace[i] == 1 )
      {
        continue;
      }
      real64 const dRes_dTemp = stack.orientation[i] * stack.dFaceResidual_dTemp[i];
      m_localMatrix.template addToRow< parallelDeviceAtomic >( stack.faceCenteredEqnRowIndex[i], &tempDofCol, &dRes_dTemp, 1 );
    }

    // energy balance of the cell
    if( m_elemGhostRank[ei] < 0 )
    {
      localIndex const row = stack.cellCenteredEqnRowIndex + CellDof::temperature;
      m_localRhs[row] += stack.energyFlux;
      m_localMatrix.template addToRowBinarySearchUnsorted< serialAtomic >( row, &stack.faceDofColIndices[0],
                                                                           &stack.dEnergyFlux_dFlux[0], NUM_FACE );
      m_localMatrix.template addToRowBinarySearchUnsorted< serialAtomic >( row, &stack.faceHeatDofColIndices[0],
                                                                           &stack.dEnergyFlux_dHeatFlux[0], NUM_FACE );
      m_localMatrix.template addToRow< serialAtomic >( row, &ownHDofCol, &stack.dEnergyFlux_dOwnH, 1 );
      for( integer i = 0; i < NUM_FACE; ++i )
      {
        if( stack.nbrHDofColIndex[i] >= 0 )
        {
          m_localMatrix.template addToRow< serialAtomic >( row, &stack.nbrHDofColIndex[i], &stack.dEnergyFlux_dNbrH[i], 1 );
        }
      }
    }

    // Fourier law, multiplied by sigma_f: the contributions of the two cells add up to the face equation
    for( integer i = 0; i < NUM_FACE; ++i )
    {
      localIndex const kf = m_elemToFaces[ei][i];
      if( m_faceGhostRank[kf] >= 0 || stack.isCondensedHeatFace[i] == 1 )
      {
        continue;
      }
      localIndex const row = stack.faceHeatEqnRowIndex[i];
      real64 const sigma = stack.orientation[i];

      if( stack.isEssentialHeatFace[i] == 1 )
      {
        // Neumann condition: q_f - sigma g |f| = 0
        real64 const one = 1.0;
        real64 const target = m_heatBc[kf].essentialFlux( sigma, m_faceArea[kf] );
        RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[row], m_faceHeatFlux[kf] - target );
        m_localMatrix.template addToRow< parallelDeviceAtomic >( row, &stack.faceHeatDofColIndices[i], &one, 1 );
        continue;
      }

      RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[row], sigma * stack.heatResidual[i] );
      real64 const dRes_dTemp = -sigma;
      m_localMatrix.template addToRow< parallelDeviceAtomic >( row, &tempDofCol, &dRes_dTemp, 1 );
      real64 dRes_dHeatFlux[NUM_FACE]{};
      for( integer j = 0; j < NUM_FACE; ++j )
      {
        dRes_dHeatFlux[j] = sigma * stack.dHeatResidual_dHeatFlux[i][j];
      }
      m_localMatrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( row, &stack.faceHeatDofColIndices[0],
                                                                                   &dRes_dHeatFlux[0], NUM_FACE );
    }
  }

protected:

  arrayView1d< globalIndex const > const m_faceHeatDofNumber;
  arrayView1d< real64 const > const m_faceArea;
  arrayView1d< real64 const > const m_faceHeatFlux;
  mixedMimeticBoundary::FaceBoundaryView const m_heatBc;
  arrayView1d< integer const > const m_isEnthalpyBcFace;
  arrayView1d< real64 const > const m_bcEnthalpy;
  arrayView1d< real64 const > const m_elemTemp;
  arrayView1d< real64 const > const m_elemEnthalpy;
  arrayView3d< real64 const > const m_effCond;
  SortedArrayView< localIndex const > const m_regionFilter;
  ElementViewConst< arrayView1d< real64 const > > const m_enthalpyAcc;
  ElementViewConst< arrayView1d< globalIndex const > > const m_elemDofNumberAcc;
  ElementViewConst< arrayView1d< globalIndex const > > const m_elemLocalToGlobalAcc;
};

/**
 * @class ElementBasedAssemblyKernelFactory
 */
class ElementBasedAssemblyKernelFactory
{
public:

  template< typename POLICY >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   real64 const lengthTolerance,
                   string const elemDofKey,
                   string const faceDofKey,
                   string const faceHeatDofKey,
                   NodeManager const & nodeManager,
                   FaceManager const & faceManager,
                   ElementRegionManager const & elemManager,
                   CellElementSubRegion const & subRegion,
                   mimeticInnerProduct::MimeticInnerProductBase const & mimeticInnerProductBase,
                   constitutive::SingleFluidBase const & fluid,
                   constitutive::PermeabilityBase const & permeability,
                   constitutive::SinglePhaseThermalConductivityBase const & conductivity,
                   SortedArrayView< localIndex const > const & regionFilter,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const enthalpyAcc =
      elemManager.constructArrayViewAccessor< real64, 1 >( fields::mixedMimetic::enthalpy::key() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const elemDofNumberAcc =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const elemLocalToGlobalAcc =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( ObjectManagerBase::viewKeyStruct::localToGlobalMapString() );

    mixedMimeticInnerProductDispatch( mimeticInnerProductBase,
                                      [&] ( auto const mimeticInnerProduct )
    {
      using IP = TYPEOFREF( mimeticInnerProduct );

      mixedMimeticKernels::internal::kernelLaunchSelectorFaceSwitch( subRegion.numFacesPerElement(), [&] ( auto NUM_FACES )
      {
        ElementBasedAssemblyKernel< NUM_FACES, IP >
        kernel( rankOffset, lengthTolerance, elemDofKey, faceDofKey, faceHeatDofKey, nodeManager, faceManager, elemManager,
                subRegion, fluid, permeability, conductivity, enthalpyAcc, elemDofNumberAcc, elemLocalToGlobalAcc,
                regionFilter, dt, localMatrix, localRhs );
        ElementBasedAssemblyKernel< NUM_FACES, IP >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    } );
  }
};

/******************************** TpfaCondensedFluxKernel ********************************/

/**
 * @class TpfaCondensedFluxKernel
 * @brief Contributions of the condensed faces, with L and R the cells of the face and Phi = p - rho ( gamma_K - gamma_f )
 *        the potential: two-point fluxes m_f = ( Phi_L - Phi_R ) / A and q_f = ( T_L - T_R ) / A_q, terms
 *        +/- dt m_f and +/- dt ( m_f h_f^up + q_f ) in the balances of L and R, and the face equations
 *        A m_f = Phi_L - Phi_R, A_q q_f = T_L - T_R.
 */
class TpfaCondensedFluxKernel : public singlePhaseMixedMFDKernels::TpfaCondensedFluxKernel
{
public:

  using Base = singlePhaseMixedMFDKernels::TpfaCondensedFluxKernel;
  using ThermalDerivOffset = constitutive::singlefluid::DerivativeOffsetC< 1 >;
  using TwoPoint = Base::TwoPointFlux< 2 >;

  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using ConductivityAccessors =
    StencilMaterialAccessors< constitutive::SinglePhaseThermalConductivityBase,
                              fields::thermalconductivity::effectiveConductivity >;

  TpfaCondensedFluxKernel( globalIndex const rankOffset,
                           real64 const & lengthTolerance,
                           string const faceDofKey,
                           string const faceHeatDofKey,
                           NodeManager const & nodeManager,
                           FaceManager const & faceManager,
                           Base::DofNumberAccessor const & elemDofNumber,
                           Base::DofNumberAccessor const & elemLocalToGlobal,
                           Base::GhostRankAccessor const & elemGhostRank,
                           ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const & elemCenter,
                           Base::FlowAccessors const & flowAccessors,
                           Base::FluidAccessors const & fluidAccessors,
                           Base::PermeabilityAccessors const & permAccessors,
                           ConductivityAccessors const & condAccessors,
                           ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const & temperature,
                           ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const & enthalpy,
                           SortedArrayView< localIndex const > const & regionFilter,
                           real64 const & dt,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs )
    : Base( rankOffset, lengthTolerance, faceDofKey, nodeManager, faceManager, elemDofNumber, elemLocalToGlobal,
            elemGhostRank, elemCenter, flowAccessors, fluidAccessors, permAccessors, regionFilter, dt, localMatrix, localRhs ),
    m_faceHeatDofNumber( faceManager.getReference< array1d< globalIndex > >( faceHeatDofKey ) ),
    m_faceArea( faceManager.faceArea() ),
    m_faceHeatFlux( faceManager.getField< fields::mixedMimetic::faceHeatFlux >() ),
    m_heatBc{ faceManager.getField< fields::mixedMimetic::heatBoundaryType >(),
             faceManager.getField< fields::mixedMimetic::heatBoundaryValue >(),
             faceManager.getField< fields::mixedMimetic::heatBoundaryCoefficient >() },
    m_isEnthalpyBcFace( faceManager.getField< fields::mixedMimetic::isEnthalpyBcFace >() ),
    m_bcEnthalpy( faceManager.getField< fields::mixedMimetic::bcEnthalpy >() ),
    m_effCond( condAccessors.get( fields::thermalconductivity::effectiveConductivity {} ) ),
    m_temp( temperature.toNestedViewConst() ),
    m_enthalpy( enthalpy.toNestedViewConst() )
  {}

  GEOS_HOST_DEVICE
  void compute( localIndex const kf ) const
  {
    if( m_faceStencilLabel[kf] != 0 )
    {
      return;
    }

    localIndex er[2]{}, esr[2]{}, ei[2]{};
    integer const numElems = gatherCells( kf, er, esr, ei );
    if( numElems == 0 )
    {
      return;
    }

    // an operator is condensed on the face unless its flux is prescribed
    bool const flowCondensed = !m_flowBc[kf].isEssential();
    bool const heatCondensed = !m_heatBc[kf].isEssential();

    globalIndex elemDofCols[2][CellDof::num];
    localIndex cellRow[2];
    bool cellOwned[2];
    for( integer k = 0; k < numElems; ++k )
    {
      for( integer d = 0; d < CellDof::num; ++d )
      {
        elemDofCols[k][d] = m_elemDofNumber[er[k]][esr[k]][ei[k]] + d;
      }
      cellRow[k] = m_elemDofNumber[er[k]][esr[k]][ei[k]] - m_rankOffset;
      cellOwned[k] = m_elemGhostRank[er[k]][esr[k]][ei[k]] < 0;
    }
    // columns ( p_L, T_L, p_R, T_R )
    globalIndex ptCols[4];
    for( integer k = 0; k < numElems; ++k )
    {
      ptCols[2 * k] = elemDofCols[k][CellDof::pressure];
      ptCols[2 * k + 1] = elemDofCols[k][CellDof::temperature];
    }

    // mass flux m_f( p_L, T_L, p_R, T_R )
    if( flowCondensed )
    {
      TwoPoint tp;
      computeTwoPointMassFlux< 2 >( kf, er, esr, ei, numElems, tp );

      // h_f^up = h_L if m_f >= 0, else h_R, or the prescribed value on an inflow boundary face
      real64 hUp;
      integer hUpCell = -1;
      if( tp.flux >= 0.0 )
      {
        hUp = m_enthalpy[er[0]][esr[0]][ei[0]];
        hUpCell = 0;
      }
      else if( numElems == 2 )
      {
        hUp = m_enthalpy[er[1]][esr[1]][ei[1]];
        hUpCell = 1;
      }
      else
      {
        hUp = ( m_isEnthalpyBcFace[kf] == 1 ) ? m_bcEnthalpy[kf] : m_enthalpy[er[0]][esr[0]][ei[0]];
        hUpCell = ( m_isEnthalpyBcFace[kf] == 1 ) ? -1 : 0;
      }

      // balances of L and R: +/- dt m_f (mass), +/- dt m_f h_f^up (energy)
      for( integer k = 0; k < numElems; ++k )
      {
        if( !cellOwned[k] )
        {
          continue;
        }
        real64 const sign = ( k == 0 ) ? 1.0 : -1.0;
        real64 massVals[4], energyVals[4];
        for( integer j = 0; j < numElems; ++j )
        {
          for( integer d = 0; d < 2; ++d )
          {
            massVals[2 * j + d] = sign * m_dt * tp.dFlux[j][d];
            energyVals[2 * j + d] = sign * m_dt * tp.dFlux[j][d] * hUp;
          }
        }
        localIndex const massRow = cellRow[k] + CellDof::pressure;
        RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[massRow], sign * m_dt * tp.flux );
        m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( massRow, ptCols, massVals, 2 * numElems );

        localIndex const energyRow = cellRow[k] + CellDof::temperature;
        RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[energyRow], sign * m_dt * tp.flux * hUp );
        m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( energyRow, ptCols, energyVals, 2 * numElems );
        if( hUpCell >= 0 )
        {
          real64 const dh = sign * m_dt * tp.flux;
          m_localMatrix.addToRow< parallelDeviceAtomic >( energyRow, &elemDofCols[hUpCell][CellDof::enthalpy], &dh, 1 );
        }
      }

      // face equation: A m_f - ( Phi_L - Phi_R ) = 0
      if( m_faceGhostRank[kf] < 0 )
      {
        localIndex const row = m_faceDofNumber[kf] - m_rankOffset;
        globalIndex const faceDofCol = m_faceDofNumber[kf];
        real64 const mCurrent = m_faceFlux[kf];
        m_localRhs[row] += tp.A * mCurrent - tp.pot;
        m_localMatrix.addToRow< serialAtomic >( row, &faceDofCol, &tp.A, 1 );
        real64 vals[4];
        for( integer k = 0; k < numElems; ++k )
        {
          for( integer d = 0; d < 2; ++d )
          {
            vals[2 * k + d] = tp.dA[k][d] * mCurrent - tp.dPot[k][d];
          }
        }
        m_localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( row, ptCols, vals, 2 * numElems );
      }
    }

    // heat flux: A_q q_f = T_L - T_R, A_q = sum_K 1 / t_K; on a boundary face T_R = g and A_q includes 1 / ( alpha |f| )
    if( heatCondensed )
    {
      real64 Aq = 0.0;
      for( integer k = 0; k < numElems; ++k )
      {
        real64 const cond[3] = { m_effCond[er[k]][esr[k]][ei[k]][0][0],
                                 m_effCond[er[k]][esr[k]][ei[k]][0][1],
                                 m_effCond[er[k]][esr[k]][ei[k]][0][2] };
        Aq += 1.0 / oneSidedTrans( kf, er[k], esr[k], ei[k], cond );
      }
      real64 const TL = m_temp[er[0]][esr[0]][ei[0]];
      real64 TR;
      if( numElems == 2 )
      {
        TR = m_temp[er[1]][esr[1]][ei[1]];
      }
      else
      {
        TR = m_heatBc[kf].value;
        Aq += m_heatBc[kf].resistance( m_faceArea[kf] );
      }
      real64 const potq = TL - TR;
      real64 const qFlux = potq / Aq;
      real64 const dq_dT[2] = { 1.0 / Aq, -1.0 / Aq };
      globalIndex tCols[2] = { elemDofCols[0][CellDof::temperature], elemDofCols[1][CellDof::temperature] };

      // energy balances of L and R: +/- dt q_f
      for( integer k = 0; k < numElems; ++k )
      {
        if( !cellOwned[k] )
        {
          continue;
        }
        real64 const sign = ( k == 0 ) ? 1.0 : -1.0;
        localIndex const energyRow = cellRow[k] + CellDof::temperature;
        RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[energyRow], sign * m_dt * qFlux );
        real64 vals[2] = { sign * m_dt * dq_dT[0], sign * m_dt * dq_dT[1] };
        m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( energyRow, tCols, vals, numElems );
      }

      // face equation: A_q q_f - ( T_L - T_R ) = 0
      if( m_faceGhostRank[kf] < 0 )
      {
        localIndex const row = m_faceHeatDofNumber[kf] - m_rankOffset;
        globalIndex const faceHeatDofCol = m_faceHeatDofNumber[kf];
        m_localRhs[row] += Aq * m_faceHeatFlux[kf] - potq;
        m_localMatrix.addToRow< serialAtomic >( row, &faceHeatDofCol, &Aq, 1 );
        real64 vals[2] = { -1.0, 1.0 };
        m_localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( row, tCols, vals, numElems );
      }
    }
  }

  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numFaces,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;
    forAll< POLICY >( numFaces, [=] GEOS_HOST_DEVICE ( localIndex const kf )
    {
      kernelComponent.compute( kf );
    } );
  }

protected:

  arrayView1d< globalIndex const > const m_faceHeatDofNumber;
  arrayView1d< real64 const > const m_faceArea;
  arrayView1d< real64 const > const m_faceHeatFlux;
  mixedMimeticBoundary::FaceBoundaryView const m_heatBc;
  arrayView1d< integer const > const m_isEnthalpyBcFace;
  arrayView1d< real64 const > const m_bcEnthalpy;
  ElementViewConst< arrayView3d< real64 const > > const m_effCond;
  ElementViewConst< arrayView1d< real64 const > > const m_temp;
  ElementViewConst< arrayView1d< real64 const > > const m_enthalpy;
};

/**
 * @class TpfaCondensedFluxKernelFactory
 */
class TpfaCondensedFluxKernelFactory
{
public:

  template< typename POLICY >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   real64 const lengthTolerance,
                   string const elemDofKey,
                   string const faceDofKey,
                   string const faceHeatDofKey,
                   string const solverName,
                   NodeManager const & nodeManager,
                   FaceManager const & faceManager,
                   ElementRegionManager const & elemManager,
                   SortedArrayView< localIndex const > const & regionFilter,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    using Base = singlePhaseMixedMFDKernels::TpfaCondensedFluxKernel;
    Base::DofNumberAccessor elemDofNumber = elemManager.constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
    elemDofNumber.setName( solverName + "/accessors/" + elemDofKey );
    Base::DofNumberAccessor const elemLocalToGlobal =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( ObjectManagerBase::viewKeyStruct::localToGlobalMapString() );
    Base::GhostRankAccessor const elemGhostRank =
      elemManager.constructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
    ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const elemCenter =
      elemManager.constructViewAccessor< array2d< real64 >, arrayView2d< real64 const > >( ElementSubRegionBase::viewKeyStruct::elementCenterString() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const temperature =
      elemManager.constructArrayViewAccessor< real64, 1 >( fields::flow::temperature::key() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const enthalpy =
      elemManager.constructArrayViewAccessor< real64, 1 >( fields::mixedMimetic::enthalpy::key() );
    Base::FlowAccessors flowAccessors( elemManager, solverName );
    Base::FluidAccessors fluidAccessors( elemManager, solverName );
    Base::PermeabilityAccessors permAccessors( elemManager, solverName );
    TpfaCondensedFluxKernel::ConductivityAccessors condAccessors( elemManager, solverName );

    TpfaCondensedFluxKernel kernel( rankOffset, lengthTolerance, faceDofKey, faceHeatDofKey,
                                    nodeManager, faceManager,
                                    elemDofNumber, elemLocalToGlobal, elemGhostRank, elemCenter,
                                    flowAccessors, fluidAccessors, permAccessors, condAccessors,
                                    temperature, enthalpy,
                                    regionFilter, dt, localMatrix, localRhs );
    TpfaCondensedFluxKernel::launch< POLICY >( faceManager.size(), kernel );
  }
};

/******************************** EnthalpyClosureKernel ********************************/

/**
 * @class EnthalpyClosureKernel
 * @brief Closure h_K - h_eos( p_K, T_K ) = 0.
 */
class EnthalpyClosureKernel
{
public:

  using ThermalDerivOffset = constitutive::singlefluid::DerivativeOffsetC< 1 >;

  template< typename POLICY >
  static void
  launch( globalIndex const rankOffset,
          string const elemDofKey,
          ElementSubRegionBase const & subRegion,
          constitutive::SingleFluidBase const & fluid,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
  {
    arrayView1d< globalIndex const > const elemDofNumber = subRegion.getReference< array1d< globalIndex > >( elemDofKey );
    arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const enthalpy = subRegion.getField< fields::mixedMimetic::enthalpy >();
    arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const eosEnthalpy = fluid.enthalpy();
    arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const dEosEnthalpy = fluid.dEnthalpy();

    forAll< POLICY >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( elemGhostRank[ei] >= 0 )
      {
        return;
      }
      localIndex const row = elemDofNumber[ei] - rankOffset + CellDof::enthalpy;
      globalIndex const cols[3] = { elemDofNumber[ei] + CellDof::pressure,
                                    elemDofNumber[ei] + CellDof::temperature,
                                    elemDofNumber[ei] + CellDof::enthalpy };
      real64 const vals[3] = { -dEosEnthalpy[ei][0][ThermalDerivOffset::dP],
                               -dEosEnthalpy[ei][0][ThermalDerivOffset::dT],
                               1.0 };
      localRhs[row] += enthalpy[ei] - eosEnthalpy[ei][0];
      localMatrix.template addToRow< serialAtomic >( row, cols, vals, 3 );
    } );
  }
};

} // namespace thermalSinglePhaseMixedMFDKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMALMIXEDMFDKERNELS_HPP
