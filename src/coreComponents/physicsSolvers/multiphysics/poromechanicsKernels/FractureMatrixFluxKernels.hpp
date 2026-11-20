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
 * @file FractureMatrixFluxKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_FRACTUREMATRIXFLUXKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_FRACTUREMATRIXFLUXKERNELS_HPP

#include "constitutive/permeability/PermeabilityFields.hpp"
#include "finiteVolume/FaceElementToCellStencil.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/FluxComputeKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/ThermalFluxComputeKernel.hpp"

namespace geos
{

namespace fractureMatrixFluxKernels
{

/**
 * @class FractureMatrixFluxKernel
 * @brief Single-phase flux between the conforming fracture elements and the matrix cells (FaceElementToCellStencil),
 *        which also computes the derivative of the flux w.r.t. the normal jump of the fracture element.
 *
 * The fracture-side half-transmissibility t1 = k_f,n( Ha ) * 2 A / Ha depends on the hydraulic aperture Ha
 * (see FaceElementToCellStencilWrapper::addHydraulicApertureContribution). The fluxes are linear in their
 * transmissibility T = m t0 t1 / ( t0 + t1 ), hence
 *   dFlux/dNormalJump = flux * dln(T)/dln(t1) * dln(t1)/dNormalJump.
 * The mass flux and the advective (enthalpy-carried) energy flux follow the hydraulic transmissibility, the
 * conductive energy flux the thermal one, whose fracture side depends on the aperture through the weight only.
 *
 * @tparam FLUX_KERNEL the single-phase flux kernel on a FaceElementToCellStencilWrapper
 *         (singlePhaseFVMKernels::FluxComputeKernel or thermalSinglePhaseFVMKernels::FluxComputeKernel)
 */
template< typename FLUX_KERNEL >
class FractureMatrixFluxKernel : public FLUX_KERNEL
{
public:

  using Base = FLUX_KERNEL;
  using typename Base::StackVariables;
  using Base::numEqn;
  using Base::m_dt;
  using Base::m_permeability;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;

  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using FracturePermeabilityAccessors = StencilMaterialAccessors< constitutive::PermeabilityBase,
                                                                  fields::permeability::dPerm_dDispJump >;
  using FractureApertureAccessors = StencilAccessors< fields::flow::hydraulicAperture,
                                                      fields::flow::dHydraulicAperture_dNormalJump >;

  /// Index, in the connection, of the fracture element (the cell is 0, see FaceElementToCellStencil)
  static constexpr localIndex fractureIndex = 1;

  /// Index, in the connection, of the matrix cell
  static constexpr localIndex cellIndex = 0;

  /**
   * @brief Constructor
   * @param[in] fracturePermeabilityAccessors accessors to the derivative of the fracture permeability wrt the jump
   * @param[in] fractureApertureAccessors accessors to the hydraulic aperture and its derivative wrt the normal jump
   * @param[in] dFlux_dNormalJump derivative, for each connection and equation, of the flux of the cell row wrt the
   *            normal jump of the fracture element (the flux of the fracture row is its opposite)
   * @param[in] baseArgs the arguments of the constructor of FLUX_KERNEL
   */
  template< typename ... BASE_ARGS >
  FractureMatrixFluxKernel( FracturePermeabilityAccessors const & fracturePermeabilityAccessors,
                            FractureApertureAccessors const & fractureApertureAccessors,
                            arrayView2d< real64 > const & dFlux_dNormalJump,
                            BASE_ARGS && ... baseArgs )
    : Base( std::forward< BASE_ARGS >( baseArgs )... ),
    m_dPerm_dDispJump( fracturePermeabilityAccessors.get( fields::permeability::dPerm_dDispJump {} ) ),
    m_hydraulicAperture( fractureApertureAccessors.get( fields::flow::hydraulicAperture {} ) ),
    m_dHydraulicAperture_dNormalJump( fractureApertureAccessors.get( fields::flow::dHydraulicAperture_dNormalJump {} ) ),
    m_dFlux_dNormalJump( dFlux_dNormalJump )
  {}

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian, and their derivatives wrt the normal jump
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack ) const
  {
    Base::computeFlux( iconn, stack );

    localIndex const er  = m_seri( iconn, fractureIndex );
    localIndex const esr = m_sesri( iconn, fractureIndex );
    localIndex const ei  = m_sei( iconn, fractureIndex );

    // dln(t1)/dNormalJump: the stencil weight is inversely proportional to the hydraulic aperture,
    // and the (normal) permeability may depend on it
    real64 const hydraulicAperture = m_hydraulicAperture[er][esr][ei];
    real64 const dLogAperture_dNormalJump = ( hydraulicAperture > 0.0 )
      ? m_dHydraulicAperture_dNormalJump[er][esr][ei] / hydraulicAperture
      : 0.0;
    real64 const perm = m_permeability[er][esr][ei][0][2];
    real64 const dLogPerm_dNormalJump = ( perm > 0.0 ) ? m_dPerm_dDispJump[er][esr][ei][0][2][0] / perm : 0.0;

    real64 trans = 0.0;
    real64 dLogTrans_dLogFractureTrans = 0.0;
    m_stencilWrapper.computeWeightAndFractureLogDerivative( iconn, m_permeability, trans, dLogTrans_dLogFractureTrans );
    real64 const dLogTrans_dNormalJump = dLogTrans_dLogFractureTrans * ( dLogPerm_dNormalJump - dLogAperture_dNormalJump );

    // stack.localFlux holds dt * flux on the cell row (k = 0), and its opposite on the fracture row
    m_dFlux_dNormalJump( iconn, 0 ) = stack.localFlux[cellIndex * numEqn] * dLogTrans_dNormalJump;

    if constexpr ( numEqn > 1 )
    {
      // energy = advective part (proportional to the hydraulic transmissibility) + conductive part
      integer constexpr energyEqnIndex = numEqn - 1;

      real64 thermalTrans = 0.0;
      real64 dLogThermalTrans_dLogFractureTrans = 0.0;
      m_stencilWrapper.computeWeightAndFractureLogDerivative( iconn, this->m_thermalConductivity, thermalTrans, dLogThermalTrans_dLogFractureTrans );

      localIndex const erc  = m_seri( iconn, cellIndex );
      localIndex const esrc = m_sesri( iconn, cellIndex );
      localIndex const eic  = m_sei( iconn, cellIndex );
      real64 const conductiveFlux = m_dt * thermalTrans * ( this->m_temp[erc][esrc][eic] - this->m_temp[er][esr][ei] );
      real64 const advectiveFlux = stack.localFlux[cellIndex * numEqn + energyEqnIndex] - conductiveFlux;

      // the conductivity does not depend on the aperture: only the weight does
      m_dFlux_dNormalJump( iconn, energyEqnIndex ) = advectiveFlux * dLogTrans_dNormalJump
                                                     - conductiveFlux * dLogThermalTrans_dLogFractureTrans * dLogAperture_dNormalJump;
    }
  }

protected:

  /// View on the derivative of the fracture permeability wrt the displacement jump
  ElementViewConst< arrayView4d< real64 const > > const m_dPerm_dDispJump;

  /// Views on the hydraulic aperture and its derivative wrt the normal jump
  ElementViewConst< arrayView1d< real64 const > > const m_hydraulicAperture;
  ElementViewConst< arrayView1d< real64 const > > const m_dHydraulicAperture_dNormalJump;

  /// Derivative, for each connection and equation, of the flux of the cell row wrt the normal jump
  arrayView2d< real64 > const m_dFlux_dNormalJump;
};

/**
 * @class FractureMatrixFluxKernelFactory
 * @tparam IS_THERMAL whether the flow is thermal (mass and energy balances) or not (mass balance only)
 */
template< bool IS_THERMAL >
struct FractureMatrixFluxKernelFactory
{
  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[out] dFlux_dNormalJump derivative, for each connection and equation, of the flux of the cell row
   *             wrt the normal jump of the fracture element
   */
  template< typename POLICY >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   string const & dofKey,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   FaceElementToCellStencilWrapper const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs,
                   arrayView2d< real64 > const & dFlux_dNormalJump )
  {
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

    if constexpr ( IS_THERMAL )
    {
      using FluxKernelType = thermalSinglePhaseFVMKernels::FluxComputeKernel< 2, 2, FaceElementToCellStencilWrapper >;
      using KernelType = FractureMatrixFluxKernel< FluxKernelType >;

      typename KernelType::SinglePhaseFlowAccessors flowAccessors( elemManager, solverName );
      typename KernelType::ThermalSinglePhaseFlowAccessors thermalFlowAccessors( elemManager, solverName );
      typename KernelType::SinglePhaseFluidAccessors fluidAccessors( elemManager, solverName );
      typename KernelType::ThermalSinglePhaseFluidAccessors thermalFluidAccessors( elemManager, solverName );
      typename KernelType::PermeabilityAccessors permAccessors( elemManager, solverName );
      typename KernelType::ThermalConductivityAccessors thermalConductivityAccessors( elemManager, solverName );
      typename KernelType::FracturePermeabilityAccessors fracturePermAccessors( elemManager, solverName );
      typename KernelType::FractureApertureAccessors fractureApertureAccessors( elemManager, solverName );

      KernelType kernel( fracturePermAccessors, fractureApertureAccessors, dFlux_dNormalJump,
                         rankOffset, stencilWrapper, dofNumberAccessor,
                         flowAccessors, thermalFlowAccessors, fluidAccessors, thermalFluidAccessors,
                         permAccessors, thermalConductivityAccessors,
                         dt, localMatrix, localRhs );
      KernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
    }
    else
    {
      using FluxKernelType = singlePhaseFVMKernels::FluxComputeKernel< 1, 1, FaceElementToCellStencilWrapper >;
      using KernelType = FractureMatrixFluxKernel< FluxKernelType >;

      typename KernelType::SinglePhaseFlowAccessors flowAccessors( elemManager, solverName );
      typename KernelType::SinglePhaseFluidAccessors fluidAccessors( elemManager, solverName );
      typename KernelType::PermeabilityAccessors permAccessors( elemManager, solverName );
      typename KernelType::FracturePermeabilityAccessors fracturePermAccessors( elemManager, solverName );
      typename KernelType::FractureApertureAccessors fractureApertureAccessors( elemManager, solverName );

      KernelType kernel( fracturePermAccessors, fractureApertureAccessors, dFlux_dNormalJump,
                         rankOffset, stencilWrapper, dofNumberAccessor,
                         flowAccessors, fluidAccessors, permAccessors,
                         dt, localMatrix, localRhs );
      KernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
    }
  }
};

} // namespace fractureMatrixFluxKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_FRACTUREMATRIXFLUXKERNELS_HPP
