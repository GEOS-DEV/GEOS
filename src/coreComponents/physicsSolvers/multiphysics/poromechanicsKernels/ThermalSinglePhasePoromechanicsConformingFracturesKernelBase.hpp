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
 * @file ThermalSinglePhasePoromechanicsConformingFracturesKernelBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICSCONFORMINGFRACTURESKERNELBASE_HPP
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICSCONFORMINGFRACTURESKERNELBASE_HPP

#include "physicsSolvers/fluidFlow/kernels/singlePhase/FluxComputeKernelBase.hpp"

namespace geos
{

namespace thermalConformingFracturesKernels
{

/**
 * @class ThermalConformingFracturesConnectorBasedAssemblyKernel
 * @tparam BASE_KERNEL the (Lagrange-multiplier- or Augmented-Lagrangian-Multiplier-formulation) isothermal
 *   ConnectorBasedAssemblyKernel this thermal layer extends. The thermal logic (enthalpy-carried advective
 *   flux, conductive flux, energy-balance row/Jacobian assembly, dR_dAper energy-block write) is identical
 *   across both contact formulations; only the underlying mass-flux/aperture-sensitivity kernel differs.
 */
template< integer NUM_EQN, integer NUM_DOF, template< integer, integer > class BASE_KERNEL >
class ThermalConformingFracturesConnectorBasedAssemblyKernel : public BASE_KERNEL< NUM_EQN, NUM_DOF >
{
public:

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using SinglePhaseFVMAbstractBase = singlePhaseFVMKernels::FluxComputeKernelBase;
  using DofNumberAccessor = SinglePhaseFVMAbstractBase::DofNumberAccessor;
  using SinglePhaseFlowAccessors = SinglePhaseFVMAbstractBase::SinglePhaseFlowAccessors;
  using SinglePhaseFluidAccessors = SinglePhaseFVMAbstractBase::SinglePhaseFluidAccessors;
  using PermeabilityAccessors = SinglePhaseFVMAbstractBase::PermeabilityAccessors;
  using SinglePhaseFVMAbstractBase::m_dt;
  using SinglePhaseFVMAbstractBase::m_rankOffset;
  using SinglePhaseFVMAbstractBase::m_dofNumber;
  using SinglePhaseFVMAbstractBase::m_gravCoef;
  using SinglePhaseFVMAbstractBase::m_mob;
  using SinglePhaseFVMAbstractBase::m_dens;
  using SinglePhaseFVMAbstractBase::m_dDens;
  using SinglePhaseFVMAbstractBase::m_dMob;

  using SinglePhaseFVMBase = singlePhaseFVMKernels::FluxComputeKernel< NUM_EQN, NUM_DOF, SurfaceElementStencilWrapper >;
  using SinglePhaseFVMBase::numDof;
  using SinglePhaseFVMBase::numEqn;
  using SinglePhaseFVMBase::maxNumElems;
  using SinglePhaseFVMBase::maxNumConns;
  using SinglePhaseFVMBase::maxStencilSize;
  using SinglePhaseFVMBase::m_stencilWrapper;
  using SinglePhaseFVMBase::m_seri;
  using SinglePhaseFVMBase::m_sesri;
  using SinglePhaseFVMBase::m_sei;
  using SinglePhaseFVMBase::m_ghostRank;

  using Base = BASE_KERNEL< NUM_EQN, NUM_DOF >;

  using ThermalSinglePhaseFlowAccessors =
    StencilAccessors< fields::flow::temperature >;

  using ThermalSinglePhaseFluidAccessors =
    StencilMaterialAccessors< constitutive::SingleFluidBase,
                              fields::singlefluid::enthalpy,
                              fields::singlefluid::dEnthalpy >;

  using ThermalConductivityAccessors =
    StencilMaterialAccessors< constitutive::SinglePhaseThermalConductivityBase,
                              fields::thermalconductivity::effectiveConductivity >;

  /// Fracture permeability derivative, hydraulic aperture and its derivative wrt the normal jump, held by the mass-flux base kernel
  using FracturePermeabilityAccessors = typename Base::FracturePermeabilityAccessors;
  using FractureApertureAccessors = typename Base::FractureApertureAccessors;
  using Base::m_dR_dAper;


  ThermalConformingFracturesConnectorBasedAssemblyKernel( globalIndex const rankOffset,
                                                          SurfaceElementStencilWrapper const & stencilWrapper,
                                                          DofNumberAccessor const & flowDofNumberAccessor,
                                                          SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                                                          ThermalSinglePhaseFlowAccessors const & thermalSinglePhaseFlowAccessors,
                                                          SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                                                          ThermalSinglePhaseFluidAccessors const & thermalSinglePhaseFluidAccessors,
                                                          PermeabilityAccessors const & permeabilityAccessors,
                                                          FracturePermeabilityAccessors const & edfmPermeabilityAccessors,
                                                          ThermalConductivityAccessors const & thermalConductivityAccessors,
                                                          FractureApertureAccessors const & fractureApertureAccessors,
                                                          real64 const & dt,
                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                          arrayView1d< real64 > const & localRhs,
                                                          CRSMatrixView< real64, localIndex const > const & dR_dAper,
                                                          localIndex const dR_dAperEnergyOffset )
    : Base( rankOffset,
            stencilWrapper,
            flowDofNumberAccessor,
            singlePhaseFlowAccessors,
            singlePhaseFluidAccessors,
            permeabilityAccessors,
            edfmPermeabilityAccessors,
            fractureApertureAccessors,
            dt,
            localMatrix,
            localRhs,
            dR_dAper ),
    m_temp( thermalSinglePhaseFlowAccessors.get( fields::flow::temperature {} ) ),
    m_enthalpy( thermalSinglePhaseFluidAccessors.get( fields::singlefluid::enthalpy {} ) ),
    m_dEnthalpy( thermalSinglePhaseFluidAccessors.get( fields::singlefluid::dEnthalpy {} ) ),
    m_thermalConductivity( thermalConductivityAccessors.get( fields::thermalconductivity::effectiveConductivity {} ) ),
    m_dR_dAperEnergyOffset( dR_dAperEnergyOffset )
  {}


  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    /**
     * @brief Constructor for the stack variables
     * @param[in] size size of the stencil for this connection
     * @param[in] numElems number of elements for this connection
     */
    GEOS_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : Base::StackVariables( size, numElems ),
      dEnergyFlux_dAperture( numElems, size )
    {}
    using SinglePhaseFVMBase::StackVariables::stencilSize;
    using SinglePhaseFVMBase::StackVariables::numFluxElems;
    using SinglePhaseFVMBase::StackVariables::transmissibility;
    using SinglePhaseFVMBase::StackVariables::dTrans_dPres;
    using SinglePhaseFVMBase::StackVariables::dofColIndices;
    using SinglePhaseFVMBase::StackVariables::localFlux;
    using SinglePhaseFVMBase::StackVariables::localFluxJacobian;
    using Base::StackVariables::dFlux_dAperture;
    using Base::StackVariables::dLogAperture_dNormalJump;

    // Thermal transmissibility (no derivative wrt conductivity: conductivity is treated explicitly)

    real64 thermalTransmissibility[maxNumConns][2]{};

    /// t_j * d(thermalTransmissibility)/d(t_j) of each connection, for every half-transmissibility
    /// t_j = weight_j * conductivity_j of the stencil (unsigned: derivative of entry [0]). Since the stencil weight is
    /// linear in the hydraulic aperture, multiplying by d(ln hydraulicAperture_j)/d(normalJump_j) gives the aperture
    /// derivative. At a fracture junction it also involves the elements the connection does not link.
    real64 dThermalTrans_dLogHalfTrans[maxNumConns][maxNumElems]{};

    // The energy flux and its p/T derivatives are per connection (local to computeFlux), so that
    // stencils with more than two elements (fracture junctions) are handled correctly.

    /// Derivatives of the (advective + conductive) energy flux wrt aperture (i.e. wrt the normal
    /// component of the displacement jump, same convention as dFlux_dAperture), same shape/convention as
    /// dFlux_dAperture.
    stackArray2d< real64, maxNumElems *maxStencilSize > dEnergyFlux_dAperture;

  };

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the flux
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack ) const

  {
    // Thermal transmissibilities of all the connections of this face (conductivity treated explicitly).
    // Their log-derivatives wrt the half-transmissibility of every element of the stencil (t_j = weight_j * conductivity_j)
    // are used below for the aperture derivative of the conductive flux.
    m_stencilWrapper.computeWeightsAndLogDerivatives( iconn,
                                                      m_thermalConductivity,
                                                      stack.thermalTransmissibility,
                                                      stack.dThermalTrans_dLogHalfTrans );

    // The base computeFlux computes the mass flux and its derivatives connection by connection, and calls the
    // lambda below for each connection (pair k[0], k[1] of the stencil), with the quantities of that connection
    // only (potGrad, upwinding, dFlux_dAper, ...). The energy flux of the connection (advective + conductive)
    // and its derivatives are computed and assembled there, so that junctions (> 2 elements) are handled.
    Base::computeFlux( iconn, stack, [&] ( localIndex const (&k)[2],
                                           localIndex const (&seri)[2],
                                           localIndex const (&sesri)[2],
                                           localIndex const (&sei)[2],
                                           localIndex const,
                                           localIndex const connectionIndex,
                                           real64 const & alpha,
                                           real64 const & mobility,
                                           real64 const & potGrad,
                                           real64 const & massFlux,
                                           real64 const & dMassFlux_dTrans,
                                           real64 const (&dMassFlux_dP)[2],
                                           real64 const (&dMassFlux_dAper)[maxNumElems] )
    {
      real64 const trans[2] = { stack.transmissibility[connectionIndex][0], stack.transmissibility[connectionIndex][1] };

      real64 energyFlux = 0.0;
      real64 dEnergyFlux_dTrans = 0.0;
      real64 dEnergyFlux_dP[2]{};
      real64 dEnergyFlux_dT[2]{};
      real64 dMassFlux_dT[2]{};

      // ***** Advective (enthalpy-carried) part *****
      singlePhaseFluxKernelsHelper::computeEnthalpyFlux( seri, sesri, sei,
                                                         trans,
                                                         m_enthalpy,
                                                         m_dEnthalpy,
                                                         m_gravCoef,
                                                         m_dDens,
                                                         m_dMob,
                                                         alpha,
                                                         mobility,
                                                         potGrad,
                                                         massFlux,
                                                         dMassFlux_dTrans,
                                                         dMassFlux_dP,
                                                         dMassFlux_dT,
                                                         energyFlux,
                                                         dEnergyFlux_dTrans,
                                                         dEnergyFlux_dP,
                                                         dEnergyFlux_dT );

      // Advective contribution to dEnergyFlux/dAperture: energyFlux = massFlux * enthalpyTimesMobWeight,
      // and enthalpy does not depend on aperture (for now). dMassFlux_dAper is this connection's only
      // (already multiplied by dt), with the same sign convention as the mass flux, for every element of the stencil.
      {
        real64 enthalpyTimesMobWeight = 0.0;
        if( alpha <= 0.0 || alpha >= 1.0 )
        {
          localIndex const k_up = 1 - localIndex( fmax( fmin( alpha, 1.0 ), 0.0 ) );
          enthalpyTimesMobWeight = m_enthalpy[seri[k_up]][sesri[k_up]][sei[k_up]][0];
        }
        else
        {
          real64 const mobWeights[2] = { alpha, 1.0 - alpha };
          for( integer ke = 0; ke < 2; ++ke )
          {
            enthalpyTimesMobWeight += mobWeights[ke] * m_enthalpy[seri[ke]][sesri[ke]][sei[ke]][0];
          }
        }

        for( localIndex j = 0; j < stack.numFluxElems; ++j )
        {
          stack.dEnergyFlux_dAperture[k[0]][j] += dMassFlux_dAper[j] * enthalpyTimesMobWeight;
          stack.dEnergyFlux_dAperture[k[1]][j] -= dMassFlux_dAper[j] * enthalpyTimesMobWeight;
        }
      }

      // add dMassFlux_dT to localFluxJacobian
      for( integer ke = 0; ke < 2; ++ke )
      {
        localIndex const localDofIndexTemp = k[ke] * numDof + numDof - 1;
        stack.localFluxJacobian[k[0]*numEqn][localDofIndexTemp] += m_dt * dMassFlux_dT[ke];
        stack.localFluxJacobian[k[1]*numEqn][localDofIndexTemp] -= m_dt * dMassFlux_dT[ke];
      }

      // ***** Conductive part (explicit conductivity) *****
      real64 const thermalTrans[2] = { stack.thermalTransmissibility[connectionIndex][0], stack.thermalTransmissibility[connectionIndex][1] };
      singlePhaseFluxKernelsHelper::computeConductiveFlux( seri, sesri, sei, m_temp, thermalTrans, energyFlux, dEnergyFlux_dT );

      // Conductive contribution to dEnergyFlux/dAperture (same convention as dFlux_dAperture, i.e. wrt the
      // normal jump): conductiveFlux = Trans_th * ( T_0 - T_1 ), and Trans_th depends linearly on the hydraulic aperture.
      {
        real64 const dConductiveFlux_dTrans = m_temp[seri[0]][sesri[0]][sei[0]] - m_temp[seri[1]][sesri[1]][sei[1]];

        for( localIndex j = 0; j < stack.numFluxElems; ++j )
        {
          real64 const dCondFlux_dAper = m_dt * dConductiveFlux_dTrans
                                         * stack.dThermalTrans_dLogHalfTrans[connectionIndex][j] * stack.dLogAperture_dNormalJump[j];
          stack.dEnergyFlux_dAperture[k[0]][j] += dCondFlux_dAper;
          stack.dEnergyFlux_dAperture[k[1]][j] -= dCondFlux_dAper;
        }
      }

      // ***** Assemble the energy flux of this connection and its derivatives *****
      stack.localFlux[k[0]*numEqn + numEqn - 1] += m_dt * energyFlux;
      stack.localFlux[k[1]*numEqn + numEqn - 1] -= m_dt * energyFlux;

      for( integer ke = 0; ke < 2; ++ke )
      {
        integer const localDofIndexPres = k[ke] * numDof;
        stack.localFluxJacobian[k[0]*numEqn + numEqn - 1][localDofIndexPres] += m_dt * dEnergyFlux_dP[ke];
        stack.localFluxJacobian[k[1]*numEqn + numEqn - 1][localDofIndexPres] -= m_dt * dEnergyFlux_dP[ke];
        integer const localDofIndexTemp = localDofIndexPres + 1;
        stack.localFluxJacobian[k[0]*numEqn + numEqn - 1][localDofIndexTemp] += m_dt * dEnergyFlux_dT[ke];
        stack.localFluxJacobian[k[1]*numEqn + numEqn - 1][localDofIndexTemp] -= m_dt * dEnergyFlux_dT[ke];
      }
    } );
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack ) const
  {
    // Call Base::complete to assemble the mass balance equations
    // In the lambda, add contribution to residual and jacobian into the energy balance equation
    Base::complete( iconn, stack, [&] ( integer const i,
                                        localIndex const localRow )
    {
      // The no. of fluxes is equal to the no. of equations in m_localRhs and m_localMatrix
      RAJA::atomicAdd( parallelDeviceAtomic{}, &SinglePhaseFVMAbstractBase::m_localRhs[localRow + numEqn-1], stack.localFlux[i * numEqn + numEqn-1] );

      SinglePhaseFVMAbstractBase::m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow + numEqn-1,
                                                                                                      stack.dofColIndices.data(),
                                                                                                      stack.localFluxJacobian[i * numEqn + numEqn-1].dataIfContiguous(),
                                                                                                      stack.stencilSize * numDof );

      // dEnergyFlux/dAperture (advective + conductive) -> energy-balance block of dR_dAper, appended after
      // all mass-balance rows
      if( m_dR_dAperEnergyOffset >= 0 )
      {
        localIndex const energyRow = m_dR_dAperEnergyOffset + LvArray::integerConversion< localIndex >( m_sei( iconn, i ) );
        m_dR_dAper.template addToRowBinarySearch< parallelDeviceAtomic >( energyRow,
                                                                          stack.localColIndices.data(),
                                                                          stack.dEnergyFlux_dAperture[i].dataIfContiguous(),
                                                                          stack.stencilSize );
      }

    } );
  }


private:

  /// Views on temperature
  ElementViewConst< arrayView1d< real64 const > > const m_temp;

  /// Views on enthalpies
  ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const m_enthalpy;

  /// Views on enthalpies
  ElementViewConst< arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > > const m_dEnthalpy;

  /// View on thermal conductivity
  ElementViewConst< arrayView3d< real64 const > > m_thermalConductivity;

  /// Row offset, in m_dR_dAper, of the energy-balance block.
  //  Negative (sentinel) means the caller did not build one - skip the energy-row write in complete().
  localIndex const m_dR_dAperEnergyOffset;

};



/**
 * @class ThermalConformingFracturesConnectorBasedAssemblyKernelFactory
 */
template< template< integer, integer > class BASE_KERNEL >
struct ThermalConformingFracturesConnectorBasedAssemblyKernelFactory
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
   */
  template< typename POLICY >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   string const & flowDofKey,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   SurfaceElementStencilWrapper const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs,
                   CRSMatrixView< real64, localIndex const > const & dR_dAper,
                   localIndex const dR_dAperEnergyOffset = -1 )
  {
    integer constexpr NUM_DOF = 2;   // pressure + temperature
    integer constexpr NUM_EQN = 2;   // mass balance + energy balance


    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > flowDofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( flowDofKey );
    flowDofNumberAccessor.setName( solverName + "/accessors/" + flowDofKey );

    using kernelType = ThermalConformingFracturesConnectorBasedAssemblyKernel< NUM_EQN, NUM_DOF, BASE_KERNEL >;
    typename kernelType::SinglePhaseFlowAccessors flowAccessors( elemManager, solverName );
    typename kernelType::ThermalSinglePhaseFlowAccessors thermalFlowAccessors( elemManager, solverName );

    typename kernelType::SinglePhaseFluidAccessors fluidAccessors( elemManager, solverName );
    typename kernelType::ThermalSinglePhaseFluidAccessors thermalFluidAccessors( elemManager, solverName );

    typename kernelType::PermeabilityAccessors permAccessors( elemManager, solverName );
    typename kernelType::FracturePermeabilityAccessors edfmPermAccessors( elemManager, solverName );
    typename kernelType::ThermalConductivityAccessors thermalConductivityAccessors( elemManager, solverName );
    typename kernelType::FractureApertureAccessors fractureApertureAccessors( elemManager, solverName );

    kernelType kernel( rankOffset, stencilWrapper,
                       flowDofNumberAccessor,
                       flowAccessors, thermalFlowAccessors, fluidAccessors, thermalFluidAccessors,
                       permAccessors, edfmPermAccessors, thermalConductivityAccessors, fractureApertureAccessors,
                       dt, localMatrix, localRhs, dR_dAper, dR_dAperEnergyOffset );

    kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
  }
};

} // namespace thermalConformingFracturesKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICSCONFORMINGFRACTURESKERNELBASE_HPP
