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
 * @file SinglePhasePoromechanicsConformingFractures.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP

#include "physicsSolvers/fluidFlow/SinglePhaseFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/FluxKernelsHelper.hpp"

namespace geosx
{

namespace singlePhasePoromechanicsConformingFracturesKernels
{

using namespace fluxKernelsHelper;

template< integer NUM_DOF >
class ConnectorBasedAssemblyKernel : public singlePhaseFVMKernels::FaceBasedAssemblyKernel< NUM_DOF, SurfaceElementStencilWrapper >
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

  using AbstractBase = singlePhaseFVMKernels::FaceBasedAssemblyKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using SinglePhaseFlowAccessors = AbstractBase::SinglePhaseFlowAccessors;
  using SinglePhaseFluidAccessors = AbstractBase::SinglePhaseFluidAccessors;
  using PermeabilityAccessors = StencilMaterialAccessors< PermeabilityBase,
                                                          fields::permeability::permeability,
                                                          fields::permeability::dPerm_dPressure,
                                                          fields::permeability::dPerm_dDispJump >;

  using AbstractBase::m_dt;
  using AbstractBase::m_rankOffset;
  using AbstractBase::m_dofNumber;
  using AbstractBase::m_permeability;
  using AbstractBase::m_dPerm_dPres;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_pres;
  using AbstractBase::m_mob;
  using AbstractBase::m_dMob_dPres;
  using AbstractBase::m_dens;
  using AbstractBase::m_dDens_dPres;

  using Base = singlePhaseFVMKernels::FaceBasedAssemblyKernel< NUM_DOF, SurfaceElementStencilWrapper >;
  using Base::numDof;
  using Base::numEqn;
  using Base::maxNumElems;
  using Base::maxNumConns;
  using Base::maxStencilSize;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;

  ConnectorBasedAssemblyKernel( globalIndex const rankOffset,
                                STENCILWRAPPER const & stencilWrapper,
                                DofNumberAccessor const & flowDofNumberAccessor,
                                SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                                ThermalSinglePhaseFlowAccessors const & thermalSinglePhaseFlowAccessors,
                                SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                                ThermalSinglePhaseFluidAccessors const & thermalSinglePhaseFluidAccessors,
                                PermeabilityAccessors const & permeabilityAccessors,
                                ThermalConductivityAccessors const & thermalConductivityAccessors,
                                real64 const & dt,
                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                arrayView1d< real64 > const & localRhs,
                                CRSMatrixView< real64, localIndex const > const & dR_dAper )
    : Base( rankOffset,
            stencilWrapper,
            dofNumberAccessor,
            singlePhaseFlowAccessors,
            singlePhaseFluidAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs ),
    m_dR_dAper( dR_dAper )
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
    GEOSX_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : Base::StackVariables( size, numElems )
    {}

    /// Derivatives of transmissibility with respect to pressure
  };

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the flux
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] NoOpFunc the function used to customize the computation of the flux
   */
  template< typename FUNC = singlePhaseBaseKernels::NoOpFunc >
  GEOSX_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack,
                    FUNC && kernelOp = singlePhaseBaseKernels::NoOpFunc{} ) const

  {

    stencilWrapper.computeWeights( iconn, 
                                   m_permeability,
                                   m_dPerm_dPres,
                                   m_dPerm_dDispJump,
                                   stack.transmissibility,
                                   stack.dTrans_dPres,
                                   stack.dTrans_dDispJump );


    localIndex k[2];
    localIndex connectionIndex = 0;
    for( k[0]=0; k[0]<stack.numFluxElems; ++k[0] )
    {
      for( k[1]=k[0]+1; k[1]<stack.numFluxElems; ++k[1] )
      {
        real64 fluxVal = 0.0;
        real64 dFlux_dTrans = 0.0;
        real64 const trans[2] = { stack.transmissibility[connectionIndex][0], stack.transmissibility[connectionIndex][1] };
        real64 const dTrans[2] = { stack.dTrans_dPres[connectionIndex][0], stack.dTrans_dPres[connectionIndex][1] };
        real64 const dFlux_dP[2] = { 0.0, 0.0 };
        localIndex const regionIndex[2]    = {m_seri[k[0]], m_seri[k[1]]};
        localIndex const subRegionIndex[2] = {m_sesri[k[0]], m_sesri[k[1]]};
        localIndex const elementIndex[2]   = {m_sei[k[0]], m_sei[k[1]]};

        computeSinglePhaseFlux( regionIndex, subRegionIndex, elementIndex,
                                trans,
                                dTrans,
                                m_pres,
                                m_gravCoef,
                                m_dens,
                                m_dDens_dPres,
                                m_mob,
                                m_dMob_dPres,
                                fluxVal,
                                dFlux_dP,
                                dFlux_dTrans );

        // populate local flux vector and derivatives
        stack.localFlux[k[0]] += m_dt * fluxVal;
        stack.localFlux[k[1]] -= m_dt * fluxVal;

        real64 dFlux_dAper[2] = {0.0, 0.0};
        dFlux_dAper[0] =  m_dt * dFlux_dTrans * dTrans_dDispJump[connectionIndex][0][0];
        dFlux_dAper[1] = -m_dt * dFlux_dTrans * dTrans_dDispJump[connectionIndex][1][0];

        stack.localFluxJacobian[k[0]][k[0]] += dFlux_dP[0] * m_dt;
        stack.localFluxJacobian[k[0]][k[1]] += dFlux_dP[1] * m_dt;
        stack.localFluxJacobian[k[1]][k[0]] -= dFlux_dP[0] * m_dt;
        stack.localFluxJacobian[k[1]][k[1]] -= dFlux_dP[1] * m_dt;

        stack.dFlux_dAperture[k[0]][k[0]] += dFlux_dAper[0];
        stack.dFlux_dAperture[k[0]][k[1]] += dFlux_dAper[1];
        stack.dFlux_dAperture[k[1]][k[0]] -= dFlux_dAper[0];
        stack.dFlux_dAperture[k[1]][k[1]] -= dFlux_dAper[1];

        kernelOp( );
        connectionIndex++;
      }
    }
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  template< typename FUNC = singlePhaseBaseKernels::NoOpFunc >
  GEOSX_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack,
                 FUNC && kernelOp = singlePhaseBaseKernels::NoOpFunc{} ) const
  {
    // Call Base::complete to assemble the mass balance equations
    // In the lambda, fill the dR_dAper matrix
    Base::complete( iconn, stack, [&] ( integer const i,
                                        localIndex const localRow )
    {
      m_dR_dAper.addToRowBinarySearch< parallelDeviceAtomic >( m_sei( iconn, i ),
                                                               stack.localColIndices.data(),
                                                               stack.dFlux_dAper[i].dataIfContiguous(),
                                                               stack.stencilSize );
      // call the lambda to assemble additional terms, such as thermal terms
      kernelOp( i, localRow );
    } );
  }

private:

  CRSMatrixView< real64, localIndex const > m_dR_dAper;

};



/**
 * @class FaceBasedAssemblyKernelFactory
 */
class ConnectorBasedAssemblyKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @tparam STENCILWRAPPER the type of the stencil wrapper
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
                   string const & pressureDofKey,
                   string const & dispJumpDofKey,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   SurfaceElementStencilWrapper const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs,
                   CRSMatrixView< real64, localIndex const > const & dR_dAper )
  {
    integer constexpr NUM_DOF = 1; // pressure

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > flowDofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( pressureDofKey );
    dofNumberAccessor.setName( solverName + "/accessors/" + pressureDofKey );

    using kernelType = ConnectorBasedAssemblyKernel< NUM_DOF, SurfaceElementStencilWrapper >;
    typename kernelType::SinglePhaseFlowAccessors flowAccessors( elemManager, solverName );
    typename kernelType::SinglePhaseFluidAccessors fluidAccessors( elemManager, solverName );
    typename kernelType::PermeabilityAccessors permAccessors( elemManager, solverName );

    kernelType kernel( rankOffset, stencilWrapper, flowDofNumberAccessor,
                       flowAccessors, fluidAccessors, permAccessors,
                       dt, localMatrix, localRhs, dR_dAper );

    kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
  }
};

} // namespace SinglePhasePoromechanicsConformingFracturesKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP
