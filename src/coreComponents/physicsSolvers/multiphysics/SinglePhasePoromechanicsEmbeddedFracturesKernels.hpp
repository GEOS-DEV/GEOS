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
 * @file singlePhasePoromechanicsEmbeddedFractureKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSEMBEDDEDFRACTURESKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSEMBEDDEDFRACTURESKERNELS_HPP

#include "physicsSolvers/fluidFlow/SinglePhaseFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/FluxKernelsHelper.hpp"

namespace geosx
{

namespace singlePhasePoromechanicsEmbeddedFracturesKernels
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
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_mob;
  using AbstractBase::m_dens;

  using Base = singlePhaseFVMKernels::FaceBasedAssemblyKernel< NUM_DOF, SurfaceElementStencilWrapper >;
  using Base::numDof;
  using Base::maxNumElems;
  using Base::maxNumConns;
  using Base::maxStencilSize;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;


  /// Compute time value for the number of equations
  static constexpr integer numEqn = numDof - 3;

  ConnectorBasedAssemblyKernel( globalIndex const rankOffset,
                                STENCILWRAPPER const & stencilWrapper,
                                DofNumberAccessor const & flowDofNumberAccessor,
                                DofNumberAccessor const & dispDofNumberAccessor,
                                SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                                ThermalSinglePhaseFlowAccessors const & thermalSinglePhaseFlowAccessors,
                                SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                                ThermalSinglePhaseFluidAccessors const & thermalSinglePhaseFluidAccessors,
                                PermeabilityAccessors const & permeabilityAccessors,
                                ThermalConductivityAccessors const & thermalConductivityAccessors,
                                real64 const & dt,
                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                arrayView1d< real64 > const & localRhs )
    : Base( rankOffset,
            stencilWrapper,
            dofNumberAccessor,
            singlePhaseFlowAccessors,
            singlePhaseFluidAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs ),
    m_dispJumpDofNumber( dispDofNumberAccessor.toNestedViewConst() ),
    m_dPerm_dDispJump( thermalSinglePhaseFlowAccessors.get( fields::flow::temperature {} ) )
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
    real64 dTrans_dDispJump[maxNumConns][2][3]{};
  };

  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] iconn the connection index
   * @param[in] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void setup( localIndex const iconn,
              StackVariables & stack ) const
  {
    // set degrees of freedom indices for this face
    for( integer i = 0; i < stack.stencilSize; ++i )
    {
      localIndex localDofIndex = numDof * i;
      dofColIndices[ localDofIndex ]     = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];
      dofColIndices[ localDofIndex + 1 ] = m_dispJumpDofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];
      dofColIndices[ localDofIndex + 2 ] = m_dispJumpDofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )] + 1;
      dofColIndices[ localDofIndex + 3 ] = m_dispJumpDofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )] + 2;
    }
  }

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


    real64 fluxVal = 0.0;
    real64 dFlux_dTrans = 0.0;

    /// EDFM connections are always only between 2 elements. There are no star connections.
    real64 trans[2] = {stack.transmissibility[0][0], stack.transmissibility[0][1]};
    real64 dTrans[2] = { stack.dTrans_dPres[0][0], stack.dTrans_dPres[0][1] };
    real64 dFlux_dP[2] = {0.0, 0.0};
    localIndex const regionIndex[2]    = {m_seri[iconn][0], m_seri[iconn][1]};
    localIndex const subRegionIndex[2] = {m_sesri[iconn][0], m_sesri[iconn][1]};
    localIndex const elementIndex[2]   = {m_sei[iconn][0], m_sei[iconn][1]};


    computeSinglePhaseFlux( regionIndex, subRegionIndex, elementIndex,
                            trans,
                            dTrans,
                            pres,
                            gravCoef,
                            dens,
                            dDens_dPres,
                            mob,
                            dMob_dPres,
                            fluxVal,
                            dFlux_dP,
                            dFlux_dTrans );



    // populate local flux vector and derivatives
    stack.localFlux[0] =  m_dt * fluxVal;
    stack.localFlux[1] = -m_dt * fluxVal;

    real64 dFlux_dDispJump[2][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    for( localIndex i=0; i < 3; i++ )
    {
      dFlux_dDispJump[0][i] = m_dt * dFlux_dTrans * dTrans_dDispJump[0][0][i];
      dFlux_dDispJump[1][i] = -m_dt * dFlux_dTrans * dTrans_dDispJump[0][1][i];
    }
    for( localIndex ke = 0; ke < 2; ++ke )
    {
      localIndex const dofIndex = 4*ke;

      stack.localFluxJacobian[0][dofIndex]   =  m_dt * dFlux_dP[ke];
      stack.localFluxJacobian[0][dofIndex+1] =  m_dt * dFlux_dDispJump[ke][0];
      stack.localFluxJacobian[0][dofIndex+2] =  m_dt * dFlux_dDispJump[ke][1];
      stack.localFluxJacobian[0][dofIndex+3] =  m_dt * dFlux_dDispJump[ke][2];

      stack.localFluxJacobian[1][dofIndex]   = -m_dt * dFlux_dP[ke];
      stack.localFluxJacobian[1][dofIndex+1] = -m_dt * dFlux_dDispJump[ke][0];
      stack.localFluxJacobian[1][dofIndex+2] = -m_dt * dFlux_dDispJump[ke][1];
      stack.localFluxJacobian[1][dofIndex+3] = -m_dt * dFlux_dDispJump[ke][2];
    }

    kernelOp( );
  }


private:

  ElementViewConst< arrayView1d< globalIndex const > > const m_dispJumpDofNumber;

  ElementViewConst< arrayView4d< real64 const > > const m_dPerm_dDispJump;

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
  template< typename POLICY, typename STENCILWRAPPER >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   string const & pressureDofKey,
                   string const & dispJumpDofKey,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    integer constexpr NUM_DOF = 4; // pressure + jumps

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > pressureDofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( pressureDofKey );
    dofNumberAccessor.setName( solverName + "/accessors/" + pressureDofKey );

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dispJumpDofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( dispJumpDofKey );
    dofNumberAccessor.setName( solverName + "/accessors/" + dispJumpDofKey );

    using kernelType = ConnectorBasedAssemblyKernel< NUM_DOF, STENCILWRAPPER >;
    typename kernelType::SinglePhaseFlowAccessors flowAccessors( elemManager, solverName );
    typename kernelType::SinglePhaseFluidAccessors fluidAccessors( elemManager, solverName );
    typename kernelType::PermeabilityAccessors permAccessors( elemManager, solverName );

    kernelType kernel( rankOffset, stencilWrapper, pressureDofNumberAccessor,
                       flowAccessors, fluidAccessors, permAccessors,
                       dt, localMatrix, localRhs );

    kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
  }
};



} // namespace SinglePhasePoromechanicsFluxKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSEMBEDDEDFRACTURSKERNELS_HPP
