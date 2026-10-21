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
 * @file SinglePhasePoromechanicsConformingFractures.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP

#include "physicsSolvers/multiphysics/poromechanicsKernels/ConformingFracturesFluxKernelBase.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/FluxComputeKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/FluxKernelsHelper.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geos
{

namespace singlePhasePoromechanicsConformingFracturesKernels
{

template< integer NUM_EQN, integer NUM_DOF >
class ConnectorBasedAssemblyKernel
  : public poromechanicsConformingFracturesKernels::
    ConformingFracturesFluxKernelBase< singlePhaseFVMKernels::FluxComputeKernel< NUM_EQN, NUM_DOF, SurfaceElementStencilWrapper > >
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

  using AbstractBase = singlePhaseFVMKernels::FluxComputeKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using SinglePhaseFlowAccessors = AbstractBase::SinglePhaseFlowAccessors;
  using SinglePhaseFluidAccessors = AbstractBase::SinglePhaseFluidAccessors;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;

  using AbstractBase::m_dt;
  using AbstractBase::m_rankOffset;
  using AbstractBase::m_dofNumber;
  using AbstractBase::m_permeability;
  using AbstractBase::m_dPerm_dPres;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_pres;
  using AbstractBase::m_mob;
  using AbstractBase::m_dMob;
  using AbstractBase::m_dMob_dPres;
  using AbstractBase::m_dens;
  using AbstractBase::m_dDens;


  using Base = poromechanicsConformingFracturesKernels::
                 ConformingFracturesFluxKernelBase< singlePhaseFVMKernels::FluxComputeKernel< NUM_EQN, NUM_DOF, SurfaceElementStencilWrapper > >;
  using FracturePermeabilityAccessors = typename Base::FracturePermeabilityAccessors;
  using FractureApertureAccessors = typename Base::FractureApertureAccessors;
  using Base::numDof;
  using Base::numEqn;
  using Base::maxNumElems;
  using Base::maxNumConns;
  using Base::maxStencilSize;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;
  using Base::m_dR_dAper;

  ConnectorBasedAssemblyKernel( globalIndex const rankOffset,
                                SurfaceElementStencilWrapper const & stencilWrapper,
                                DofNumberAccessor const & flowDofNumberAccessor,
                                SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                                SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                                PermeabilityAccessors const & permeabilityAccessors,
                                FracturePermeabilityAccessors const & fracturePermeabilityAccessors,
                                FractureApertureAccessors const & fractureApertureAccessors,
                                real64 const & dt,
                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                arrayView1d< real64 > const & localRhs,
                                CRSMatrixView< real64, localIndex const > const & dR_dAper )
    : Base( fracturePermeabilityAccessors,
            fractureApertureAccessors,
            dR_dAper,
            rankOffset,
            stencilWrapper,
            flowDofNumberAccessor,
            singlePhaseFlowAccessors,
            singlePhaseFluidAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs )
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
      dFlux_dAperture( numElems, size )
    {}

    stackArray2d< real64, maxNumElems * maxStencilSize > dFlux_dAperture;
  };

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the flux
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] NoOpFunc the function used to customize the computation of the flux
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack,
                    FUNC && kernelOp = NoOpFunc{} ) const

  {
    Base::computeTransmissibilityAndApertureDerivatives( iconn, stack );

    localIndex k[2];
    localIndex connectionIndex = 0;
    for( k[0]=0; k[0]<stack.numFluxElems; ++k[0] )
    {
      for( k[1]=k[0]+1; k[1]<stack.numFluxElems; ++k[1] )
      {
        real64 fluxVal = 0.0;
        real64 dFlux_dTrans = 0.0;
        real64 alpha = 0.0;
        real64 mobility = 0.0;
        real64 potGrad = 0.0;
        real64 const trans[2] = { stack.transmissibility[connectionIndex][0], stack.transmissibility[connectionIndex][1] };
        real64 const dTrans[2] = { stack.dTrans_dPres[connectionIndex][0], stack.dTrans_dPres[connectionIndex][1] };
        real64 dFlux_dP[2] = { 0.0, 0.0 };
        localIndex const regionIndex[2]    = {m_seri[iconn][k[0]], m_seri[iconn][k[1]]};
        localIndex const subRegionIndex[2] = {m_sesri[iconn][k[0]], m_sesri[iconn][k[1]]};
        localIndex const elementIndex[2]   = {m_sei[iconn][k[0]], m_sei[iconn][k[1]]};

        singlePhaseFluxKernelsHelper::computeSinglePhaseFlux( regionIndex, subRegionIndex, elementIndex,
                                                              trans,
                                                              dTrans,
                                                              m_pres,
                                                              m_gravCoef,
                                                              m_dens,
                                                              m_dDens,
                                                              m_mob,
                                                              m_dMob,
                                                              alpha,
                                                              mobility,
                                                              potGrad,
                                                              fluxVal,
                                                              dFlux_dP,
                                                              dFlux_dTrans );

        // populate local flux vector and derivatives
        stack.localFlux[k[0]* numDof] += m_dt * fluxVal;
        stack.localFlux[k[1]* numDof] -= m_dt * fluxVal;

        stack.localFluxJacobian[k[0]*numEqn][k[0]* numDof] += dFlux_dP[0] * m_dt;
        stack.localFluxJacobian[k[0]*numEqn][k[1]* numDof] += dFlux_dP[1] * m_dt;
        stack.localFluxJacobian[k[1]*numEqn][k[0]* numDof] -= dFlux_dP[0] * m_dt;
        stack.localFluxJacobian[k[1]*numEqn][k[1]* numDof] -= dFlux_dP[1] * m_dt;

        // dFlux/dNormalJump_j for every element j of the stencil (at a junction, also those this connection does not link)
        real64 dFlux_dAper[maxNumElems]{};
        Base::accumulateFluxApertureDerivative( stack, k, connectionIndex, dFlux_dTrans, stack.dFlux_dAperture, dFlux_dAper );

        // per-connection quantities (incl. this connection's aperture derivative of the flux, before it is
        // accumulated over the connections of the stencil), e.g. for the thermal terms at fracture junctions
        kernelOp( k, regionIndex, subRegionIndex, elementIndex, iconn, connectionIndex,
                  alpha, mobility, potGrad, fluxVal, dFlux_dTrans, dFlux_dP, dFlux_dAper );
        connectionIndex++;
      }
    }
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack,
                 FUNC && kernelOp = NoOpFunc{} ) const
  {
    // Call Base::complete to assemble the mass balance equations
    // In the lambda, fill the dR_dAper matrix
    Base::complete( iconn, stack, [&] ( integer const i,
                                        localIndex const localRow )
    {

      localIndex const row = LvArray::integerConversion< localIndex >( m_sei( iconn, i ) );

      m_dR_dAper.template addToRowBinarySearch< parallelDeviceAtomic >( row,
                                                                        stack.localColIndices.data(),
                                                                        stack.dFlux_dAperture[i].dataIfContiguous(),
                                                                        stack.stencilSize );
      // call the lambda to assemble additional terms, such as thermal terms
      kernelOp( i, localRow );
    } );
  }

};



/**
 * @class ConnectorBasedAssemblyKernelFactory
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
                   string const & dofKey,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   SurfaceElementStencilWrapper const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs,
                   CRSMatrixView< real64, localIndex const > const & dR_dAper )
  {
    integer constexpr NUM_DOF = 1; // pressure
    integer constexpr NUM_EQN = 1;

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > flowDofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    flowDofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

    using kernelType = ConnectorBasedAssemblyKernel< NUM_EQN, NUM_DOF >;
    typename kernelType::SinglePhaseFlowAccessors flowAccessors( elemManager, solverName );
    typename kernelType::SinglePhaseFluidAccessors fluidAccessors( elemManager, solverName );
    typename kernelType::PermeabilityAccessors permAccessors( elemManager, solverName );
    typename kernelType::FracturePermeabilityAccessors fracPermAccessors( elemManager, solverName );
    typename kernelType::FractureApertureAccessors fractureApertureAccessors( elemManager, solverName );

    kernelType kernel( rankOffset, stencilWrapper, flowDofNumberAccessor,
                       flowAccessors, fluidAccessors, permAccessors, fracPermAccessors, fractureApertureAccessors,
                       dt, localMatrix, localRhs, dR_dAper );

    kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
  }
};

} // namespace SinglePhasePoromechanicsConformingFracturesKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP
