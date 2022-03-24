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
 * @file StabilizedCompositionalMultiphaseFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_STABILIZEDCOMPOSITIONALMULTIPHASEFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_STABILIZEDCOMPOSITIONALMULTIPHASEFVMKERNELS_HPP

#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVMKernels.hpp"

namespace geosx
{

namespace stabilizedCompositionalMultiphaseFVMKernels
{

using namespace constitutive;


/******************************** FaceBasedAssemblyKernel ********************************/

/**
 * @class FaceBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NUM_COMP, integer NUM_DOF, typename STENCILWRAPPER >
class FaceBasedAssemblyKernel : public compositionalMultiphaseFVMKernels::FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >
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

  using AbstractBase = compositionalMultiphaseFVMKernels::FaceBasedAssemblyKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using CompFlowAccessors = AbstractBase::CompFlowAccessors;
  using MultiFluidAccessors = AbstractBase::MultiFluidAccessors;
  using CapPressureAccessors = AbstractBase::CapPressureAccessors;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;

  using StabCompFlowAccessors = StencilAccessors< extrinsicMeshData::flow::phaseVolumeFractionOld,
                                                  extrinsicMeshData::flow::phaseDensityOld,
                                                  extrinsicMeshData::flow::phaseComponentFractionOld >;

  using AbstractBase::m_dt;
  using AbstractBase::m_numPhases;
  using AbstractBase::m_hasCapPressure;
  using AbstractBase::m_rankOffset;
  using AbstractBase::m_ghostRank;
  using AbstractBase::m_dofNumber;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_phaseMob;
  using AbstractBase::m_dPhaseMassDens;
  using AbstractBase::m_phaseCompFrac;
  using AbstractBase::m_dPhaseCompFrac;
  using AbstractBase::m_dCompFrac_dCompDens;
  using AbstractBase::m_dPhaseCapPressure_dPhaseVolFrac;

  using Base = compositionalMultiphaseFVMKernels::FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
  using Base::numComp;
  using Base::numDof;
  using Base::numEqn;
  using Base::maxNumElems;
  using Base::maxNumConns;
  using Base::maxStencilSize;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;


  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] hasCapPressure flag specifying whether capillary pressure is used or not
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor accessor for the dofs numbers
   * @param[in] compFlowAccessor accessor for wrappers registered by the solver
   * @param[in] multiFluidAccessor accessor for wrappers registered by the multifluid model
   * @param[in] capPressureAccessors accessor for wrappers registered by the cap pressure model
   * @param[in] permeabilityAccessors accessor for wrappers registered by the permeability model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FaceBasedAssemblyKernel( integer const numPhases,
                           globalIndex const rankOffset,
                           integer const hasCapPressure,
                           STENCILWRAPPER const & stencilWrapper,
                           DofNumberAccessor const & dofNumberAccessor,
                           CompFlowAccessors const & compFlowAccessors,
                           StabCompFlowAccessors const & stabCompFlowAccessors,
                           MultiFluidAccessors const & multiFluidAccessors,
                           CapPressureAccessors const & capPressureAccessors,
                           PermeabilityAccessors const & permeabilityAccessors,
                           real64 const & dt,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs )
    : Base( numPhases,
            rankOffset,
            hasCapPressure,
            stencilWrapper,
            dofNumberAccessor,
            compFlowAccessors,
            multiFluidAccessors,
            capPressureAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs ), 
    m_phaseMassDensOld(stabCompFlowAccessors.get(extrinsicMeshData::flow::phaseDensityOld {}) ),
    m_phaseCompFracOld(stabCompFlowAccessors.get(extrinsicMeshData::flow::phaseComponentFractionOld {})),
    m_phaseVolFracOld(stabCompFlowAccessors.get(extrinsicMeshData::flow::phaseVolumeFractionOld {})),
    m_stabWeights(stencilWrapper.getStabWeights())

  {}

  struct StackVariables : public Base::StackVariables
  {
public:

    GEOSX_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : Base::StackVariables( size, numElems ),
      stabFlux( 0.0 ),
      dStabFlux_dP( size )
    {}

    using Base::StackVariables::stencilSize;
    using Base::StackVariables::numFluxElems;
    using Base::StackVariables::transmissibility;
    using Base::StackVariables::dTrans_dPres;
    using Base::StackVariables::dofColIndices;
    using Base::StackVariables::localFlux;
    using Base::StackVariables::localFluxJacobian;

    real64 stabFlux;
    stackArray1d< real64, maxStencilSize > dStabFlux_dP;

  };

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack ) const
  {
    // using Deriv = multifluid::DerivativeOffset;


    // ***********************************************
    // First, we call the base computeFlux to compute:
    //  1) compFlux and its derivatives,
    //
    // Computing stabilization flux requires quantities already computed in the base computeFlux,
    // such as potGrad, phaseFlux, and the indices of the upwind cell
    // We use the lambda below (called **inside** the phase loop of the base computeFlux) to access these variables
    Base::computeFlux( iconn, stack, [=] GEOSX_HOST_DEVICE ( integer const ip,
                                                             localIndex const k_up,
                                                             localIndex const er_up,
                                                             localIndex const esr_up,
                                                             localIndex const ei_up,
                                                             real64 const & potGrad,
                                                             real64 const & phaseFlux,
                                                             real64 const (&dPhaseFlux_dP)[maxStencilSize],
                                                             real64 const (&dPhaseFlux_dC)[maxStencilSize][numComp] )
    {
      GEOSX_UNUSED_VAR( ip, k_up, er_up, esr_up, ei_up, potGrad, phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC );

      // // We are in the loop over phases, ip provides the current phase index.
      
      // real64 dPresGrad{};

      // real64 tauStab = 0.0; // need to get correct value


      // // compute potential difference MPFA-style
      // for( integer i = 0; i < stack.stencilSize; ++i )
      // {
      //   localIndex const er  = m_seri( iconn, i );
      //   localIndex const esr = m_sesri( iconn, i );
      //   localIndex const ei  = m_sei( iconn, i );

      //   dPresGrad += tauStab*m_stabWeights[iconn][i]*m_dPres[er][esr][ei]; // need to figure out indices into stabWeigts and tau

      // }

      // // modify flux with stabilization
      // // multiply dPresGrad with upwind quantities


    } ); // end call to Base::computeFlux

    
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack ) const
  {
    // Step 1: assemble the component mass balance equations (i = 0 to i = numDof-2)
    Base::complete( iconn, stack );

  }

protected:

  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const  m_phaseMassDensOld;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_COMP > > const m_phaseCompFracOld;
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const m_phaseVolFracOld;

  typename STENCILWRAPPER::WeightContainerViewConstType m_stabWeights;

};

/**
 * @class FaceBasedAssemblyKernelFactory
 */
class FaceBasedAssemblyKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @tparam STENCILWRAPPER the type of the stencil wrapper
   * @param[in] numComps the number of fluid components
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] hasCapPressure flag specifying whether capillary pressure is used or not
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY, typename STENCILWRAPPER >
  static void
  createAndLaunch( integer const numComps,
                   integer const numPhases,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   integer const hasCapPressure,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    compositionalMultiphaseBaseKernels::
      internal::kernelLaunchSelectorCompSwitch( numComps, [&] ( auto NC )
    {
      integer constexpr NUM_COMP = NC();
      integer constexpr NUM_DOF = NC()+2;

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

      using KERNEL_TYPE = FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
      typename KERNEL_TYPE::CompFlowAccessors compFlowAccessors( elemManager, solverName );
      typename KERNEL_TYPE::MultiFluidAccessors multiFluidAccessors( elemManager, solverName );
      typename KERNEL_TYPE::StabCompFlowAccessors stabCompFlowAccessors( elemManager, solverName );
      typename KERNEL_TYPE::CapPressureAccessors capPressureAccessors( elemManager, solverName );
      typename KERNEL_TYPE::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );

      KERNEL_TYPE kernel( numPhases, rankOffset, hasCapPressure, stencilWrapper, dofNumberAccessor,
                          compFlowAccessors, stabCompFlowAccessors, multiFluidAccessors,
                          capPressureAccessors, permeabilityAccessors,
                          dt, localMatrix, localRhs );
      KERNEL_TYPE::template launch< POLICY >( stencilWrapper.size(), kernel );
    } );
  }
};


} // namespace stabilizedCompositionalMultiphaseFVMKernels

} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_STABILIZEDCOMPOSITIONALMULTIPHASEFVMKERNELS_HPP