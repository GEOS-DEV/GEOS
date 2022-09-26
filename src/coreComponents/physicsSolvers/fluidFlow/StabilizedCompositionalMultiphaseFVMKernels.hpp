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

#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseFVMKernels.hpp"

#include "constitutive/solid/SolidBase.hpp"
#include "constitutive/solid/SolidExtrinsicData.hpp"
#include "constitutive/solid/porosity/PorosityExtrinsicData.hpp"

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
class FaceBasedAssemblyKernel : public isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >
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

  using AbstractBase = isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using CompFlowAccessors = AbstractBase::CompFlowAccessors;
  using MultiFluidAccessors = AbstractBase::MultiFluidAccessors;
  using CapPressureAccessors = AbstractBase::CapPressureAccessors;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;

  using StabCompFlowAccessors =
    StencilAccessors< extrinsicMeshData::flow::macroElementIndex,
                      extrinsicMeshData::flow::pressure_n >;

  using StabMultiFluidAccessors =
    StencilMaterialAccessors< MultiFluidBase,
                              extrinsicMeshData::multifluid::phaseDensity_n,
                              extrinsicMeshData::multifluid::phaseCompFraction_n >;

  using SolidAccessors =
    StencilMaterialAccessors< SolidBase,
                              extrinsicMeshData::solid::bulkModulus,
                              extrinsicMeshData::solid::shearModulus >;

  using PorosityAccessors =
    StencilMaterialAccessors< PorosityBase,
                              extrinsicMeshData::porosity::biotCoefficient >;

  using RelPermAccessors =
    StencilMaterialAccessors< RelativePermeabilityBase, extrinsicMeshData::relperm::phaseRelPerm_n >;

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
  using AbstractBase::m_pres;

  using Base = isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
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
   * @param[in] stabCompFlowAccessor accessor for wrappers registered by the solver needed for stabilization
   * @param[in] multiFluidAccessor accessor for wrappers registered by the multifluid model
   * @param[in] stabMultiFluidAccessor accessor for wrappers registered by the multifluid model needed for stabilization
   * @param[in] capPressureAccessors accessor for wrappers registered by the cap pressure model
   * @param[in] permeabilityAccessors accessor for wrappers registered by the permeability model
   * @param[in] solidAccessors accessor for wrappers registered by the solid model needed for stabilization
   * @param[in] porosityAccessors accessor for wrappers registered by the porosity model needed for stabilization
   * @param[in] relPermAccessors accessor for wrappers registered by the relative permeability model needed for stabilization
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
                           StabMultiFluidAccessors const & stabMultiFluidAccessors,
                           CapPressureAccessors const & capPressureAccessors,
                           PermeabilityAccessors const & permeabilityAccessors,
                           SolidAccessors const & solidAccessors,
                           PorosityAccessors const & porosityAccessors,
                           RelPermAccessors const & relPermAccessors,
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
    m_pres_n( stabCompFlowAccessors.get( extrinsicMeshData::flow::pressure_n {} ) ),
    m_phaseDens_n( stabMultiFluidAccessors.get( extrinsicMeshData::multifluid::phaseDensity_n {} ) ),
    m_phaseCompFrac_n( stabMultiFluidAccessors.get( extrinsicMeshData::multifluid::phaseCompFraction_n {} ) ),
    m_phaseRelPerm_n( relPermAccessors.get( extrinsicMeshData::relperm::phaseRelPerm_n {} ) ),
    m_macroElementIndex( stabCompFlowAccessors.get( extrinsicMeshData::flow::macroElementIndex {} ) ),
    m_bulkModulus( solidAccessors.get( extrinsicMeshData::solid::bulkModulus {} ) ),
    m_shearModulus( solidAccessors.get( extrinsicMeshData::solid::shearModulus {} ) ),
    m_biotCoefficient( porosityAccessors.get( extrinsicMeshData::porosity::biotCoefficient {} ) )
  {}

  struct StackVariables : public Base::StackVariables
  {
public:

    GEOSX_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : Base::StackVariables( size, numElems ),
      stabFlux( numElems * numEqn ),
      dStabFlux_dP( numElems * numEqn, size * numDof )
    {}

    using Base::StackVariables::stencilSize;
    using Base::StackVariables::numFluxElems;
    using Base::StackVariables::transmissibility;
    using Base::StackVariables::dTrans_dPres;
    using Base::StackVariables::dofColIndices;
    using Base::StackVariables::localFlux;
    using Base::StackVariables::localFluxJacobian;

    /// Stabilization transmissibility
    real64 stabTransmissibility[maxNumConns][2]{};

    /// Stabilization flux and derivatives
    stackArray1d< real64, maxNumElems * numEqn > stabFlux;
    stackArray2d< real64, maxNumElems * numEqn * maxStencilSize * numDof > dStabFlux_dP;

  };

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian, including stabilization
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack ) const
  {

    m_stencilWrapper.computeStabilizationWeights( iconn,
                                                  stack.stabTransmissibility );

    // ***********************************************
    // First, we call the base computeFlux to compute:
    //  1) compFlux and its derivatives,
    //
    // We use the lambda below (called **inside** the phase loop of the base computeFlux) to compute stabilization terms
    Base::computeFlux( iconn, stack, [&] ( integer const ip,
                                           localIndex const k_up,
                                           localIndex const er_up,
                                           localIndex const esr_up,
                                           localIndex const ei_up,
                                           real64 const & potGrad,
                                           real64 const & phaseFlux,
                                           real64 const (&dPhaseFlux_dP)[maxStencilSize],
                                           real64 const (&dPhaseFlux_dC)[maxStencilSize][numComp] )
    {
      GEOSX_UNUSED_VAR( k_up, potGrad, phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC, er_up, esr_up, ei_up );

      // We are in the loop over phases, ip provides the current phase index.

      real64 dPresGradStab = 0.0;
      real64 tauStab = 0.0;
      integer stencilMacroElements[2]{};

      // compute potential difference MPFA-style
      for( integer i = 0; i < stack.stencilSize; ++i )
      {
        localIndex const er  = m_seri( iconn, i );
        localIndex const esr = m_sesri( iconn, i );
        localIndex const ei  = m_sei( iconn, i );

        stencilMacroElements[i] = m_macroElementIndex[er][esr][ei];

        tauStab = 9.0 * ( m_biotCoefficient[er][esr][ei] * m_biotCoefficient[er][esr][ei] )
                  / ( 32.0 * ( 10.0 * m_shearModulus[er][esr][ei] / 3.0 + m_bulkModulus[er][esr][ei] ) );

        dPresGradStab += tauStab * stack.stabTransmissibility[0][i] * (m_pres[er][esr][ei] - m_pres_n[er][esr][ei]); // jump in dp, not p
      }

      // modify stabilization flux
      // multiply dPresGrad with upwind, lagged quantities

      integer const k_up_stab = (dPresGradStab >= 0) ? 0 : 1;

      localIndex const er_up_stab   = m_seri( iconn, k_up_stab );
      localIndex const esr_up_stab  = m_sesri( iconn, k_up_stab );
      localIndex const ei_up_stab   = m_sei( iconn, k_up_stab );

      bool const areInSameMacroElement = stencilMacroElements[0] == stencilMacroElements[1];
      bool const isStabilizationActive = stencilMacroElements[0] >= 0 && stencilMacroElements[1] >= 0;
      if( isStabilizationActive && areInSameMacroElement )
      {

        for( integer ic = 0; ic < numComp; ++ic )
        {
          real64 const laggedUpwindCoef = m_phaseDens_n[er_up_stab ][esr_up_stab ][ei_up_stab ][0][ip]
                                          * m_phaseCompFrac_n[er_up_stab ][esr_up_stab ][ei_up_stab ][0][ip][ic]
                                          * m_phaseRelPerm_n[er_up_stab ][esr_up_stab ][ei_up_stab][0][ip];
          stack.stabFlux[ic] += dPresGradStab * laggedUpwindCoef;

          for( integer ke = 0; ke < stack.stencilSize; ++ke )
          {
            stack.dStabFlux_dP[ke][ic] += tauStab * stack.stabTransmissibility[0][ke] * laggedUpwindCoef;
          }
        }
      }


    } ); // end call to Base::computeFlux

    // populate local flux vector and derivatives
    for( integer ic = 0; ic < numComp; ++ic )
    {
      stack.localFlux[ic]          +=  stack.stabFlux[ic];
      stack.localFlux[numEqn + ic] += -stack.stabFlux[ic];

      for( integer ke = 0; ke < stack.stencilSize; ++ke )
      {
        localIndex const localDofIndexPres = ke * numDof;
        stack.localFluxJacobian[ic][localDofIndexPres]          +=  stack.dStabFlux_dP[ke][ic];
        stack.localFluxJacobian[numEqn + ic][localDofIndexPres] += -stack.dStabFlux_dP[ke][ic];
      }
    }

  }

protected:

  /// Views on flow properties at the previous converged time step
  ElementViewConst< arrayView1d< real64 const > > const m_pres_n;
  ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const m_phaseDens_n;
  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const m_phaseCompFrac_n;
  ElementViewConst< arrayView3d< real64 const, relperm::USD_RELPERM > > const m_phaseRelPerm_n;

  /// Views on the macroelement indices
  ElementViewConst< arrayView1d< integer const > > const m_macroElementIndex;

  /// Views on the rock/porosity properties
  ElementViewConst< arrayView1d< real64 const > > const m_bulkModulus;
  ElementViewConst< arrayView1d< real64 const > > const m_shearModulus;
  ElementViewConst< arrayView1d< real64 const > > const m_biotCoefficient;

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
    isothermalCompositionalMultiphaseBaseKernels::
      internal::kernelLaunchSelectorCompSwitch( numComps, [&] ( auto NC )
    {
      integer constexpr NUM_COMP = NC();
      integer constexpr NUM_DOF = NC()+1;

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

      using KERNEL_TYPE = FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
      typename KERNEL_TYPE::CompFlowAccessors compFlowAccessors( elemManager, solverName );
      typename KERNEL_TYPE::MultiFluidAccessors multiFluidAccessors( elemManager, solverName );
      typename KERNEL_TYPE::StabCompFlowAccessors stabCompFlowAccessors( elemManager, solverName );
      typename KERNEL_TYPE::StabMultiFluidAccessors stabMultiFluidAccessors( elemManager, solverName );
      typename KERNEL_TYPE::CapPressureAccessors capPressureAccessors( elemManager, solverName );
      typename KERNEL_TYPE::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );
      typename KERNEL_TYPE::SolidAccessors solidAccessors( elemManager, solverName );
      typename KERNEL_TYPE::PorosityAccessors porosityAccessors( elemManager, solverName );
      typename KERNEL_TYPE::RelPermAccessors relPermAccessors( elemManager, solverName );

      KERNEL_TYPE kernel( numPhases, rankOffset, hasCapPressure, stencilWrapper, dofNumberAccessor,
                          compFlowAccessors, stabCompFlowAccessors, multiFluidAccessors, stabMultiFluidAccessors,
                          capPressureAccessors, permeabilityAccessors, solidAccessors, porosityAccessors, relPermAccessors,
                          dt, localMatrix, localRhs );
      KERNEL_TYPE::template launch< POLICY >( stencilWrapper.size(), kernel );
    } );
  }
};


} // namespace stabilizedCompositionalMultiphaseFVMKernels

} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_STABILIZEDCOMPOSITIONALMULTIPHASEFVMKERNELS_HPP
