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
 * @file ThermalSourceFluxComputeKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_THERMALSOURCEFLUXCOMPUTEKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_THERMALSOURCEFLUXCOMPUTEKERNELS_HPP

#include "physicsSolvers/fluidFlow/kernels/singlePhase/reactive/SourceFluxComputeKernel.hpp"

namespace geos
{

namespace thermalSinglePhaseReactiveBaseKernels
{

/******************************** SourceFluxComputeKernel ********************************/

/**
 * @class SourceFluxComputeKernel
 * @brief Define the interface for the assembly kernel in charge of source flux
 */
template< integer NUM_DOF, integer NUM_SPECIES, typename BASE_FLUID_TYPE >
class SourceFluxComputeKernel : public singlePhaseReactiveBaseKernels::SourceFluxComputeKernel< NUM_DOF, NUM_SPECIES, BASE_FLUID_TYPE >
{

public:

  using Base = singlePhaseReactiveBaseKernels::SourceFluxComputeKernel< NUM_DOF, NUM_SPECIES, BASE_FLUID_TYPE >;
  using Base::numSpecies;
  using Base::numDof;
  using Base::numEqn;
  using Base::m_sizeScalingFactor;
  using Base::m_solventDensity;
  using Base::m_primarySpeciesAggregateConcentration;
  using Base::m_density;
  using Base::m_dDensity;
  using Base::m_localMatrix;
  using Base::m_localRhs;

  using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< 1 >;

  SourceFluxComputeKernel( globalIndex const rankOffset,
                           arrayView1d< globalIndex const > const dofNumber,
                           arrayView1d< integer const > const elemGhostRank,
                           arrayView1d< real64 const > const rhsContributionArrayView,
                           real64 const sizeScalingFactor,
                           constitutive::reactivefluid::ReactiveSinglePhaseFluid< BASE_FLUID_TYPE > const & fluid,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs,
                           RAJA::ReduceSum< parallelDeviceReduce, real64 > massProd )
    :
    Base( rankOffset,
          dofNumber,
          elemGhostRank,
          rhsContributionArrayView,
          sizeScalingFactor,
          fluid,
          localMatrix,
          localRhs,
          massProd ),
    m_enthalpy( fluid.enthalpy() ),
    m_dEnthalpy( fluid.dEnthalpy() )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables : Base::StackVariables
  {
public:

    GEOS_HOST_DEVICE
    StackVariables()
      : Base::StackVariables()
    {}

    using Base::StackVariables::massRowIndex;
    using Base::StackVariables::dofIndices;
    using Base::StackVariables::localSpeciesJacobian;
    using Base::StackVariables::totalInflowMass;

    /// Storage for the element local residual vector for energy row
    real64 localEnergyRhs=0.0;

    /// Storage for the element local Jacobian matrix for energy row
    real64 localEnergyJacobian[numDof]{};

  };

  /**
   * @brief Compute the local source flux contributions to the residual and Jacobian
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] kernelOp the function used to customize the kernel
   */
  GEOS_HOST_DEVICE
  void computeSourceFlux( localIndex const ei,
                          StackVariables & stack ) const
  {
    Base::computeSourceFlux( ei, stack );

    real64 const scaledInflowMass = stack.totalInflowMass / m_sizeScalingFactor;

    stack.localEnergyRhs += m_enthalpy[ei][0] * scaledInflowMass;
    stack.localEnergyJacobian[0] = scaledInflowMass * m_dEnthalpy[ei][0][DerivOffset::dP];
    stack.localEnergyJacobian[numDof-numSpecies-1] = scaledInflowMass * m_dEnthalpy[ei][0][DerivOffset::dT];

    for( integer i = 0; i < numSpecies; ++i )
    {
      stack.localSpeciesJacobian[i][numDof-numSpecies-1] += -m_primarySpeciesAggregateConcentration[ei][0][i] * m_solventDensity * m_dDensity[ei][0][DerivOffset::dT] /
                                                            (m_density[ei][0] * m_density[ei][0]) *
                                                            scaledInflowMass;
    }
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void complete( localIndex const ei,
                 StackVariables & stack ) const
  {
    Base::complete( ei, stack );

    if( stack.totalInflowMass > 0.0 )
    {
      globalIndex const energyRowIndex = stack.massRowIndex + 1;
      m_localRhs[energyRowIndex] += stack.localEnergyRhs;

      m_localMatrix.template addToRow< serialAtomic >( energyRowIndex,
                                                       stack.dofIndices,
                                                       stack.localEnergyJacobian,
                                                       numDof );
    }
  }

protected:

  /// Views on enthalpies
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const m_enthalpy;
  arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const m_dEnthalpy;

};

/**
 * @class SourceFluxComputeKernelFactory
 */
class SourceFluxComputeKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numSpecies the number of primary species
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofNumber the degress of freedom numbers
   * @param[in] elemGhostRank the array of element ghost rank
   * @param[in] targetSet the target set array
   * @param[in] rhsContributionArrayView the rhs contribution array
   * @param[in] sizeScalingFactor the size scaling factor
   * @param[in] fluid the fluid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[inout] massProd the total mass produced
   */
  template< typename POLICY, typename BASE_FLUID_TYPE >
  static void
  createAndLaunch( integer const numSpecies,
                   globalIndex const rankOffset,
                   arrayView1d< globalIndex const > const dofNumber,
                   arrayView1d< integer const > const elemGhostRank,
                   SortedArrayView< localIndex const > const targetSet,
                   arrayView1d< real64 const > const rhsContributionArrayView,
                   real64 const sizeScalingFactor,
                   constitutive::reactivefluid::ReactiveSinglePhaseFluid< BASE_FLUID_TYPE > const & fluid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs,
                   RAJA::ReduceSum< parallelDeviceReduce, real64 > massProd )
  {
    singlePhaseReactiveBaseKernels::
      internal::kernelLaunchSelectorCompSwitch( numSpecies, [&] ( auto NS )
    {
      integer constexpr NUM_SPECIES = NS();
      integer constexpr NUM_DOF = 2+NS();

      SourceFluxComputeKernel< NUM_DOF, NUM_SPECIES, BASE_FLUID_TYPE > kernel( rankOffset, dofNumber, elemGhostRank, rhsContributionArrayView, sizeScalingFactor, fluid, localMatrix, localRhs,
                                                                               massProd );
      SourceFluxComputeKernel< NUM_DOF, NUM_SPECIES, BASE_FLUID_TYPE >::template launch< POLICY >( targetSet, kernel );
    } );
  }
};
} // namespace thermalSinglePhaseReactiveBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_THERMALSOURCEFLUXCOMPUTEKERNELS_HPP
