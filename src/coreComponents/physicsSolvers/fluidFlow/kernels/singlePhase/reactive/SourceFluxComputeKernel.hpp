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
 * @file SourceFluxComputeKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_SOURCEFLUXCOMPUTEKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_SOURCEFLUXCOMPUTEKERNELS_HPP

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/reactivefluid/ReactiveSinglePhaseFluid.hpp"
#include "constitutive/fluid/reactivefluid/ReactiveFluidLayouts.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidUtils.hpp"
#include "codingUtilities/Utilities.hpp"

#include "physicsSolvers/fluidFlow/kernels/singlePhase/reactive/KernelLaunchSelectors.hpp"

namespace geos
{

namespace singlePhaseReactiveBaseKernels
{

/******************************** SourceFluxComputeKernel ********************************/

/**
 * @class SourceFluxComputeKernel
 * @brief Define the interface for the assembly kernel in charge of source flux
 */
template< integer NUM_DOF, integer NUM_SPECIES, typename BASE_FLUID_TYPE >
class SourceFluxComputeKernel
{

public:

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;

  /// Compute time value for the number of equations
  static constexpr integer numEqn = NUM_DOF;

  /// Compile time value for the number of primary species
  static constexpr integer numSpecies = NUM_SPECIES;

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
    m_rankOffset( rankOffset ),
    m_dofNumber( dofNumber ),
    m_elemGhostRank( elemGhostRank ),
    m_rhsContributionArrayView( rhsContributionArrayView ),
    m_sizeScalingFactor( sizeScalingFactor ),
    m_solventDensity( fluid.solventDensity() ),
    m_primarySpeciesAggregateConcentration( fluid.primarySpeciesAggregateConcentration() ),
    m_dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations( fluid.dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations() ),
    m_density( fluid.density() ),
    m_dDensity( fluid.dDensity() ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs ),
    m_massProd( massProd )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {
public:

    /// Index of the local row corresponding to this element
    localIndex massRowIndex = -1;

    /// Index of the matrix row/column corresponding to the dof in this element
    globalIndex dofIndices[numDof]{};

    /// Storage for the element local residual vector for species rows
    real64 localSpeciesRhs[numSpecies]{};

    /// Storage for the element local Jacobian matrix for species rows
    real64 localSpeciesJacobian[numSpecies][numDof]{};

    real64 totalInflowMass = 0.0;

  };

  /**
   * @brief Getter for the ghost rank of an element
   * @param[in] ei the element index
   * @return the ghost rank of the element
   */
  GEOS_HOST_DEVICE
  integer elemGhostRank( localIndex const ei ) const
  { return m_elemGhostRank( ei ); }

  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] ei the element index
   * @param[in] a the target set index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const ei,
              localIndex const a,
              StackVariables & stack ) const
  {
    // set row index and degrees of freedom indices for this element
    stack.massRowIndex = m_dofNumber[ei] - m_rankOffset;
    for( integer idof = 0; idof < numDof; ++idof )
    {
      stack.dofIndices[idof] = m_dofNumber[ei] + idof - m_rankOffset;
    }

    stack.totalInflowMass = m_rhsContributionArrayView[a];
  }

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
    real64 const scaledInflowMass = stack.totalInflowMass / m_sizeScalingFactor;

    for( integer i = 0; i < numSpecies; ++i )
    {
      stack.localSpeciesRhs[i] += m_primarySpeciesAggregateConcentration[ei][0][i] * m_solventDensity / m_density[ei][0] * scaledInflowMass;
      stack.localSpeciesJacobian[i][0] += -m_primarySpeciesAggregateConcentration[ei][0][i] * m_solventDensity * m_dDensity[ei][0][DerivOffset::dP] / (m_density[ei][0] * m_density[ei][0]) *
                                          scaledInflowMass;

      for( integer j = 0; j < numSpecies; ++j )
      {
        stack.localSpeciesJacobian[i][j+numDof-numSpecies] += m_dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations[ei][0][i][j] * m_solventDensity / m_density[ei][0] *
                                                              scaledInflowMass;
      }
    }
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void complete( localIndex const GEOS_UNUSED_PARAM( ei ),
                 StackVariables & stack ) const
  {
    m_massProd += stack.totalInflowMass / m_sizeScalingFactor;
    m_localRhs[stack.massRowIndex] += stack.totalInflowMass / m_sizeScalingFactor;

    if( stack.totalInflowMass > 0.0 )
    {
      for( integer i = 0; i < numSpecies; ++i )
      {
        globalIndex const speciesRowBeginIndex = stack.massRowIndex + numEqn - numSpecies;
        m_localRhs[speciesRowBeginIndex + i] += stack.localSpeciesRhs[i];

        // add contribution to global residual and jacobian (no need for atomics here)
        m_localMatrix.template addToRow< serialAtomic >( speciesRowBeginIndex+i,
                                                         stack.dofIndices,
                                                         stack.localSpeciesJacobian[i],
                                                         numDof );
      }
    }
  }

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] targetSet the target set for source flux
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack variables
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( SortedArrayView< localIndex const > const targetSet,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;

    forAll< POLICY >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      // we need to filter out ghosts here, because targetSet may contain them
      localIndex const ei = targetSet[a];

      if( kernelComponent.elemGhostRank( ei ) >= 0 )
      {
        return;
      }

      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( ei, a, stack );
      kernelComponent.computeSourceFlux( ei, stack );
      kernelComponent.complete( ei, stack );
    } );
  }

protected:

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// View on the dof numbers
  arrayView1d< globalIndex const > const m_dofNumber;

  /// View on the ghost ranks
  arrayView1d< integer const > const m_elemGhostRank;

  /// View on the rhs contribution
  arrayView1d< real64 const > const m_rhsContributionArrayView;
  /// size scaling factor
  real64 const m_sizeScalingFactor;

  /// Solvent density [kg/m³] used to convert molality [mol/kg] to molarity [mol/m³]
  real64 const m_solventDensity;

  // View on the total concentration of ions that contain the primary species
  arrayView3d< real64 const, constitutive::reactivefluid::USD_SPECIES > const m_primarySpeciesAggregateConcentration;
  // View on the derivatives of total ion concentration for the primary species wrt log of primary species concentration
  arrayView4d< real64 const, constitutive::reactivefluid::USD_SPECIES_DC > const m_dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations;

  // View on the fluid density
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const m_density;
  // View on the derivatives of fluid density
  arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const m_dDensity;

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

  RAJA::ReduceSum< parallelDeviceReduce, real64 > m_massProd;

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
    internal::kernelLaunchSelectorCompSwitch( numSpecies, [&] ( auto NS )
    {
      integer constexpr NUM_SPECIES = NS();
      integer constexpr NUM_DOF = 1+NS();

      SourceFluxComputeKernel< NUM_DOF, NUM_SPECIES, BASE_FLUID_TYPE > kernel( rankOffset, dofNumber, elemGhostRank, rhsContributionArrayView, sizeScalingFactor, fluid, localMatrix, localRhs,
                                                                               massProd );
      SourceFluxComputeKernel< NUM_DOF, NUM_SPECIES, BASE_FLUID_TYPE >::template launch< POLICY >( targetSet, kernel );
    } );
  }
};
} // namespace singlePhaseReactiveBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_SOURCEFLUXCOMPUTEKERNELS_HPP
