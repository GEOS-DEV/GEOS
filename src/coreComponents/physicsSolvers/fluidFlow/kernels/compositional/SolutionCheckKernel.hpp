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
 * @file SolutionCheckKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_SOLUTIONCHECKKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_SOLUTIONCHECKKERNEL_HPP

#include "LvArray/src/math.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/SolutionScalingAndCheckingKernelBase.hpp"
#include "physicsSolvers/fluidFlow/kernels/SolutionCheckKernelsHelpers.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseBaseKernels
{

/******************************** SolutionCheckKernel ********************************/

/**
 * @class SolutionCheckKernel
 * @brief Define the kernel for checking the updated solution
 */
class SolutionCheckKernel : public SolutionScalingAndCheckingKernelBase< integer >
{
public:

  using Base = SolutionScalingAndCheckingKernelBase< integer >;
  using Base::m_rankOffset;
  using Base::m_numComp;
  using Base::m_dofNumber;
  using Base::m_ghostRank;
  using Base::m_localSolution;
  using Base::m_pressure;
  using Base::m_compDens;

  /**
   * @brief Create a new kernel instance
   * @param[in] allowCompDensChopping flag to allow the component density chopping
   * @param[in] scalingFactor the scaling factor
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   * @param[in] pressure the pressure vector
   * @param[in] compDens the component density vector
   */
  SolutionCheckKernel( integer const allowCompDensChopping,
                       integer const allowNegativePressure,
                       compositionalMultiphaseUtilities::ScalingType const scalingType,
                       real64 const scalingFactor,
                       arrayView1d< real64 const > const pressure,
                       arrayView2d< real64 const, compflow::USD_COMP > const compDens,
                       arrayView1d< real64 > pressureScalingFactor,
                       arrayView1d< real64 > compDensScalingFactor,
                       globalIndex const rankOffset,
                       integer const numComp,
                       string const dofKey,
                       ElementSubRegionBase const & subRegion,
                       arrayView1d< real64 const > const localSolution,
                       ElementsReporterCollector const & negPressureIds,
                       ElementsReporterCollector const & negDensityIds,
                       ElementsReporterCollector const & negTotalDensityIds )
    : Base( rankOffset,
            numComp,
            dofKey,
            subRegion,
            localSolution,
            pressure,
            compDens,
            pressureScalingFactor,
            compDensScalingFactor ),
    m_allowCompDensChopping( allowCompDensChopping ),
    m_allowNegativePressure( allowNegativePressure ),
    m_scalingFactor( scalingFactor ),
    m_scalingType( scalingType ),
    m_negPressureIds( negPressureIds ),
    m_negDensityIds( negDensityIds ),
    m_negTotalDensityIds( negTotalDensityIds )
  {}

  /**
   * @struct KernelStats
   * @brief Kernel variables located on the stack
   */
  struct KernelStats : public Base::StackVariables
  {
    GEOS_HOST_DEVICE
    KernelStats():
      Base::StackVariables()
    {}

    KernelStats( real64 _localMinVal,
                 real64 _localMinNegPres,
                 real64 _localMinNegDens,
                 real64 _localMinNegTotalDens )
      :
      Base::StackVariables( _localMinVal ),
      localMinNegPres( _localMinNegPres ),
      localMinNegDens( _localMinNegDens ),
      localMinNegTotalDens( _localMinNegTotalDens )
    { }

    real64 localMinNegPres;
    real64 localMinNegDens;
    real64 localMinNegTotalDens;
  };

  /**
   * @struct StackVariables
   * @brief Kernel variables located on the stack
   */
  struct StackVariables : public KernelStats
  {
    GEOS_HOST_DEVICE
    StackVariables():
      KernelStats()
    { }

    localIndex localNumNegPres; // 0 or 1
    localIndex localNumNegDens; // 0 -> num comp
    localIndex localNumNegTotalDens; // 0 or 1
  };

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[inout] kernelComponent the kernel component providing access to the compute function
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static KernelStats
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    using reducePolicy = ReducePolicy< POLICY >;
    using atomicPolicy = AtomicPolicy< POLICY >;

    RAJA::ReduceMin< reducePolicy, integer > globalMinVal( 1 );

    RAJA::ReduceMin< reducePolicy, real64 > minPres( 0.0 );
    RAJA::ReduceMin< reducePolicy, real64 > minDens( 0.0 );
    RAJA::ReduceMin< reducePolicy, real64 > minTotalDens( 0.0 );

    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( kernelComponent.ghostRank( ei ) >= 0 )
      {
        return;
      }

      StackVariables stack;
      kernelComponent.setup( ei, stack );
      kernelComponent.compute( ei, stack );

      globalMinVal.min( stack.localMinVal );

      minPres.min( stack.localMinNegPres );
      minDens.min( stack.localMinNegDens );
      minTotalDens.min( stack.localMinNegTotalDens );

      if( stack.localNumNegPres > 0 )
        kernelComponent.m_negPressureIds.collectElement( atomicPolicy{}, { ei, stack.localMinNegPres } );

      if( stack.localNumNegDens > 0 )
        kernelComponent.m_negDensityIds.collectElement( atomicPolicy{}, { ei, stack.localMinNegDens } );

      if( stack.localNumNegTotalDens > 0 )
        kernelComponent.m_negDensityIds.collectElement( atomicPolicy{}, { ei, stack.localMinNegTotalDens } );
    } );

    return KernelStats( globalMinVal.get(),
                        minPres.get(),
                        minDens.get(),
                        minTotalDens.get() );
  }

  GEOS_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    Base::setup( ei, stack );

    stack.localMinNegPres = 0.0;
    stack.localMinNegDens = 0.0;
    stack.localMinNegTotalDens = 0.0;

    stack.localNumNegPres = 0;
    stack.localNumNegDens = 0;
    stack.localNumNegTotalDens = 0;
  }

  /**
   * @brief Compute the local value
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void compute( localIndex const ei,
                StackVariables & stack ) const
  {
    computeSolutionCheck( ei, stack );
  }

  /**
   * @brief Compute the local value of the check
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] kernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeSolutionCheck( localIndex const ei,
                             StackVariables & stack,
                             FUNC && kernelOp = NoOpFunc{} ) const
  {
    bool const localScaling = m_scalingType == compositionalMultiphaseUtilities::ScalingType::Local;

    real64 const newPres = m_pressure[ei] + (localScaling ? m_pressureScalingFactor[ei] : m_scalingFactor) * m_localSolution[stack.localRow];
    if( newPres < 0.0 && !m_allowNegativePressure )
      stack.localMinVal = 0;

    if( newPres <= 0.0 )
    {
      stack.localNumNegPres += 1;

      if( newPres < stack.localMinNegPres )
        stack.localMinNegPres = newPres;
    }

    // if component density chopping is not allowed, the time step fails if a component density is negative
    // otherwise, we just check that the total density is positive, and negative component densities
    // will be chopped (i.e., set to zero) in ApplySystemSolution)
    real64 totalDens = 0.0;
    if( !m_allowCompDensChopping )
    {
      for( integer ic = 0; ic < m_numComp; ++ic )
      {
        real64 const newDens = m_compDens[ei][ic] + (localScaling ? m_compDensScalingFactor[ei] : m_scalingFactor) * m_localSolution[stack.localRow + ic + 1];
        totalDens += newDens;

        // we invalidate timestep if new density is negative, and we report density as too low if negative or equal to zero
        bool const isNewDensNegative = newDens < 0.0;
        bool const isNewDensTooLow = newDens <= 0.0;
        stack.localMinVal = isNewDensNegative ? 0 : stack.localMinVal;
        stack.localNumNegDens += isNewDensTooLow;
        stack.localMinNegDens = isNewDensTooLow ?
                                LvArray::math::min( stack.localMinNegDens, newDens ) :
                                stack.localMinNegDens;
      }
    }
    else
    {
      for( integer ic = 0; ic < m_numComp; ++ic )
      {
        real64 const newDens = m_compDens[ei][ic] + (localScaling ? m_compDensScalingFactor[ei] : m_scalingFactor) * m_localSolution[stack.localRow + ic + 1];
        totalDens += LvArray::math::max( newDens, 0.0 );

        // we never invalidate timestep as negatives will be chopped later in ApplySystemSolution(), and we report density as too low if
        // negative or equal to zero
        bool const isNewDensTooLow = newDens <= 0.0;
        stack.localNumNegDens += isNewDensTooLow;
        stack.localMinNegDens = isNewDensTooLow ?
                                LvArray::math::min( stack.localMinNegDens, newDens ) :
                                stack.localMinNegDens;
      }
    }
    {
      bool const isNewTotalDensTooLow = totalDens <= 0.0;
      stack.localNumNegTotalDens += isNewTotalDensTooLow;
      stack.localMinNegTotalDens = isNewTotalDensTooLow ?
                                   LvArray::math::min( stack.localMinNegTotalDens, totalDens ) :
                                   stack.localMinNegTotalDens;
    }

    kernelOp();
  }

protected:

  /// flag to allow the component density chopping
  integer const m_allowCompDensChopping;

  /// flag to allow negative pressure values
  integer const m_allowNegativePressure;

  /// scaling factor
  real64 const m_scalingFactor;

  /// scaling type (global or local)
  compositionalMultiphaseUtilities::ScalingType const m_scalingType;

  ElementsReporterCollector const m_negPressureIds;

  ElementsReporterCollector const m_negDensityIds;

  ElementsReporterCollector const m_negTotalDensityIds;

};

/**
 * @class SolutionCheckKernelFactory
 */
class SolutionCheckKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] allowCompDensChopping flag to allow the component density chopping
   * @param[in] scalingFactor the scaling factor
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   */
  template< typename POLICY >
  static SolutionCheckKernel::KernelStats
  createAndLaunch( integer const allowCompDensChopping,
                   integer const allowNegativePressure,
                   compositionalMultiphaseUtilities::ScalingType const scalingType,
                   real64 const scalingFactor,
                   arrayView1d< real64 const > const pressure,
                   arrayView2d< real64 const, compflow::USD_COMP > const compDens,
                   arrayView1d< real64 > pressureScalingFactor,
                   arrayView1d< real64 > compDensScalingFactor,
                   globalIndex const rankOffset,
                   integer const numComp,
                   string const dofKey,
                   ElementSubRegionBase & subRegion,
                   arrayView1d< real64 const > const localSolution,
                   ElementsReporterCollector const & negPressureIds,
                   ElementsReporterCollector const & negDensityIds,
                   ElementsReporterCollector const & negTotalDensityIds )
  {
    SolutionCheckKernel kernel( allowCompDensChopping, allowNegativePressure, scalingType, scalingFactor,
                                pressure, compDens, pressureScalingFactor, compDensScalingFactor, rankOffset,
                                numComp, dofKey, subRegion, localSolution, negPressureIds, negDensityIds,
                                negTotalDensityIds );
    return SolutionCheckKernel::launch< POLICY >( subRegion.size(), kernel );
  }

};

} // namespace isothermalCompositionalMultiphaseBaseKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_SOLUTIONCHECKKERNEL_HPP
