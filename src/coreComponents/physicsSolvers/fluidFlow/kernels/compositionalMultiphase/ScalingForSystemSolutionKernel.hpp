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
 * @file ScalingForSystemSolutionKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SCALINGFORSYSTEMSOLUTIONKERNEL_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SCALINGFORSYSTEMSOLUTIONKERNEL_HPP

#include "ScalingAndCheckingSystemSolutionKernelBase.hpp"


namespace geos
{

namespace isothermalCompositionalMultiphaseBaseKernels
{

/******************************** ScalingForSystemSolutionKernel ********************************/

/**
 * @class ScalingForSystemSolutionKernel
 * @brief Define the kernel for scaling the Newton update
 */
class ScalingForSystemSolutionKernel : public ScalingAndCheckingSystemSolutionKernelBase< geos::real64 >
{
public:

  using Base = ScalingAndCheckingSystemSolutionKernelBase< geos::real64 >;
  using Base::m_rankOffset;
  using Base::m_numComp;
  using Base::m_dofNumber;
  using Base::m_ghostRank;
  using Base::m_localSolution;
  using Base::m_pressure;
  using Base::m_compDens;

  /**
   * @brief Create a new kernel instance
   * @param[in] maxRelativePresChange the max allowed relative pressure change
   * @param[in] maxCompFracChange the max allowed comp fraction change
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   * @param[in] pressure the pressure vector
   * @param[in] compDens the component density vector
   */
  ScalingForSystemSolutionKernel( geos::real64 const maxRelativePresChange,
                                  geos::real64 const maxCompFracChange,
                                  geos::globalIndex const rankOffset,
                                  geos::integer const numComp,
                                  geos::string const dofKey,
                                  geos::ElementSubRegionBase const & subRegion,
                                  geos::arrayView1d< geos::real64 const > const localSolution,
                                  geos::arrayView1d< geos::real64 const > const pressure,
                                  arrayView2d< geos::real64 const, geos::compflow::USD_COMP > const compDens )
    : Base( rankOffset,
            numComp,
            dofKey,
            subRegion,
            localSolution,
            pressure,
            compDens ),
    m_maxRelativePresChange( maxRelativePresChange ),
    m_maxCompFracChange( maxCompFracChange )
  { }

  GEOS_HOST_DEVICE
  virtual void compute( localIndex const ei,
                        StackVariables & stack ) const override
  {
    computeScalingFactor( ei, stack );
  }

  /**
   * @brief Compute the local value of the scaling factor
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] kernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeScalingFactor( localIndex const ei,
                             StackVariables & stack,
                             FUNC && kernelOp = NoOpFunc{} ) const
  {
    real64 constexpr eps = minDensForDivision;

    // compute the change in pressure
    real64 const pres = m_pressure[ei];
    if( pres > eps )
    {
      real64 const absPresChange = LvArray::math::abs( m_localSolution[stack.localRow] );
      real64 const relativePresChange = absPresChange / pres;
      if( relativePresChange > m_maxRelativePresChange )
      {
        real64 const presScalingFactor = m_maxRelativePresChange / relativePresChange;

        if( stack.localMinVal > presScalingFactor )
        {
          stack.localMinVal = presScalingFactor;
        }
      }
    }

    real64 prevTotalDens = 0;
    for( integer ic = 0; ic < m_numComp; ++ic )
    {
      prevTotalDens += m_compDens[ei][ic];
    }

    // compute the change in component densities and component fractions
    for( integer ic = 0; ic < m_numComp; ++ic )
    {
      // compute scaling factor based on relative change in component densities
      real64 const absCompDensChange = LvArray::math::abs( m_localSolution[stack.localRow + ic + 1] );
      real64 const maxAbsCompDensChange = m_maxCompFracChange * prevTotalDens;

      // This actually checks the change in component fraction, using a lagged total density
      // Indeed we can rewrite the following check as:
      //    | prevCompDens / prevTotalDens - newCompDens / prevTotalDens | > maxCompFracChange
      // Note that the total density in the second term is lagged (i.e, we use prevTotalDens)
      // because I found it more robust than using directly newTotalDens (which can vary also
      // wildly when the compDens change is large)
      if( absCompDensChange > maxAbsCompDensChange && absCompDensChange > eps )
      {
        real64 const compScalingFactor = maxAbsCompDensChange / absCompDensChange;
        if( stack.localMinVal > compScalingFactor )
        {
          stack.localMinVal = compScalingFactor;
        }
      }
    }

    // compute the scaling factor for other vars, such as temperature
    kernelOp();
  }

protected:

  /// Max allowed changes in primary variables
  real64 const m_maxRelativePresChange;
  real64 const m_maxCompFracChange;

};

/**
 * @class ScalingForSystemSolutionKernelFactory
 */
class ScalingForSystemSolutionKernelFactory
{
public:

  /*
   * @brief Create and launch the kernel computing the scaling factor
   * @tparam POLICY the kernel policy
   * @param[in] maxRelativePresChange the max allowed relative pressure change
   * @param[in] maxCompFracChange the max allowed comp fraction change
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   * @return the scaling factor
   */
  template< typename POLICY >
  static real64
  createAndLaunch( real64 const maxRelativePresChange,
                   real64 const maxCompFracChange,
                   globalIndex const rankOffset,
                   integer const numComp,
                   string const dofKey,
                   ElementSubRegionBase const & subRegion,
                   arrayView1d< real64 const > const localSolution )
  {
    arrayView1d< real64 const > const pressure =
      subRegion.getField< fields::flow::pressure >();
    arrayView2d< real64 const, compflow::USD_COMP > const compDens =
      subRegion.getField< fields::flow::globalCompDensity >();
    ScalingForSystemSolutionKernel kernel( maxRelativePresChange, maxCompFracChange, rankOffset,
                                           numComp, dofKey, subRegion, localSolution, pressure, compDens );
    return ScalingForSystemSolutionKernel::launch< POLICY >( subRegion.size(), kernel );
  }
};

}

}

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SCALINGFORSYSTEMSOLUTIONKERNEL_HPP
