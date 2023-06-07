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

#ifndef GEOSX_SCALINGFORSYSTEMSOLUTIONKERNEL_HPP
#define GEOSX_SCALINGFORSYSTEMSOLUTIONKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/compositionalMultiphase/ScalingForSystemSolutionKernel.hpp"


namespace geos
{

namespace thermalCompositionalMultiphaseBaseKernels
{

/******************************** ScalingForSystemSolutionKernel ********************************/

/**
 * @class ScalingForSystemSolutionKernel
 * @brief Define the kernel for scaling the Newton update
 */
class ScalingForSystemSolutionKernel : public isothermalCompositionalMultiphaseBaseKernels::ScalingForSystemSolutionKernel
{
public:

  using Base = isothermalCompositionalMultiphaseBaseKernels::ScalingForSystemSolutionKernel;
  using Base::m_numComp;
  using Base::m_localSolution;

  /**
   * @brief Create a new kernel instance
   * @param[in] maxRelativePresChange the max allowed relative pressure change
   * @param[in] maxRelativeTempChange the max allowed relative temperature change
   * @param[in] maxCompFracChange the max allowed comp fraction change
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   * @param[in] pressure the pressure vector
   * @param[in] temperature the temperature vector
   * @param[in] compDens the component density vector
   */
  ScalingForSystemSolutionKernel( geos::real64 const maxRelativePresChange,
                                  geos::real64 const maxRelativeTempChange,
                                  geos::real64 const maxCompFracChange,
                                  geos::globalIndex const rankOffset,
                                  geos::integer const numComp,
                                  geos::string const dofKey,
                                  geos::ElementSubRegionBase const & subRegion,
                                  geos::arrayView1d< geos::real64 const > const localSolution,
                                  geos::arrayView1d< geos::real64 const > const pressure,
                                  geos::arrayView1d< geos::real64 const > const temperature,
                                  arrayView2d< geos::real64 const, geos::compflow::USD_COMP > const compDens )
    : Base( maxRelativePresChange,
            maxCompFracChange,
            rankOffset,
            numComp,
            dofKey,
            subRegion,
            localSolution,
            pressure,
            compDens ),
    m_maxRelativeTempChange( maxRelativeTempChange ),
    m_temperature( temperature )
  { }

  /**
   * @brief Compute the local value of the scaling factor
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeScalingFactor( localIndex const ei,
                             StackVariables & stack ) const
  {
    real64 constexpr eps = isothermalCompositionalMultiphaseBaseKernels::minDensForDivision;

    Base::computeScalingFactor( ei, stack, [&]()
    {
      // compute the change in temperature
      real64 const temp = m_temperature[ei];
      if( temp > eps )
      {
        real64 const absTempChange = LvArray::math::abs( m_localSolution[stack.localRow + m_numComp + 1] );
        real64 const relativeTempChange = absTempChange / temp;
        if( relativeTempChange > m_maxRelativeTempChange )
        {
          real64 const tempScalingFactor = m_maxRelativeTempChange / relativeTempChange;
          if( stack.localMinVal > tempScalingFactor )
          {
            stack.localMinVal = tempScalingFactor;
          }
        }
      }
    } );
  }

protected:

  /// Max allowed changes in primary variables
  real64 const m_maxRelativeTempChange;

  /// View on the primary variables
  arrayView1d< real64 const > const m_temperature;

};

/**
 * @class ScalingForSystemSolutionKernelFactory
 */
class ScalingForSystemSolutionKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] maxRelativePresChange the max allowed relative pressure change
   * @param[in] maxRelativeTempChange the max allowed relative temperature change
   * @param[in] maxCompFracChange the max allowed comp fraction change
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   */
  template< typename POLICY >
  static real64
  createAndLaunch( real64 const maxRelativePresChange,
                   real64 const maxRelativeTempChange,
                   real64 const maxCompFracChange,
                   globalIndex const rankOffset,
                   integer const numComp,
                   string const dofKey,
                   ElementSubRegionBase const & subRegion,
                   arrayView1d< real64 const > const localSolution )
  {
    arrayView1d< real64 const > const pressure = subRegion.getField< fields::flow::pressure >();
    arrayView1d< real64 const > const temperature = subRegion.getField< fields::flow::temperature >();
    arrayView2d< real64 const, compflow::USD_COMP > const compDens = subRegion.getField< fields::flow::globalCompDensity >();
    ScalingForSystemSolutionKernel kernel( maxRelativePresChange, maxRelativeTempChange, maxCompFracChange,
                                           rankOffset, numComp, dofKey, subRegion, localSolution,
                                           pressure, temperature, compDens );
    return ScalingForSystemSolutionKernel::launch< POLICY >( subRegion.size(), kernel );
  }

};

}

}

#endif //GEOSX_SCALINGFORSYSTEMSOLUTIONKERNEL_HPP
