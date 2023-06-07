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
 * @file SolutionCheckKernel.hpp
 */

#ifndef GEOSX_SOLUTIONCHECKKERNEL_HPP
#define GEOSX_SOLUTIONCHECKKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/compositionalMultiphase/SolutionCheckKernel.hpp"


namespace geos
{

namespace thermalCompositionalMultiphaseBaseKernels
{

/******************************** SolutionCheckKernel ********************************/

/**
 * @class SolutionCheckKernel
 * @brief Define the kernel for checking the updated solution
 */
class SolutionCheckKernel : public isothermalCompositionalMultiphaseBaseKernels::SolutionCheckKernel
{
public:

  using Base = isothermalCompositionalMultiphaseBaseKernels::SolutionCheckKernel;
  using Base::m_numComp;
  using Base::m_localSolution;
  using Base::m_scalingFactor;

  static geos::real64 constexpr minTemperature = 273.15;

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
   * @param[in] temperature the temperature vector
   * @param[in] compDens the component density vector
   */
  SolutionCheckKernel( geos::integer const allowCompDensChopping,
                       geos::real64 const scalingFactor,
                       geos::globalIndex const rankOffset,
                       geos::integer const numComp,
                       geos::string const dofKey,
                       geos::ElementSubRegionBase const & subRegion,
                       geos::arrayView1d< geos::real64 const > const localSolution,
                       geos::arrayView1d< geos::real64 const > const pressure,
                       geos::arrayView1d< geos::real64 const > const temperature,
                       arrayView2d< geos::real64 const, geos::compflow::USD_COMP > const compDens )
    : Base( allowCompDensChopping,
            scalingFactor,
            rankOffset,
            numComp,
            dofKey,
            subRegion,
            localSolution,
            pressure,
            compDens ),
    m_temperature( temperature )
  { }

  /**
   * @brief Compute the local value of the solution check
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeSolutionCheck( localIndex const ei,
                             StackVariables & stack ) const
  {
    Base::computeSolutionCheck( ei, stack, [&]()
    {
      // compute the change in temperature
      real64 const newTemp = m_temperature[ei] + m_scalingFactor * m_localSolution[stack.localRow + m_numComp + 1];
      if( newTemp < minTemperature )
      {
        stack.localMinVal = 0;
      }
    } );
  }

protected:

  /// View on the primary variables
  arrayView1d< real64 const > const m_temperature;

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
  static integer
  createAndLaunch( integer const allowCompDensChopping,
                   real64 const scalingFactor,
                   globalIndex const rankOffset,
                   integer const numComp,
                   string const dofKey,
                   ElementSubRegionBase const & subRegion,
                   arrayView1d< real64 const > const localSolution )
  {
    arrayView1d< real64 const > const pressure =
      subRegion.getField< fields::flow::pressure >();
    arrayView1d< real64 const > const temperature =
      subRegion.getField< fields::flow::temperature >();
    arrayView2d< real64 const, compflow::USD_COMP > const compDens =
      subRegion.getField< fields::flow::globalCompDensity >();
    SolutionCheckKernel kernel( allowCompDensChopping, scalingFactor,
                                rankOffset, numComp, dofKey, subRegion, localSolution,
                                pressure, temperature, compDens );
    return SolutionCheckKernel::launch< POLICY >( subRegion.size(), kernel );
  }

};

}

}

#endif //GEOSX_SOLUTIONCHECKKERNEL_HPP
