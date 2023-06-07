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

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SOLUTIONCHECKKERNEL_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SOLUTIONCHECKKERNEL_HPP

#include "ScalingAndCheckingSystemSolutionKernelBase.hpp"


namespace geos
{

namespace isothermalCompositionalMultiphaseBaseKernels
{

/******************************** SolutionCheckKernel ********************************/

/**
 * @class SolutionCheckKernel
 * @brief Define the kernel for checking the updated solution
 */
class SolutionCheckKernel : public ScalingAndCheckingSystemSolutionKernelBase< geos::integer >
{
public:

  using Base = ScalingAndCheckingSystemSolutionKernelBase< geos::integer >;
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
  SolutionCheckKernel( geos::integer const allowCompDensChopping,
                       geos::real64 const scalingFactor,
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
    m_allowCompDensChopping( allowCompDensChopping ),
    m_scalingFactor( scalingFactor )
  { }

  GEOS_HOST_DEVICE
  virtual void compute( localIndex const ei,
                        StackVariables & stack ) const override
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
    real64 const newPres = m_pressure[ei] + m_scalingFactor * m_localSolution[stack.localRow];
    if( newPres < 0 )
    {
      stack.localMinVal = 0;
    }

    // if component density chopping is not allowed, the time step fails if a component density is negative
    // otherwise, we just check that the total density is positive, and negative component densities
    // will be chopped (i.e., set to zero) in ApplySystemSolution)
    if( !m_allowCompDensChopping )
    {
      for( integer ic = 0; ic < m_numComp; ++ic )
      {
        real64 const newDens = m_compDens[ei][ic] + m_scalingFactor * m_localSolution[stack.localRow + ic + 1];
        if( newDens < 0 )
        {
          stack.localMinVal = 0;
        }
      }
    }
    else
    {
      real64 totalDens = 0.0;
      for( integer ic = 0; ic < m_numComp; ++ic )
      {
        real64 const newDens = m_compDens[ei][ic] + m_scalingFactor * m_localSolution[stack.localRow + ic + 1];
        totalDens += ( newDens > 0.0 ) ? newDens : 0.0;
      }
      if( totalDens < 0 )
      {
        stack.localMinVal = 0;
      }
    }

    kernelOp();
  }

protected:

  /// flag to allow the component density chopping
  integer const m_allowCompDensChopping;

  /// scaling factor
  real64 const m_scalingFactor;

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
    arrayView2d< real64 const, compflow::USD_COMP > const compDens =
      subRegion.getField< fields::flow::globalCompDensity >();
    SolutionCheckKernel kernel( allowCompDensChopping, scalingFactor, rankOffset,
                                numComp, dofKey, subRegion, localSolution, pressure, compDens );
    return SolutionCheckKernel::launch< POLICY >( subRegion.size(), kernel );
  }

};

}

}

#endif //GEOSX_SOLUTIONCHECKKERNEL_HPP
