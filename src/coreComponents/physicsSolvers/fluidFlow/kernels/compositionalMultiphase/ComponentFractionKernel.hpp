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
 * @file ComponentFractionKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPONENTFRACTIONKERNEL_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPONENTFRACTIONKERNEL_HPP

#include "PropertyKernelBase.hpp"


namespace geos
{

namespace isothermalCompositionalMultiphaseBaseKernels
{

/******************************** ComponentFractionKernel ********************************/

/**
 * @class ComponentFractionKernel
 * @tparam NUM_COMP number of fluid components
 * @brief Define the interface for the update kernel in charge of computing the phase volume fractions
 */
template< integer NUM_COMP >
class ComponentFractionKernel : public PropertyKernelBase< NUM_COMP >
{
public:

  using Base = PropertyKernelBase< NUM_COMP >;
  using Base::numComp;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  ComponentFractionKernel( geos::ObjectManagerBase & subRegion )
    : Base(),
    m_compDens( subRegion.getField< geos::fields::flow::globalCompDensity >() ),
    m_compFrac( subRegion.getField< geos::fields::flow::globalCompFraction >() ),
    m_dCompFrac_dCompDens( subRegion.getField< geos::fields::flow::dGlobalCompFraction_dGlobalCompDensity >() )
  { }

  /**
   * @brief Compute the phase volume fractions in an element
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[in] phaseVolFractionKernelOp the function used to customize the kernel
   */
  template< typename FUNC = geos::NoOpFunc >
  GEOS_HOST_DEVICE
  void compute( geos::localIndex const ei,
                FUNC && compFractionKernelOp = geos::NoOpFunc{} ) const
  {
    arraySlice1d< geos::real64 const, geos::compflow::USD_COMP - 1 > const compDens = m_compDens[ei];
    arraySlice1d< real64, compflow::USD_COMP - 1 > const compFrac = m_compFrac[ei];
    arraySlice2d< real64, compflow::USD_COMP_DC - 1 > const dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];

    real64 totalDensity = 0.0;

    for( integer ic = 0; ic < numComp; ++ic )
    {
      totalDensity += compDens[ic];
    }

    real64 const totalDensityInv = 1.0 / totalDensity;

    for( integer ic = 0; ic < numComp; ++ic )
    {
      compFrac[ic] = compDens[ic] * totalDensityInv;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dCompFrac_dCompDens[ic][jc] = -compFrac[ic] * totalDensityInv;
      }
      dCompFrac_dCompDens[ic][ic] += totalDensityInv;
    }

    compFractionKernelOp( compFrac, dCompFrac_dCompDens );
  }

protected:

  // inputs

  // Views on component densities
  arrayView2d< real64 const, compflow::USD_COMP > m_compDens;

  // outputs

  // Views on component fraction
  arrayView2d< real64, compflow::USD_COMP > m_compFrac;
  arrayView3d< real64, compflow::USD_COMP_DC > m_dCompFrac_dCompDens;

};

/**
 * @class ComponentFractionKernelFactory
 */
class ComponentFractionKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComp the number of fluid components
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComp,
                   ObjectManagerBase & subRegion )
  {
    internal::kernelLaunchSelectorCompSwitch( numComp, [&]( auto NC )
    {
      integer constexpr NUM_COMP = NC();
      ComponentFractionKernel< NUM_COMP > kernel( subRegion );
      ComponentFractionKernel< NUM_COMP >::template launch< POLICY >( subRegion.size(), kernel );
    } );
  }

};

}

}

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPONENTFRACTIONKERNEL_HPP
