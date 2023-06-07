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
 * @file PropertyKernelBase.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_PROPERTYKERNELBASE_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_PROPERTYKERNELBASE_HPP


namespace geos
{

namespace isothermalCompositionalMultiphaseBaseKernels
{

/******************************** PropertyKernelBase ********************************/

/**
 * @class PropertyKernelBase
 * @tparam NUM_COMP number of fluid components
 * @brief Define the base interface for the property update kernels
 */
template< integer NUM_COMP >
class PropertyKernelBase
{
public:

  /// Compile time value for the number of components
  static constexpr geos::integer numComp = NUM_COMP;

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[inout] kernelComponent the kernel component providing access to the compute function
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( geos::localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    geos::forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( geos::localIndex const ei )
    {
      kernelComponent.compute( ei );
    } );
  }

  /**
   * @brief Performs the kernel launch on a sorted array
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] targetSet the indices of the elements in which we compute the property
   * @param[inout] kernelComponent the kernel component providing access to the compute function
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( geos::SortedArrayView< geos::localIndex const > const & targetSet,
          KERNEL_TYPE const & kernelComponent )
  {
    geos::forAll< POLICY >( targetSet.size(), [=] GEOS_HOST_DEVICE ( geos::localIndex const i )
    {
      geos::localIndex const ei = targetSet[i];
      kernelComponent.compute( ei );
    } );
  }

};

}

}

#endif //GEOSX_PROPERTYKERNELBASE_HPP
