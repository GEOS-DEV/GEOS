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
 * @file FluidUpdateKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_KERNELS_COMPOSITIONAL_FLUIDUPDATEKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_KERNELS_COMPOSITIONAL_FLUIDUPDATEKERNEL_HPP

#include "common/DataTypes.hpp"

namespace geos
{
namespace constitutive
{
class MultiFluidBase;
}
namespace thermalCompositionalMultiphaseBaseKernels
{

/******************************** FluidUpdateKernel ********************************/

struct FluidUpdate
{
  /**
   * @brief Update fluid properties
   * @details Will update the properties on the fluid model based on the provided primary variables.
   *          This will use the execution policy that is specified by the fluid model itself.
   * @param[in] size - The number of points to update
   * @param[inout] fluid - The fluid model object with the properties
   * @param[in] pressure - The input pressure values
   * @param[in] temperature - The input temperature values
   * @param[in] componentFractions - The input total composition values
   */
  static void update( localIndex const size,
                      constitutive::MultiFluidBase & fluid,
                      arrayView1d< real64 const > const & pressure,
                      arrayView1d< real64 const > const & temperature,
                      arrayView2d< real64 const, compflow::USD_COMP > const & componentFractions );

  /**
   * @brief Update fluid properties
   * @details Will update the properties on the fluid model based on the provided primary variables.
   *          This will use the execution policy specified regardless of the policy by the fluid model.
   * @tparam POLICY - The execution policy to use
   * @param[in] size - The number of points to update
   * @param[inout] fluid - The fluid model object with the properties
   * @param[in] pressure - The input pressure values
   * @param[in] temperature - The input temperature values
   * @param[in] componentFractions - The input total composition values
   */
  template< typename POLICY >
  static void localUpdate( localIndex const size,
                           constitutive::MultiFluidBase & fluid,
                           arrayView1d< real64 const > const & pressure,
                           arrayView1d< real64 const > const & temperature,
                           arrayView2d< real64 const, compflow::USD_COMP > const & componentFractions );

  /**
   * @brief Update fluid properties
   * @details Will update the properties on the fluid model based on the provided primary variables.
   *          This will use the execution policy that is specified by the fluid model itself.
   * @param[in] targetSet - The list of points (cells) to update
   * @param[inout] fluid - The fluid model object with the properties
   * @param[in] pressure - The input pressure values
   * @param[in] temperature - The input temperature values
   * @param[in] componentFractions - The input total composition values
   */
  static void update( SortedArrayView< localIndex const > const & targetSet,
                      constitutive::MultiFluidBase & fluid,
                      arrayView1d< real64 const > const & pressure,
                      arrayView1d< real64 const > const & temperature,
                      arrayView2d< real64 const, compflow::USD_COMP > const & componentFractions );
};

template< typename POLICY, typename FLUID_WRAPPER >
struct FluidUpdateKernel
{
  static void
  launch( localIndex const size,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pressure,
          arrayView1d< real64 const > const & temperature,
          arrayView2d< real64 const, compflow::USD_COMP > const & componentFractions );

  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pressure,
          arrayView1d< real64 const > const & temperature,
          arrayView2d< real64 const, compflow::USD_COMP > const & componentFractions );
};

} // namespace thermalCompositionalMultiphaseBaseKernels
} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_KERNELS_COMPOSITIONAL_FLUIDUPDATEKERNEL_HPP
