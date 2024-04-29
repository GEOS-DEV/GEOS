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
 * @file RachfordRice.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_RACHFORDRICE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_RACHFORDRICE_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

namespace constitutive
{

struct RachfordRice
{
public:

  /// Tolerance of the SSI loop
  static constexpr real64 SSITolerance = 1e-3;
  /// Tolerance of the Newton loop
  static constexpr real64 newtonTolerance = 1e-12;
  /// Max number of SSI iterations
  static constexpr integer maxSSIIterations = 200;
  /// Max number of Newton iterations
  static constexpr integer maxNewtonIterations = 30;
  /// Epsilon used in the calculations
  static constexpr real64 epsilon = LvArray::NumericLimits< real64 >::epsilon;

  /**
   * @brief Function solving the Rachford-Rice equation
   * @input[in] kValues the array fo K-values
   * @input[in] feed the component fractions
   * @input[in] presentComponentIds the indices of components with a non-zero fractions
   * @return the gas mole fraction
   **/
  GEOS_HOST_DEVICE
  static real64
  solve( arrayView1d< real64 const > const kValues,
         arrayView1d< real64 const > const feed,
         arrayView1d< integer const > const presentComponentIds );

private:

  /**
   * @brief Function evaluating the Rachford-Rice function
   * @input[in] kValues the array fo K-values
   * @input[in] feed the component fractions
   * @input[in] presentComponentIds the indices of components with a non-zero fractions
   * @input[in] x the value at which the Rachford-Rice function is evaluated
   * @return the value of the Rachford-Rice function at x
   **/
  GEOS_HOST_DEVICE
  static real64
  evaluate( arrayView1d< real64 const > const kValues,
            arrayView1d< real64 const > const feed,
            arrayView1d< integer const > const presentComponentIds,
            real64 const & x );

  /**
   * @brief Function evaluating the derivative of the Rachford-Rice function
   * @input[in] kValues the array fo K-values
   * @input[in] feed the component fractions
   * @input[in] presentComponentIds the indices of components with a non-zero fractions
   * @input[in] x the value at which the derivative of the Rachford-Rice function is evaluated
   * @return the value of the derivative of the Rachford-Rice function at x
   **/
  GEOS_HOST_DEVICE
  static real64
  evaluateDerivative( arrayView1d< real64 const > const kValues,
                      arrayView1d< real64 const > const feed,
                      arrayView1d< integer const > const presentComponentIds,
                      real64 const & x );

};

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_RACHFORDRICE_HPP_
