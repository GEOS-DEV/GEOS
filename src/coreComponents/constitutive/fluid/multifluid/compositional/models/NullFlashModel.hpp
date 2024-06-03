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
 * @file NullFlashModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NULLFLASHMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NULLFLASHMODEL_HPP_

#include "NullModel.hpp"
#include "ModelParameters.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class NullFlashModelUpdate final : public FunctionBaseUpdate
{
public:
  NullFlashModelUpdate() = default;

  // Mark as a 2-phase flash
  GEOS_HOST_DEVICE
  static constexpr integer getNumberOfPhases() { return 2; }
};

class NullFlashModel : public FunctionBase
{
public:
  NullFlashModel( string const & name,
                              ComponentProperties const & componentProperties )
                              :FunctionBase(name, componentProperties)
                              {}

  FunctionType functionType() const override
  {
    return FunctionType::FLASH;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = NullFlashModelUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;
};

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NULLFLASHMODEL_HPP_
