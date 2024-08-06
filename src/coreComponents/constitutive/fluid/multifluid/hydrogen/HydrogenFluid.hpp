/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HydrogenFluid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_HYDROGEN_HYDROGENFLUID_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_HYDROGEN_HYDROGENFLUID_HPP_

#include "HydrogenFluidUpdate.hpp"

#include <memory>

namespace geos
{

namespace constitutive
{

class HydrogenFluid : public MultiFluidBase
{
public:
  using exec_policy = parallelDevicePolicy<>;

  HydrogenFluid( string const & name,
                 Group * const parent );

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                Group * const parent ) const override;

  static string catalogName();

  virtual string getCatalogName() const override { return catalogName(); }

  static constexpr bool isThermalType()
  {
    return false;
  }

  static constexpr integer min_n_components = 3;
  static constexpr integer max_n_components = 4;

  bool isThermal() const override
  {
    return isThermalType();
  }

  using KernelWrapper = HydrogenFluidUpdate;

  integer getWaterPhaseIndex() const override final;
  void checkTablesParameters( real64 pressure, real64 temperature ) const override;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

protected:
  void postInputInitialization() override;

private:
  integer m_gasPhaseIndex{-1};
  integer m_watPhaseIndex{-1};
  integer m_h2ComponentIndex{-1};
  integer m_h2oComponentIndex{-1};

};

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_HYDROGEN_HYDROGENFLUID_HPP_
