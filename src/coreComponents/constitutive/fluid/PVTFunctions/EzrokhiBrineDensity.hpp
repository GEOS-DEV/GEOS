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
 * @file EzrokhiBrineViscosity.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_EZROKHIBRINEDENSITY_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_EZROKHIBRINEDENSITY_HPP_

#include "EzrokhiBrineViscosity.hpp"

#include "PVTFunctionBase.hpp"

namespace geosx
{

namespace constitutive
{

namespace PVTProps
{

class EzrokhiBrineDensity : public PVTFunctionBase
{
public:

  EzrokhiBrineDensity( string const & name,
                       string_array const & inputPara,
                       string_array const & componentNames,
                       array1d< real64 > const & componentMolarWeight );

  virtual ~EzrokhiBrineDensity() override = default;

  static string catalogName() { return "EzrokhiBrineDensity"; }

  virtual string getCatalogName() const override final { return catalogName(); }

  virtual PVTFunctionType functionType() const override
  {
    return PVTFunctionType::DENSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = EzrokhiBrineViscosityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

private:

  void makeCoefficients( string_array const & inputPara );

  /// Index of the CO2 phase
  integer m_CO2Index;

  /// Index of the water phase
  integer m_waterIndex;

  real64 m_coef0;

  real64 m_coef1;

  real64 m_coef2;

  real64 m_coef3;

};


} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_EZROKHIBRINEVISCOSITY_HPP_
