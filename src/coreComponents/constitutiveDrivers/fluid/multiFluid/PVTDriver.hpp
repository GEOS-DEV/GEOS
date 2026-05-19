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
 * @file PVTDriver.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_PVTDRIVER_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_PVTDRIVER_HPP_

#include "constitutiveDrivers/ConstitutiveDriver.hpp"

namespace geos
{

namespace constitutive
{
class MultiFluidBase;
}

/**
 * @class PVTDriver
 *
 * Class to allow for testing PVT behavior without the
 * complexity of setting up a full simulation.
 *
 */
class PVTDriver : public ConstitutiveDriver
{
public:
  PVTDriver( const string & name,
             Group * const parent );

  static string catalogName() { return "PVTDriver"; }

  void postInputInitialization() override;

  using ConstitutiveDriver::execute;

  bool execute() override;

  void getColumnNames( string_array & columnNames ) const override;

  /**
   * @brief Run test using loading protocol in table
   * @param i Fluid constitutive model
   * @param table Table with input/output time history
   */
  template< typename FLUID_TYPE >
  void runTest( FLUID_TYPE & fluid, arrayView2d< real64 > const & table );

private:
  /**
   * @brief Get the fluid model from the catalog
   */
  constitutive::MultiFluidBase & getFluid();
  constitutive::MultiFluidBase const & getFluid() const;

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct : ConstitutiveDriver::viewKeyStruct
  {
    constexpr static char const * fluidNameString() { return "fluid"; }
    constexpr static char const * pressureFunctionString() { return "pressureControl"; }
    constexpr static char const * temperatureFunctionString() { return "temperatureControl"; }
    constexpr static char const * feedString() { return "feedComposition"; }
    constexpr static char const * outputMassDensityString() { return "outputMassDensity"; }
    constexpr static char const * outputCompressibilityString() { return "outputCompressibility"; }
    constexpr static char const * outputPhaseCompositionString() { return "outputPhaseComposition"; }
  };

  integer m_numPhases;     ///< Number of fluid phases
  integer m_numComponents; ///< Number of fluid components

  string m_fluidName;                   ///< Fluid identifier
  string m_pressureFunctionName;        ///< Time-dependent function controlling pressure
  string m_temperatureFunctionName;     ///< Time-dependent function controlling temperature
  integer m_outputMassDensity{0};       ///< Flag to indicate that the mass density of each phase should be output
  integer m_outputCompressibility{0};   ///< Flag to indicate that the total compressibility should be output
  integer m_outputPhaseComposition{0};  ///< Flag to indicate that phase compositions should be output
  array1d< real64 > m_feed;             ///< User specified feed composition

  enum columnKeys { PRES = 1, TEMP }; ///< Enumeration of "input" column keys for readability
};

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_FLUID_PVTDRIVER_HPP_ */
