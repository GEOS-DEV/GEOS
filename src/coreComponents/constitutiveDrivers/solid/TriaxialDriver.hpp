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
 * @file TriaxialDriver.hpp
 */

#ifndef SRC_CORECOMPONENTS_CONSTITUTIVEDRIVERS_SOLID_TRIAXIALDRIVER_HPP
#define SRC_CORECOMPONENTS_CONSTITUTIVEDRIVERS_SOLID_TRIAXIALDRIVER_HPP

#include "constitutiveDrivers/ConstitutiveDriver.hpp"

namespace geos
{

namespace constitutive
{
class SolidBase;
}

/**
 * @class TriaxialDriver
 *
 * Class to allow for triaxial (and similar) tests of the solid constitutive models without the
 * complexity of setting up a single element test.
 *
 */
class TriaxialDriver : public ConstitutiveDriver
{
public:

  enum class Mode
  {
    MixedControl,
    StrainControl,
    StressControl,
  };

  TriaxialDriver( const string & name,
                  Group * const parent );
  ~TriaxialDriver() override = default;

  static string catalogName() { return "TriaxialDriver"; }

  void postInputInitialization() override;

  using ConstitutiveDriver::execute;

  bool execute() override;

  void getColumnNames( string_array & columnNames ) const override;

  /**
   * @brief Run a strain-controlled test using loading protocol in table
   * @param solid Solid constitutive model
   * @param table Table with stress / strain time history
   */
  template< typename SOLID_TYPE >
  void runStrainControlTest( SOLID_TYPE & solid, arrayView2d< real64 > const & table );

  /**
   * @brief Run a stress-controlled test using loading protocol in table
   * @param solid Solid constitutive model
   * @param table Table with stress / strain time history
   */
  template< typename SOLID_TYPE >
  void runStressControlTest( SOLID_TYPE & solid, arrayView2d< real64 > const & table );

  /**
   * @brief Run a mixed stress/strain-controlled test using loading protocol in table
   * @param solid Solid constitutive model
   * @param table Table with stress / strain time history
   */
  template< typename SOLID_TYPE >
  void runMixedControlTest( SOLID_TYPE & solid, arrayView2d< real64 > const & table );

protected:
  /**
   * @brief Validate results by checking residual and removing erroneous data
   */
  bool validateResults() override;

private:
  /**
   * @brief Get the material model from the catalog
   */
  constitutive::SolidBase & getSolid();
  constitutive::SolidBase const & getSolid() const;

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    constexpr static char const * solidMaterialNameString() { return "material"; }
    constexpr static char const * modeString() { return "mode"; }
    constexpr static char const * axialFunctionString() { return "axialControl"; }
    constexpr static char const * radialFunctionString() { return "radialControl"; }
    constexpr static char const * initialStressString() { return "initialStress"; }
  };

  string m_solidMaterialName;  ///< Material identifier
  Mode m_mode;                 ///< Test mode: strainControl, stressControl, mixedControl
  string m_axialFunctionName;  ///< Time-dependent function controlling axial stress or strain (depends on test mode)
  string m_radialFunctionName; ///< Time-dependent function controlling radial stress or strain (depends on test mode)
  real64 m_initialStress;      ///< Initial stress value (scalar used to set an isotropic stress state)

  enum columnKeys { EPS0 = 1, EPS1, EPS2, SIG0, SIG1, SIG2, ITER, NORM }; ///< Enumeration of column keys

  static constexpr integer m_maxIter = 25;   ///< Max Newton iterations for mixed-control tests
  static constexpr integer m_maxCuts = 8;    ///< Max backtracking cuts in line search algorithm
  static constexpr real64 m_newtonTol = 1e-6;   ///< Newton tolerance for mixed-control tests
};

/// Declare strings associated with enumeration values.
ENUM_STRINGS( TriaxialDriver::Mode,
              "mixedControl",
              "strainControl",
              "stressControl" );

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_CONSTITUTIVEDRIVERS_SOLID_TRIAXIALDRIVER_HPP */
