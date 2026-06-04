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

#ifndef GEOS_CONSTITUTIVEDRIVERS_CONTACT_FRICTIONDRIVER_HPP
#define GEOS_CONSTITUTIVEDRIVERS_CONTACT_FRICTIONDRIVER_HPP

#include "constitutiveDrivers/ConstitutiveDriver.hpp"
#include "physicsSolvers/solidMechanics/contact/ContactSolverBase.hpp"
// #include "physicsSolvers/solidMechanics/contact/SolidMechanicsAugmentedLagrangianContact.hpp"

namespace geos
{
namespace constitutive
{
class FrictionBase;
}

class FrictionDriver : public ConstitutiveDriver
{
public:
  FrictionDriver( const string & name,
                  Group * const parent );

  static string catalogName()
  { return "FrictionDriver"; }

  void postInputInitialization() override;

  bool execute() override;

  void getColumnNames( string_array & columnNames ) const override;

  template< typename FRICTION_TYPE >
  void
  runTest( FRICTION_TYPE & friction,
           const arrayView2d< real64, 1 > & table );

private:
  /**
   * @brief Get the friction model from the catalog
   */
  constitutive::FrictionBase & getFriction();
  constitutive::FrictionBase const & getFriction() const;

  void initializeTable();

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct : ConstitutiveDriver::viewKeyStruct
  {
    constexpr static char const * frictionNameString()
    { return "friction"; }

    constexpr static char const * contactNameString()
    { return "contact"; }

    constexpr static char const * jumpFunctionString()
    { return "jumpControl"; }

    constexpr static char const * dJumpFunctionString()
    { return "dJumpControl"; }

    constexpr static char const * tractionFunctionString()
    { return "tractionControl"; }

    constexpr static char const * thetaString()
    { return "xTiltAngle";}

    constexpr static char const * phiString()
    { return "yTiltAngle";}
  };

  // Time is defined in base class
  enum columnKeys { NJUMP = 1, SLIP0, SLIP1, NDJUMP, DSLIP0, DSLIP1, NTRAC, STRAC0, STRAC1, FS, NEWTRAC, SNEWTRAC0, SNEWTRAC1, TLIM };

  string m_jumpFunctionName; ///<
  string m_dJumpFunctionName; ///<
  string m_tractionFunctionName; ///<

  real64 m_theta{0.0}; ///< x-tilt of fault
  real64 m_phi{0.0};  ///< y-tilt of fault

  string m_frictionName;               ///< frictionType identifier
};

}

#endif //GEOS_CONSTITUTIVEDRIVERS_CONTACT_FRICTIONDRIVER_HPP
