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


// Source includes
#include "codingUtilities/UnitTestUtilities.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/CubicEOSPhaseModel.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geos;
using namespace geos::testing;
using namespace geos::constitutive;

static constexpr real64 relTol = 1e-5;

TEST( CubicEOSTest, testCubicEOSTwoComponentsSRK )
{
  constexpr integer numComps = 2;

  real64 pressure = 0.0;
  real64 temperature = 0.0;
  array1d< real64 > composition( numComps );
  array1d< real64 > criticalPressure( numComps );
  array1d< real64 > criticalTemperature( numComps );
  array1d< real64 > omega( numComps );
  real64 binaryInteractionCoefficients = 0.0; // not implemented yet
  array1d< real64 > logFugacityCoefficients( numComps );
  array1d< real64 > expectedLogFugacityCoefficients( numComps );

  criticalPressure[0] = 12.96e5;
  criticalPressure[1] = 45.99e5;

  criticalTemperature[0] = 33.15;
  criticalTemperature[1] = 190.6;

  omega[0] = -0.219;
  omega[1] = 0.0114;

  ///////////////////////////////////////////

  pressure = 1e6;
  temperature = 350;
  composition[0] = 0.1;
  composition[1] = 0.9;

  expectedLogFugacityCoefficients[0] = 0.0126163;
  expectedLogFugacityCoefficients[1] = -0.00820777;

  CubicEOSPhaseModel< SoaveRedlichKwongEOS >::
  compute( numComps,
           pressure, temperature, composition,
           criticalPressure, criticalTemperature, omega, binaryInteractionCoefficients,
           logFugacityCoefficients );

  checkRelativeError( logFugacityCoefficients[0], expectedLogFugacityCoefficients[0], relTol );
  checkRelativeError( logFugacityCoefficients[1], expectedLogFugacityCoefficients[1], relTol );

  ///////////////////////////////////////////

  pressure = 5e7;
  temperature = 350;
  composition[0] = 0.1;
  composition[1] = 0.9;

  expectedLogFugacityCoefficients[0] = 0.481514;
  expectedLogFugacityCoefficients[1] = -0.0701117;

  CubicEOSPhaseModel< SoaveRedlichKwongEOS >::
  compute( numComps,
           pressure, temperature, composition,
           criticalPressure, criticalTemperature, omega, binaryInteractionCoefficients,
           logFugacityCoefficients );

  checkRelativeError( logFugacityCoefficients[0], expectedLogFugacityCoefficients[0], relTol );
  checkRelativeError( logFugacityCoefficients[1], expectedLogFugacityCoefficients[1], relTol );

  ///////////////////////////////////////////

  pressure = 1e6;
  temperature = 350;
  composition[0] = 0.5;
  composition[1] = 0.5;

  expectedLogFugacityCoefficients[0] = 0.00721367;
  expectedLogFugacityCoefficients[1] = -0.00589892;

  CubicEOSPhaseModel< SoaveRedlichKwongEOS >::
  compute( numComps,
           pressure, temperature, composition,
           criticalPressure, criticalTemperature, omega, binaryInteractionCoefficients,
           logFugacityCoefficients );

  checkRelativeError( logFugacityCoefficients[0], expectedLogFugacityCoefficients[0], relTol );
  checkRelativeError( logFugacityCoefficients[1], expectedLogFugacityCoefficients[1], relTol );

  ///////////////////////////////////////////

  pressure = 5e7;
  temperature = 350;
  composition[0] = 0.5;
  composition[1] = 0.5;

  expectedLogFugacityCoefficients[0] = 0.334027;
  expectedLogFugacityCoefficients[1] = -0.00629384;

  CubicEOSPhaseModel< SoaveRedlichKwongEOS >::
  compute( numComps,
           pressure, temperature, composition,
           criticalPressure, criticalTemperature, omega, binaryInteractionCoefficients,
           logFugacityCoefficients );

  checkRelativeError( logFugacityCoefficients[0], expectedLogFugacityCoefficients[0], relTol );
  checkRelativeError( logFugacityCoefficients[1], expectedLogFugacityCoefficients[1], relTol );

}

TEST( CubicEOSTest, testCubicEOSFourComponentsPR )
{
  constexpr integer numComps = 4;

  real64 pressure = 0.0;
  real64 temperature = 0.0;
  array1d< real64 > composition( numComps );
  array1d< real64 > criticalPressure( numComps );
  array1d< real64 > criticalTemperature( numComps );
  array1d< real64 > omega( numComps );
  real64 binaryInteractionCoefficients = 0.0; // not implemented yet
  array1d< real64 > logFugacityCoefficients( numComps );
  array1d< real64 > expectedLogFugacityCoefficients( numComps );

  criticalPressure[0] = 34e5;
  criticalPressure[1] = 25.3e5;
  criticalPressure[2] = 14.6e5;
  criticalPressure[3] = 220.5e5;

  criticalTemperature[0] = 126.2;
  criticalTemperature[1] = 622.0;
  criticalTemperature[2] = 782.0;
  criticalTemperature[3] = 647.0;

  omega[0] = 0.04;
  omega[1] = 0.443;
  omega[2] = 0.816;
  omega[3] = 0.344;

  ///////////////////////////////////////////

  pressure = 1e7;
  temperature = 297.15;
  composition[0] = 0.0569514;
  composition[1] = 0.104818;
  composition[2] = 0.104822;
  composition[3] = 0.733409;

  expectedLogFugacityCoefficients[0] = 2.8298;
  expectedLogFugacityCoefficients[1] = -8.88628;
  expectedLogFugacityCoefficients[2] = -17.0201;
  expectedLogFugacityCoefficients[3] = -5.33003;

  CubicEOSPhaseModel< PengRobinsonEOS >::
  compute( numComps,
           pressure, temperature, composition,
           criticalPressure, criticalTemperature, omega, binaryInteractionCoefficients,
           logFugacityCoefficients );

  checkRelativeError( logFugacityCoefficients[0], expectedLogFugacityCoefficients[0], relTol );
  checkRelativeError( logFugacityCoefficients[1], expectedLogFugacityCoefficients[1], relTol );
  checkRelativeError( logFugacityCoefficients[2], expectedLogFugacityCoefficients[2], relTol );
  checkRelativeError( logFugacityCoefficients[3], expectedLogFugacityCoefficients[3], relTol );

  ///////////////////////////////////////////

  pressure = 1e5;
  temperature = 297.15;
  composition[0] = 0.00185559;
  composition[1] = 0.332324;
  composition[2] = 0.664862;
  composition[3] = 0.000958244;

  expectedLogFugacityCoefficients[0] = 6.28652;
  expectedLogFugacityCoefficients[1] = -5.83771;
  expectedLogFugacityCoefficients[2] = -16.638;
  expectedLogFugacityCoefficients[3] = 0.361984;

  CubicEOSPhaseModel< PengRobinsonEOS >::
  compute( numComps,
           pressure, temperature, composition,
           criticalPressure, criticalTemperature, omega, binaryInteractionCoefficients,
           logFugacityCoefficients );

  checkRelativeError( logFugacityCoefficients[0], expectedLogFugacityCoefficients[0], relTol );
  checkRelativeError( logFugacityCoefficients[1], expectedLogFugacityCoefficients[1], relTol );
  checkRelativeError( logFugacityCoefficients[2], expectedLogFugacityCoefficients[2], relTol );
  checkRelativeError( logFugacityCoefficients[3], expectedLogFugacityCoefficients[3], relTol );

  ///////////////////////////////////////////

  pressure = 4.78429e+06;
  temperature = 297.15;
  composition[0] = 0.0566196;
  composition[1] = 0.31411;
  composition[2] = 0.628223;
  composition[3] = 0.001047;

  expectedLogFugacityCoefficients[0] = 2.49484;
  expectedLogFugacityCoefficients[1] = -9.36508;
  expectedLogFugacityCoefficients[2] = -19.8123;
  expectedLogFugacityCoefficients[3] = -3.42481;

  CubicEOSPhaseModel< PengRobinsonEOS >::
  compute( numComps,
           pressure, temperature, composition,
           criticalPressure, criticalTemperature, omega, binaryInteractionCoefficients,
           logFugacityCoefficients );

  checkRelativeError( logFugacityCoefficients[0], expectedLogFugacityCoefficients[0], relTol );
  checkRelativeError( logFugacityCoefficients[1], expectedLogFugacityCoefficients[1], relTol );
  checkRelativeError( logFugacityCoefficients[2], expectedLogFugacityCoefficients[2], relTol );
  checkRelativeError( logFugacityCoefficients[3], expectedLogFugacityCoefficients[3], relTol );

}

TEST( CubicEOSTest, testCubicEOSFourComponentsSRK )
{
  integer const numComps = 4;

  real64 pressure = 0.0;
  real64 temperature = 0.0;
  array1d< real64 > composition( 4 );
  array1d< real64 > criticalPressure( 4 );
  array1d< real64 > criticalTemperature( 4 );
  array1d< real64 > omega( 4 );
  real64 binaryInteractionCoefficients = 0.0; // not implemented yet
  array1d< real64 > logFugacityCoefficients( 4 );
  array1d< real64 > expectedLogFugacityCoefficients( 4 );

  criticalPressure[0] = 34e5;
  criticalPressure[1] = 25.3e5;
  criticalPressure[2] = 14.6e5;
  criticalPressure[3] = 220.5e5;

  criticalTemperature[0] = 126.2;
  criticalTemperature[1] = 622.0;
  criticalTemperature[2] = 782.0;
  criticalTemperature[3] = 647.0;

  omega[0] = 0.04;
  omega[1] = 0.443;
  omega[2] = 0.816;
  omega[3] = 0.344;

  ///////////////////////////////////////////

  pressure = 1e7;
  temperature = 297.15;
  composition[0] = 0.994214;
  composition[1] = 6.05198e-05;
  composition[2] = 5.98122e-08;
  composition[3] = 0.00572563;

  expectedLogFugacityCoefficients[0] = 0.00588361;
  expectedLogFugacityCoefficients[1] = -1.44445;
  expectedLogFugacityCoefficients[2] = -2.83086;
  expectedLogFugacityCoefficients[3] = -0.618972;

  CubicEOSPhaseModel< SoaveRedlichKwongEOS >::
  compute( numComps,
           pressure, temperature, composition,
           criticalPressure, criticalTemperature, omega, binaryInteractionCoefficients,
           logFugacityCoefficients );

  checkRelativeError( logFugacityCoefficients[0], expectedLogFugacityCoefficients[0], relTol );
  checkRelativeError( logFugacityCoefficients[1], expectedLogFugacityCoefficients[1], relTol );
  checkRelativeError( logFugacityCoefficients[2], expectedLogFugacityCoefficients[2], relTol );
  checkRelativeError( logFugacityCoefficients[3], expectedLogFugacityCoefficients[3], relTol );

  ///////////////////////////////////////////

  pressure = 1e5;
  temperature = 297.15;
  composition[0] = 0.997965;
  composition[1] = 0.000851981;
  composition[2] = 2.89283e-08;
  composition[3] = 0.00118249;

  expectedLogFugacityCoefficients[0] = -5.94544e-05;
  expectedLogFugacityCoefficients[1] = -0.0168209;
  expectedLogFugacityCoefficients[2] = -0.0334318;
  expectedLogFugacityCoefficients[3] = -0.00664411;

  CubicEOSPhaseModel< SoaveRedlichKwongEOS >::
  compute( numComps,
           pressure, temperature, composition,
           criticalPressure, criticalTemperature, omega, binaryInteractionCoefficients,
           logFugacityCoefficients );

  checkRelativeError( logFugacityCoefficients[0], expectedLogFugacityCoefficients[0], relTol );
  checkRelativeError( logFugacityCoefficients[1], expectedLogFugacityCoefficients[1], relTol );
  checkRelativeError( logFugacityCoefficients[2], expectedLogFugacityCoefficients[2], relTol );
  checkRelativeError( logFugacityCoefficients[3], expectedLogFugacityCoefficients[3], relTol );

  ///////////////////////////////////////////

  pressure = 1.83959e+06;
  temperature = 297.15;
  composition[0] = 0.0309329;
  composition[1] = 0.319683;
  composition[2] = 0.637861;
  composition[3] = 0.011523;

  expectedLogFugacityCoefficients[0] = 3.47428;
  expectedLogFugacityCoefficients[1] = -8.75355;
  expectedLogFugacityCoefficients[2] = -19.6075;
  expectedLogFugacityCoefficients[3] = -2.69792;

  CubicEOSPhaseModel< SoaveRedlichKwongEOS >::
  compute( numComps,
           pressure, temperature, composition,
           criticalPressure, criticalTemperature, omega, binaryInteractionCoefficients,
           logFugacityCoefficients );

  checkRelativeError( logFugacityCoefficients[0], expectedLogFugacityCoefficients[0], relTol );
  checkRelativeError( logFugacityCoefficients[1], expectedLogFugacityCoefficients[1], relTol );
  checkRelativeError( logFugacityCoefficients[2], expectedLogFugacityCoefficients[2], relTol );
  checkRelativeError( logFugacityCoefficients[3], expectedLogFugacityCoefficients[3], relTol );

}
