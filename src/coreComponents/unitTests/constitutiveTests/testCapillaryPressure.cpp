/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source inclues
#include "constitutiveTestHelpers.hpp"
#include "mainInterface/initialization.hpp"

using namespace geosx;
using namespace geosx::testing;
using namespace geosx::constitutive;
using namespace geosx::constitutive::cappres;
using namespace geosx::dataRepository;

CapillaryPressureBase & makeBrooksCoreyCapPressureTwoPhase( string const & name, Group & parent )
{
  BrooksCoreyCapillaryPressure & capPressure = parent.registerGroup< BrooksCoreyCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "water"; phaseNames[1] = "oil";

  array1d< real64 > & phaseMinSat = capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString() );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.1; phaseMinSat[1] = 0.05;

  array1d< real64 > & phaseCapPressureExpInv =
    capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString() );
  phaseCapPressureExpInv.resize( 2 );
  phaseCapPressureExpInv[0] = 2; phaseCapPressureExpInv[1] = 4;

  array1d< real64 > & phaseEntryPressure = capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseEntryPressureString() );
  phaseEntryPressure.resize( 2 );
  phaseEntryPressure[0] = 1; phaseEntryPressure[1] = 1;

  real64 & capPressureEpsilon = capPressure.getReference< real64 >( BrooksCoreyCapillaryPressure::viewKeyStruct::capPressureEpsilonString() );
  capPressureEpsilon = 1e-4;

  capPressure.postProcessInputRecursive();
  return capPressure;
}


CapillaryPressureBase & makeBrooksCoreyCapPressureThreePhase( string const & name, Group & parent )
{
  BrooksCoreyCapillaryPressure & capPressure = parent.registerGroup< BrooksCoreyCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 3 );
  phaseNames[0] = "water"; phaseNames[1] = "oil"; phaseNames[2] = "gas";

  array1d< real64 > & phaseMinSat = capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString() );
  phaseMinSat.resize( 3 );
  phaseMinSat[0] = 0.04; phaseMinSat[1] = 0.02; phaseMinSat[2] = 0.1;

  array1d< real64 > & phaseCapPressureExpInv =
    capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString() );
  phaseCapPressureExpInv.resize( 3 );
  phaseCapPressureExpInv[0] = 2; phaseCapPressureExpInv[1] = -3; phaseCapPressureExpInv[2] = 2.5;

  array1d< real64 > & phaseEntryPressure = capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseEntryPressureString() );
  phaseEntryPressure.resize( 3 );
  phaseEntryPressure[0] = 1; phaseEntryPressure[1] = -1; phaseEntryPressure[2] = 2;

  real64 & capPressureEpsilon = capPressure.getReference< real64 >( BrooksCoreyCapillaryPressure::viewKeyStruct::capPressureEpsilonString() );
  capPressureEpsilon = 1e-7;

  capPressure.postProcessInputRecursive();
  return capPressure;
}


CapillaryPressureBase & makeVanGenuchtenCapPressureTwoPhase( string const & name, Group & parent )
{
  VanGenuchtenCapillaryPressure & capPressure = parent.registerGroup< VanGenuchtenCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas";

  array1d< real64 > & phaseMinSat = capPressure.getReference< array1d< real64 > >( VanGenuchtenCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString() );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.04; phaseMinSat[1] = 0.1;

  array1d< real64 > & phaseCapPressureExpInv = capPressure.getReference< array1d< real64 > >(
    VanGenuchtenCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString() );
  phaseCapPressureExpInv.resize( 2 );
  phaseCapPressureExpInv[0] = 0.4; phaseCapPressureExpInv[1] = 0.5;

  array1d< real64 > & phaseCapPressureMultiplier = capPressure.getReference< array1d< real64 > >(
    VanGenuchtenCapillaryPressure::viewKeyStruct::phaseCapPressureMultiplierString() );
  phaseCapPressureMultiplier.resize( 2 );
  phaseCapPressureMultiplier[0] = 0.5; phaseCapPressureMultiplier[1] = 1;

  real64 & capPressureEpsilon = capPressure.getReference< real64 >( VanGenuchtenCapillaryPressure::viewKeyStruct::capPressureEpsilonString() );
  capPressureEpsilon = 1e-4;

  capPressure.postProcessInputRecursive();
  return capPressure;
}

CapillaryPressureBase & makeVanGenuchtenCapPressureThreePhase( string const & name, Group & parent )
{
  VanGenuchtenCapillaryPressure & capPressure = parent.registerGroup< VanGenuchtenCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas"; phaseNames[2] = "water";

  array1d< real64 > & phaseMinSat = capPressure.getReference< array1d< real64 > >( VanGenuchtenCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString() );
  phaseMinSat.resize( 3 );
  phaseMinSat[0] = 0.04; phaseMinSat[1] = 0.1; phaseMinSat[2] = 0.;

  array1d< real64 > & phaseCapPressureExpInv = capPressure.getReference< array1d< real64 > >(
    VanGenuchtenCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString() );
  phaseCapPressureExpInv.resize( 3 );
  phaseCapPressureExpInv[0] = 0.33; phaseCapPressureExpInv[1] = 0.4; phaseCapPressureExpInv[ 2 ] = 0.5;

  array1d< real64 > & phaseCapPressureMultiplier = capPressure.getReference< array1d< real64 > >(
    VanGenuchtenCapillaryPressure::viewKeyStruct::phaseCapPressureMultiplierString() );
  phaseCapPressureMultiplier.resize( 3 );
  phaseCapPressureMultiplier[0] = 0.5; phaseCapPressureMultiplier[1] = 1; phaseCapPressureMultiplier[2] = 0.2;

  real64 & capPressureEpsilon = capPressure.getReference< real64 >( VanGenuchtenCapillaryPressure::viewKeyStruct::capPressureEpsilonString() );
  capPressureEpsilon = 1e-4;

  capPressure.postProcessInputRecursive();
  return capPressure;
}


class CapillaryPressureTest : public ConstitutiveTestBase< CapillaryPressureBase >
{
public:
  void test( arraySlice1d< real64 const > const sat, real64 const eps, real64 const tol )
  {
    arrayView3d< real64 const, USD_CAPPRES > phaseCapPressure;
    arrayView4d< real64 const, USD_CAPPRES_DS > dPhaseCapPressure_dPhaseVolFraction;
    testNumericalDerivatives( m_parent,
                              *m_model,
                              sat,
                              eps,
                              tol,
                              "phaseRelPerm",
                              [&phaseCapPressure] ( CapillaryPressureBase & relPerm )
    {
      phaseCapPressure = relPerm.phaseCapPressure();
      return phaseCapPressure[ 0 ][ 0 ];
    },
                              [&dPhaseCapPressure_dPhaseVolFraction] ( CapillaryPressureBase & relPerm )
    {
      dPhaseCapPressure_dPhaseVolFraction = relPerm.dPhaseCapPressure_dPhaseVolFraction();
      return dPhaseCapPressure_dPhaseVolFraction[ 0 ][ 0 ];
    }
                              );
  }
};

TEST_F( CapillaryPressureTest, numericalDerivatives_brooksCoreyCapPressureTwoPhase )
{
  initialize( makeBrooksCoreyCapPressureTwoPhase( "capPressure", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.4;
  real64 const end_sat   = 0.6;
  real64 const dS = 1e-1;
  array1d< real64 > sat( 2 );
  sat[0] = start_sat; sat[1] = 1.0-sat[0];
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1 - sat[0];
  }
}


TEST_F( CapillaryPressureTest, numericalDerivatives_brooksCoreyCapPressureThreePhase )
{
  initialize( makeBrooksCoreyCapPressureThreePhase( "capPressure", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.4;
  real64 const end_sat   = 0.6;
  real64 const dS = 1e-1;
  array1d< real64 > sat( 3 );
  sat[0] = start_sat;
  sat[1] = 0.5*(1-sat[0]);
  sat[2] = 1.0-sat[0]-sat[1];
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 0.5 * ( 1-sat[0] );
    sat[2] = 1.0 - sat[0] - sat[1];
  }
}


TEST_F( CapillaryPressureTest, numericalDerivatives_vanGenuchtenCapPressureTwoPhase )
{
  initialize( makeVanGenuchtenCapPressureTwoPhase( "capPressure", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.4;
  real64 const end_sat   = 0.6;
  real64 const dS        = 1e-1;
  array1d< real64 > sat( 2 );
  sat[0] = start_sat; sat[1] = 1-sat[1];
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1 - sat[0];
  }

}


TEST_F( CapillaryPressureTest, numericalDerivatives_vanGenuchtenCapPressureThreePhase )
{
  initialize( makeVanGenuchtenCapPressureThreePhase( "capPressure", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.4;
  real64 const end_sat   = 0.6;
  real64 const dS        = 1e-1;
  array1d< real64 > sat( 3 );
  sat[0] = start_sat;
  sat[1] = 0.5*(1-sat[0]);
  sat[2] = 1.0-sat[0]-sat[1];
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 0.5*(1-sat[0]);
    sat[2] = 1 - sat[0] - sat[1];
  }
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
