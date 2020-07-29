/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "constitutiveTestHelpers.hpp"
#include "managers/initialization.hpp"

using namespace geosx;
using namespace geosx::testing;
using namespace geosx::constitutive;
using namespace geosx::dataRepository;

RelativePermeabilityBase * makeBrooksCoreyRelPerm( string const & name, Group & parent )
{
  auto relPerm = parent.RegisterGroup< BrooksCoreyRelativePermeability >( name );

  auto & phaseNames = relPerm->getReference< string_array >( RelativePermeabilityBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 2 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas";

  auto & phaseMinSat = relPerm->getReference< array1d< real64 > >( BrooksCoreyRelativePermeability::viewKeyStruct::phaseMinVolumeFractionString );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.1; phaseMinSat[1] = 0.15;

  auto & phaseRelPermExp = relPerm->getReference< array1d< real64 > >( BrooksCoreyRelativePermeability::viewKeyStruct::phaseRelPermExponentString );
  phaseRelPermExp.resize( 2 );
  phaseRelPermExp[0] = 2.0; phaseRelPermExp[1] = 2.0;

  auto & phaseRelPermMaxVal = relPerm->getReference< array1d< real64 > >( BrooksCoreyRelativePermeability::viewKeyStruct::phaseRelPermMaxValueString );
  phaseRelPermMaxVal.resize( 2 );
  phaseRelPermMaxVal[0] = 0.8; phaseRelPermMaxVal[1] = 0.9;

  relPerm->PostProcessInputRecursive();
  return relPerm;
}

RelativePermeabilityBase * makeBrooksCoreyBakerRelPermTwoPhase( string const & name, Group & parent )
{
  auto relPerm = parent.RegisterGroup< BrooksCoreyBakerRelativePermeability >( name );

  auto & phaseNames = relPerm->getReference< string_array >( RelativePermeabilityBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 2 );
  phaseNames[0] = "water"; phaseNames[1] = "oil";

  auto & phaseMinSat = relPerm->getReference< array1d< real64 > >( BrooksCoreyBakerRelativePermeability::viewKeyStruct::phaseMinVolumeFractionString );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.03; phaseMinSat[1] = 0.01;

  auto & waterOilRelPermExp = relPerm->getReference< array1d< real64 > >( BrooksCoreyBakerRelativePermeability::viewKeyStruct::waterOilRelPermExponentString );
  waterOilRelPermExp.resize( 2 );
  waterOilRelPermExp[0] = 1.9; waterOilRelPermExp[1] = 3.95;

  auto & waterOilRelPermMaxVal =
    relPerm->getReference< array1d< real64 > >( BrooksCoreyBakerRelativePermeability::viewKeyStruct::waterOilRelPermMaxValueString );
  waterOilRelPermMaxVal.resize( 2 );
  waterOilRelPermMaxVal[0] = 0.8; waterOilRelPermMaxVal[1] = 0.75;

  relPerm->PostProcessInputRecursive();
  return relPerm;
}


RelativePermeabilityBase * makeBrooksCoreyBakerRelPermThreePhase( string const & name, Group & parent )
{
  auto relPerm = parent.RegisterGroup< BrooksCoreyBakerRelativePermeability >( name );

  auto & phaseNames = relPerm->getReference< string_array >( RelativePermeabilityBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas"; phaseNames[2] = "water";

  auto & phaseMinSat = relPerm->getReference< array1d< real64 > >( BrooksCoreyBakerRelativePermeability::viewKeyStruct::phaseMinVolumeFractionString );
  phaseMinSat.resize( 3 );
  phaseMinSat[0] = 0.03; phaseMinSat[1] = 0.01; phaseMinSat[2] = 0.025;

  auto & waterOilRelPermExp = relPerm->getReference< array1d< real64 > >( BrooksCoreyBakerRelativePermeability::viewKeyStruct::waterOilRelPermExponentString );
  waterOilRelPermExp.resize( 2 );
  waterOilRelPermExp[0] = 2.4; waterOilRelPermExp[1] = 1.5;

  auto & waterOilRelPermMaxVal =
    relPerm->getReference< array1d< real64 > >( BrooksCoreyBakerRelativePermeability::viewKeyStruct::waterOilRelPermMaxValueString );
  waterOilRelPermMaxVal.resize( 2 );
  waterOilRelPermMaxVal[0] = 0.9; waterOilRelPermMaxVal[1] = 0.65;

  auto & gasOilRelPermExp = relPerm->getReference< array1d< real64 > >( BrooksCoreyBakerRelativePermeability::viewKeyStruct::gasOilRelPermExponentString );
  gasOilRelPermExp.resize( 2 );
  gasOilRelPermExp[0] = 1.9; gasOilRelPermExp[1] = 3.95;

  auto & gasOilRelPermMaxVal = relPerm->getReference< array1d< real64 > >( BrooksCoreyBakerRelativePermeability::viewKeyStruct::gasOilRelPermMaxValueString );
  gasOilRelPermMaxVal.resize( 2 );
  gasOilRelPermMaxVal[0] = 0.8; gasOilRelPermMaxVal[1] = 0.95;

  relPerm->PostProcessInputRecursive();
  return relPerm;
}

RelativePermeabilityBase * makeVanGenuchtenBakerRelPermTwoPhase( string const & name, Group & parent )
{
  auto relPerm = parent.RegisterGroup< VanGenuchtenBakerRelativePermeability >( name );

  auto & phaseNames = relPerm->getReference< string_array >( RelativePermeabilityBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 2 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas";

  auto & phaseMinSat = relPerm->getReference< array1d< real64 > >( VanGenuchtenBakerRelativePermeability::viewKeyStruct::phaseMinVolumeFractionString );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.02; phaseMinSat[1] = 0.05;

  auto & gasOilRelPermExpInv =
    relPerm->getReference< array1d< real64 > >( VanGenuchtenBakerRelativePermeability::viewKeyStruct::gasOilRelPermExponentInvString );
  gasOilRelPermExpInv.resize( 2 );
  gasOilRelPermExpInv[0] = 1.7; gasOilRelPermExpInv[1] = 2.15;

  auto & gasOilRelPermMaxVal = relPerm->getReference< array1d< real64 > >( VanGenuchtenBakerRelativePermeability::viewKeyStruct::gasOilRelPermMaxValueString );
  gasOilRelPermMaxVal.resize( 2 );
  gasOilRelPermMaxVal[0] = 0.5; gasOilRelPermMaxVal[1] = 0.75;

  relPerm->PostProcessInputRecursive();
  return relPerm;
}


RelativePermeabilityBase * makeVanGenuchtenBakerRelPermThreePhase( string const & name, Group & parent )
{
  auto relPerm = parent.RegisterGroup< VanGenuchtenBakerRelativePermeability >( name );

  auto & phaseNames = relPerm->getReference< string_array >( RelativePermeabilityBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas"; phaseNames[2] = "water";

  auto & phaseMinSat = relPerm->getReference< array1d< real64 > >( VanGenuchtenBakerRelativePermeability::viewKeyStruct::phaseMinVolumeFractionString );
  phaseMinSat.resize( 3 );
  phaseMinSat[0] = 0.03; phaseMinSat[1] = 0.01; phaseMinSat[2] = 0.025;

  auto & waterOilRelPermExpInv = relPerm->getReference< array1d< real64 > >(
    VanGenuchtenBakerRelativePermeability::viewKeyStruct::waterOilRelPermExponentInvString );
  waterOilRelPermExpInv.resize( 2 );
  waterOilRelPermExpInv[0] = 2.4; waterOilRelPermExpInv[1] = 2.5;

  auto & waterOilRelPermMaxVal =
    relPerm->getReference< array1d< real64 > >( VanGenuchtenBakerRelativePermeability::viewKeyStruct::waterOilRelPermMaxValueString );
  waterOilRelPermMaxVal.resize( 2 );
  waterOilRelPermMaxVal[0] = 0.9; waterOilRelPermMaxVal[1] = 0.75;

  auto & gasOilRelPermExpInv =
    relPerm->getReference< array1d< real64 > >( VanGenuchtenBakerRelativePermeability::viewKeyStruct::gasOilRelPermExponentInvString );
  gasOilRelPermExpInv.resize( 2 );
  gasOilRelPermExpInv[0] = 1.9; gasOilRelPermExpInv[1] = 3.95;

  auto & gasOilRelPermMaxVal = relPerm->getReference< array1d< real64 > >( BrooksCoreyBakerRelativePermeability::viewKeyStruct::gasOilRelPermMaxValueString );
  gasOilRelPermMaxVal.resize( 2 );
  gasOilRelPermMaxVal[0] = 0.8; gasOilRelPermMaxVal[1] = 0.75;

  relPerm->PostProcessInputRecursive();
  return relPerm;
}

class RelPermTest : public ConstitutiveTestBase< RelativePermeabilityBase >
{
public:
  void test( arraySlice1d< real64 const > const sat, real64 const eps, real64 const tol )
  {
    testNumericalDerivatives( m_parent,
                              *m_model,
                              sat,
                              eps,
                              tol,
                              "phaseRelPerm",
                              [] ( RelativePermeabilityBase & relPerm )
                              { return relPerm.phaseRelPerm()[ 0 ][ 0 ]; },
                              [] ( RelativePermeabilityBase & relPerm )
                              { return relPerm.dPhaseRelPerm_dPhaseVolFraction()[ 0 ][ 0 ]; }
                             );
  }
};

TEST_F( RelPermTest, numericalDerivatives_brooksCoreyRelPerm )
{
  initialize( makeBrooksCoreyRelPerm( "relPerm", m_parent ) );

  // TODO test over a range of values
  array1d< real64 > sat( 4 );
  sat[0] = 0.7; sat[1] = 0.3;

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  test( sat, eps, tol );
}

TEST_F( RelPermTest, numericalDerivatives_BrooksCoreyBakerRelPermTwoPhase )
{
  initialize( makeBrooksCoreyBakerRelPermTwoPhase( "relPerm", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  // TODO test over a range of values
  real64 const start_sat = 0.3;
  real64 const end_sat   = 0.7;
  real64 const dS = 1e-1;
  real64 const alpha = 0.4;
  array1d< real64 > sat( 2 );
  sat[0] = start_sat;
  sat[1] = alpha*(1.0-sat[0]);
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1-sat[0];
  }
}

TEST_F( RelPermTest, numericalDerivatives_BrooksCoreyBakerRelPermThreePhase )
{
  initialize( makeBrooksCoreyBakerRelPermThreePhase( "relPerm", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.3;
  real64 const end_sat   = 0.7;
  real64 const dS = 1e-1;
  real64 const alpha = 0.4;
  array1d< real64 > sat( 3 );
  sat[0] = start_sat;
  sat[1] = alpha*(1.0-sat[0]);
  sat[2] = (1-alpha)*(1.0-sat[0]);
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = alpha *(1-sat[0]);
    sat[2] = (1-alpha) *(1-sat[0]);
  }
}


TEST_F( RelPermTest, numericalDerivatives_VanGenuchtenBakerRelPermTwoPhase )
{
  initialize( makeVanGenuchtenBakerRelPermTwoPhase( "relPerm", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.3;
  real64 const end_sat   = 0.7;
  real64 const dS = 1e-1;
  real64 const alpha = 0.4;
  array1d< real64 > sat( 2 );
  sat[0] = start_sat;
  sat[1] = alpha*(1.0-sat[0]);
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1-sat[0];
  }
}

TEST_F( RelPermTest, numericalDerivatives_VanGenuchtenBakerRelPermThreePhase )
{
  initialize( makeVanGenuchtenBakerRelPermThreePhase( "relPerm", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.3;
  real64 const end_sat   = 0.7;
  real64 const dS = 1e-1;
  real64 const alpha = 0.4;
  array1d< real64 > sat( 3 );
  sat[0] = start_sat;
  sat[1] = alpha*(1.0-sat[0]);
  sat[2] = (1-alpha)*(1.0-sat[0]);
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = alpha *(1-sat[0]);
    sat[2] = (1-alpha) *(1-sat[0]);
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
