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

// Source includes
#include "managers/initialization.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/relativePermeability/relativePermeabilitySelector.hpp"
#include "physicsSolvers/fluidFlow/unitTests/testCompFlowUtils.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::testing;
using namespace geosx::constitutive;
using namespace geosx::dataRepository;

void testNumericalDerivatives( RelativePermeabilityBase & relPerm,
                               arraySlice1d< real64 const > const & saturation,
                               real64 const perturbParameter,
                               real64 const relTol )
{
  localIndex const NP = relPerm.numFluidPhases();
  auto const & phases = relPerm.phaseNames();

  // create a clone of the rel perm to run updates on
  std::unique_ptr< ConstitutiveBase > relPermCopyPtr = relPerm.deliverClone( "fluidCopy", nullptr );
  RelativePermeabilityBase & relPermCopy = *relPermCopyPtr->groupCast< RelativePermeabilityBase * >();

  relPerm.allocateConstitutiveData( relPerm.getParent(), 1 );
  relPermCopy.allocateConstitutiveData( relPerm.getParent(), 1 );

  arrayView3d< real64 const > const phaseRelPerm = relPerm.phaseRelPerm();
  arrayView4d< real64 const > const dPhaseRelPerm_dSat = relPerm.dPhaseRelPerm_dPhaseVolFraction();
  arrayView3d< real64 const > const phaseRelPermCopy = relPermCopy.phaseRelPerm();

  // set the fluid state to current
  constitutive::constitutiveUpdatePassThru( relPerm, [&] ( auto & castedRelPerm )
  {
    typename TYPEOFREF( castedRelPerm ) ::KernelWrapper relPermWrapper = castedRelPerm.createKernelWrapper();
    relPermWrapper.update( 0, 0, saturation );
  } );

  // update saturation and check derivatives
  auto dPhaseRelPerm_dS = invertLayout( dPhaseRelPerm_dSat[ 0 ][ 0 ], NP, NP );

  array1d< real64 > satNew( NP );
  for( localIndex jp = 0; jp < NP; ++jp )
  {
    real64 const dS = perturbParameter * (saturation[jp] + perturbParameter);
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      satNew[ip] = saturation[ip];
    }
    satNew[jp] += dS;

    constitutive::constitutiveUpdatePassThru( relPermCopy, [&] ( auto & castedRelPerm )
    {
      typename TYPEOFREF( castedRelPerm ) ::KernelWrapper relPermWrapper = castedRelPerm.createKernelWrapper();
      relPermWrapper.update( 0, 0, satNew );
    } );

    string const var = "phaseVolFrac[" + phases[jp] + "]";
    checkDerivative( phaseRelPermCopy[ 0 ][ 0 ],
                     phaseRelPerm[ 0 ][ 0 ],
                     dPhaseRelPerm_dS[ jp ].toSliceConst(),
                     dS,
                     relTol,
                     "phaseRelPerm",
                     var,
                     phases );
  }
}

RelativePermeabilityBase * makeBrooksCoreyRelPerm( string const & name, Group * parent )
{
  auto relPerm = parent->registerGroup< BrooksCoreyRelativePermeability >( name );

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

  relPerm->postProcessInputRecursive();
  return relPerm;
}

RelativePermeabilityBase * makeBrooksCoreyBakerRelPermTwoPhase( string const & name, Group * parent )
{
  auto relPerm = parent->registerGroup< BrooksCoreyBakerRelativePermeability >( name );

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

  relPerm->postProcessInputRecursive();
  return relPerm;
}

RelativePermeabilityBase * makeBrooksCoreyBakerRelPermThreePhase( string const & name, Group * parent )
{
  auto relPerm = parent->registerGroup< BrooksCoreyBakerRelativePermeability >( name );

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

  relPerm->postProcessInputRecursive();
  return relPerm;
}

RelativePermeabilityBase * makeVanGenuchtenBakerRelPermTwoPhase( string const & name, Group * parent )
{
  auto relPerm = parent->registerGroup< VanGenuchtenBakerRelativePermeability >( name );

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

  relPerm->postProcessInputRecursive();
  return relPerm;
}

RelativePermeabilityBase * makeVanGenuchtenBakerRelPermThreePhase( string const & name, Group * parent )
{
  auto relPerm = parent->registerGroup< VanGenuchtenBakerRelativePermeability >( name );

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

  relPerm->postProcessInputRecursive();
  return relPerm;
}

RelativePermeabilityBase * makeTableRelPermTwoPhase( string const & name, Group * parent )
{
  FunctionManager * functionManager = &FunctionManager::FunctionManager::instance();

  // 1) First, define the tables

  // 1D table, various interpolation methods
  localIndex Naxis = 6;

  // Setup table
  array1d< real64_array > coordinates;
  coordinates.resize( 1 );
  coordinates[0].resize( Naxis );
  coordinates[0][0] = 0.0;
  coordinates[0][1] = 0.11;
  coordinates[0][2] = 0.23;
  coordinates[0][3] = 0.54;
  coordinates[0][4] = 0.85;
  coordinates[0][5] = 1.0;

  real64_array values( Naxis );
  for( localIndex i = 0; i < coordinates[0].size(); ++i )
  {
    values[i] = coordinates[0][i]*coordinates[0][i];
  }

  TableFunction * table_w = functionManager->createChild( "TableFunction", "water_swof" )->groupCast< TableFunction * >();
  table_w->setTableCoordinates( coordinates );
  table_w->setTableValues( values );
  table_w->reInitializeFunction();

  table_w->setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction * table_o = functionManager->createChild( "TableFunction", "oil_swof" )->groupCast< TableFunction * >();
  table_o->setTableCoordinates( coordinates );
  table_o->setTableValues( values );
  table_o->reInitializeFunction();

  table_o->setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 2) Then set up the constitutive model

  auto relPerm = parent->registerGroup< TableRelativePermeability >( name );

  auto & phaseNames = relPerm->getReference< string_array >( RelativePermeabilityBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 2 );
  phaseNames[0] = "oil"; phaseNames[1] = "water";

  auto & waterOilTableNames = relPerm->getReference< array1d< string > >( TableRelativePermeability::viewKeyStruct::waterOilRelPermTableNamesString );
  waterOilTableNames.resize( 2 );
  waterOilTableNames[0] = "water_swof"; waterOilTableNames[1] = "oil_swof";

  relPerm->postProcessInputRecursive();
  return relPerm;
}

RelativePermeabilityBase * makeTableRelPermThreePhase( string const & name, Group * parent )
{
  FunctionManager * functionManager = &FunctionManager::FunctionManager::instance();

  // 1) First, define the tables

  // 1D table, various interpolation methods
  localIndex Naxis = 6;

  // 1.a) First pair of phases (ow)

  // Setup table
  array1d< real64_array > coordinates;
  coordinates.resize( 1 );
  coordinates[0].resize( Naxis );
  coordinates[0][0] = 0.0;
  coordinates[0][1] = 0.11;
  coordinates[0][2] = 0.23;
  coordinates[0][3] = 0.54;
  coordinates[0][4] = 0.85;
  coordinates[0][5] = 1.0;

  real64_array values( Naxis );
  for( localIndex i = 0; i < coordinates[0].size(); ++i )
  {
    values[i] = coordinates[0][i]*coordinates[0][i];
  }

  TableFunction * table_ow_w = functionManager->createChild( "TableFunction", "water_swof" )->groupCast< TableFunction * >();
  table_ow_w->setTableCoordinates( coordinates );
  table_ow_w->setTableValues( values );
  table_ow_w->reInitializeFunction();

  table_ow_w->setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction * table_ow_o = functionManager->createChild( "TableFunction", "oil_swof" )->groupCast< TableFunction * >();
  table_ow_o->setTableCoordinates( coordinates );
  table_ow_o->setTableValues( values );
  table_ow_o->reInitializeFunction();

  table_ow_o->setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 1.a) Second pair of phases (og)

  coordinates[0].resize( Naxis );
  coordinates[0][0] = 0.0;
  coordinates[0][1] = 0.01;
  coordinates[0][2] = 0.23;
  coordinates[0][3] = 0.44;
  coordinates[0][4] = 0.83;
  coordinates[0][5] = 1.0;

  for( localIndex i = 0; i < coordinates[0].size(); ++i )
  {
    values[i] = coordinates[0][i]*coordinates[0][i]*coordinates[0][i];
  }

  TableFunction * table_og_g = functionManager->createChild( "TableFunction", "gas_sgof" )->groupCast< TableFunction * >();
  table_og_g->setTableCoordinates( coordinates );
  table_og_g->setTableValues( values );
  table_og_g->reInitializeFunction();

  table_og_g->setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction * table_og_o = functionManager->createChild( "TableFunction", "oil_sgof" )->groupCast< TableFunction * >();
  table_og_o->setTableCoordinates( coordinates );
  table_og_o->setTableValues( values );
  table_og_o->reInitializeFunction();

  table_og_o->setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 2) Then set up the constitutive model

  auto relPerm = parent->registerGroup< TableRelativePermeability >( name );

  auto & phaseNames = relPerm->getReference< string_array >( RelativePermeabilityBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil"; phaseNames[1] = "water"; phaseNames[2] = "gas";

  auto & waterOilTableNames = relPerm->getReference< array1d< string > >( TableRelativePermeability::viewKeyStruct::waterOilRelPermTableNamesString );
  waterOilTableNames.resize( 2 );
  waterOilTableNames[0] = "water_swof"; waterOilTableNames[1] = "oil_swof";

  auto & gasOilTableNames = relPerm->getReference< array1d< string > >( TableRelativePermeability::viewKeyStruct::gasOilRelPermTableNamesString );
  gasOilTableNames.resize( 2 );
  gasOilTableNames[0] = "gas_sgof"; gasOilTableNames[1] = "oil_sgof";

  relPerm->postProcessInputRecursive();
  return relPerm;
}

TEST( testRelPerm, numericalDerivatives_brooksCoreyRelPerm )
{
  auto parent = std::make_unique< Group >( "parent", nullptr );
  parent->resize( 1 );

  RelativePermeabilityBase * relperm = makeBrooksCoreyRelPerm( "relPerm", parent.get() );

  parent->initialize( parent.get() );
  parent->initializePostInitialConditions( parent.get() );

  // TODO test over a range of values
  array1d< real64 > sat( 4 );
  sat[0] = 0.7; sat[1] = 0.3;

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  testNumericalDerivatives( *relperm, sat, eps, tol );
}

TEST( testRelPerm, numericalDerivatives_BrooksCoreyBakerRelPermTwoPhase )
{
  auto parent = std::make_unique< Group >( "parent", nullptr );
  parent->resize( 1 );

  RelativePermeabilityBase * relperm = makeBrooksCoreyBakerRelPermTwoPhase( "relPerm", parent.get() );

  parent->initialize( parent.get() );
  parent->initializePostInitialConditions( parent.get() );

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
    testNumericalDerivatives( *relperm, sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1-sat[0];
  }
}

TEST( testRelPerm, numericalDerivatives_BrooksCoreyBakerRelPermThreePhase )
{
  auto parent = std::make_unique< Group >( "parent", nullptr );
  parent->resize( 1 );

  RelativePermeabilityBase * relperm = makeBrooksCoreyBakerRelPermThreePhase( "relPerm", parent.get() );

  parent->initialize( parent.get() );
  parent->initializePostInitialConditions( parent.get() );

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
    testNumericalDerivatives( *relperm, sat, eps, tol );
    sat[0] += dS;
    sat[1] = alpha *(1-sat[0]);
    sat[2] = (1-alpha) *(1-sat[0]);
  }
}


TEST( testRelPerm, numericalDerivatives_VanGenuchtenBakerRelPermTwoPhase )
{
  auto parent = std::make_unique< Group >( "parent", nullptr );
  parent->resize( 1 );

  RelativePermeabilityBase * relperm = makeVanGenuchtenBakerRelPermTwoPhase( "relPerm", parent.get() );

  parent->initialize( parent.get() );
  parent->initializePostInitialConditions( parent.get() );

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
    testNumericalDerivatives( *relperm, sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1-sat[0];
  }
}

TEST( testRelPerm, numericalDerivatives_VanGenuchtenBakerRelPermThreePhase )
{
  auto parent = std::make_unique< Group >( "parent", nullptr );
  parent->resize( 1 );

  RelativePermeabilityBase * relperm = makeVanGenuchtenBakerRelPermThreePhase( "relPerm", parent.get() );

  parent->initialize( parent.get() );
  parent->initializePostInitialConditions( parent.get() );

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
    testNumericalDerivatives( *relperm, sat, eps, tol );
    sat[0] += dS;
    sat[1] = alpha *(1-sat[0]);
    sat[2] = (1-alpha) *(1-sat[0]);
  }
}

TEST( testRelPerm, numericalDerivatives_TableRelPermTwoPhase )
{
  auto parent = std::make_unique< Group >( "parent", nullptr );
  parent->resize( 1 );

  RelativePermeabilityBase * relperm = makeTableRelPermTwoPhase( "relPerm", parent.get() );

  parent->initialize( parent.get() );
  parent->initializePostInitialConditions( parent.get() );

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
    testNumericalDerivatives( *relperm, sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1-sat[0];
  }
}

TEST( testRelPerm, numericalDerivatives_TableRelPermThreePhase )
{
  auto parent = std::make_unique< Group >( "parent", nullptr );
  parent->resize( 1 );

  RelativePermeabilityBase * relperm = makeTableRelPermThreePhase( "relPerm", parent.get() );

  parent->initialize( parent.get() );
  parent->initializePostInitialConditions( parent.get() );

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
    testNumericalDerivatives( *relperm, sat, eps, tol );
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
