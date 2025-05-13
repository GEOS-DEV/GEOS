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
 * @file testTwoPhaseImmiscibleFluid.cpp
 */

#include "FluidModelTest.hpp"
#include "constitutive/fluid/twophaseimmisciblefluid/TwoPhaseImmiscibleFluid.hpp"
#include "common/initializeEnvironment.hpp"

using namespace geos::constitutive;

namespace geos
{
namespace testing
{

static constexpr char const * tableContentPhase0 = "# Pg(Pa) Dens(kg/m3) Visc(Pa.s)\n"
                                                   "0.22     0.00603  40203\n"
                                                   "0.3      0.04224  31311\n"
                                                   "0.5      0.15011  22423\n"
                                                   "0.6      0.22423  15011\n"
                                                   "0.8      0.31311  4224\n"
                                                   "1.0      0.40203  603";

static constexpr char const * tableContentPhase1 = "# Pg(Pa) Dens(kg/m3) Visc(Pa.s)\n"
                                                   "1.22     0.00603  0.22\n"
                                                   "1.3      0.04224  0.22\n"
                                                   "1.5      0.15011  0.22\n"
                                                   "1.6      0.22423  0.22\n"
                                                   "1.8      0.31311  0.22\n"
                                                   "2.0      0.40203  0.22";

template< bool FROM_TABLE >
class TwoPhaseImmiscibleFluidTest : public FluidModelTest< TwoPhaseImmiscibleFluid, 2 >
{
public:
  using Base = FluidModelTest< TwoPhaseImmiscibleFluid, 2 >;

public:
  TwoPhaseImmiscibleFluidTest()
  {
    if constexpr (!FROM_TABLE)
    {
      writeTableToFile( "phase0.txt", tableContentPhase0 );
      writeTableToFile( "phase1.txt", tableContentPhase1 );
    }

    Base::createFluid( getFluidName(), []( TwoPhaseImmiscibleFluid & fluid ){
      makeTwoPhaseImmiscibleFluid( fluid );
    } );
  }

  ~TwoPhaseImmiscibleFluidTest()
  {
    if constexpr (!FROM_TABLE)
    {
      removeFile( "phase0.txt" );
      removeFile( "phase1.txt" );
    }
  }

  void testDerivatives( TwoPhaseImmiscibleFluid * fluid,
                        real64 const pressure,
                        real64 const perturbParameter,
                        real64 const relTol,
                        real64 const absTol = std::numeric_limits< real64 >::max() )
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;

    integer constexpr actualDof = 1; // Pressure only

    // Number of execution points actual value plus left and right pertubations for each variable
    integer constexpr size = 2*actualDof + 1;
    m_parent.resize( size );
    fluid->allocateConstitutiveData( m_parent, 1 );

    array1d< real64 > pressureArray( size );
    array1d< real64 > deltaArray( size );

    // Name the degrees of freedom
    string_array dofNames( actualDof );
    string_array const & phaseNames = fluid->phaseNames();

    for( integer i = 0; i < size; ++i )
    {
      pressureArray[i] = pressure;
    }

    string const testValues = GEOS_FMT( "Pressure: {}", pressureArray[0] );

    real64 const dP = perturbParameter * (pressure + perturbParameter);
    deltaArray[2*Deriv::dP+1] = dP;
    deltaArray[2*Deriv::dP+2] = -dP;
    pressureArray[2*Deriv::dP+1] += deltaArray[2*Deriv::dP+1];
    pressureArray[2*Deriv::dP+2] += deltaArray[2*Deriv::dP+2];
    dofNames[Deriv::dP] = "pressure";

    auto const pressureView = pressureArray.toViewConst();

    auto fluidWrapper = fluid->createKernelWrapper();
    forAll< parallelHostPolicy >( size, [fluidWrapper,
                                         pressureView]
                                  GEOS_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pressureView[k] );
      }
    } );

    // Bring everything back to host
    forAll< serialPolicy >( 1, [fluidWrapper] ( localIndex const )
    {
      /* Do nothing */
    } );

    for( integer phaseIndex = 0; phaseIndex < numPhase; phaseIndex++ )
    {
      Base::testDerivatives( GEOS_FMT( "Phase density ({})", phaseNames[phaseIndex] ),
                             testValues,
                             fluid->phaseDensity().toViewConst(),
                             fluid->dPhaseDensity().toViewConst(),
                             deltaArray.toSliceConst(),
                             1.0,
                             dofNames,
                             relTol,
                             absTol,
                             phaseIndex );

      Base::testDerivatives( GEOS_FMT( "Phase viscosity ({})", phaseNames[phaseIndex] ),
                             testValues,
                             fluid->phaseViscosity().toViewConst(),
                             fluid->dPhaseViscosity().toViewConst(),
                             deltaArray.toSliceConst(),
                             1.0,
                             dofNames,
                             relTol,
                             absTol,
                             phaseIndex );
    }
  }   // void testDerivatives

  static string getFluidName()
  {
    return GEOS_FMT( "fluid{}", (FROM_TABLE ? "Tables" : "Files"));
  }

private:
  static TwoPhaseImmiscibleFluid * makeTwoPhaseImmiscibleFluid( TwoPhaseImmiscibleFluid & fluid );
};  // class TwoPhaseImmiscibleFluidTest

template<>
TwoPhaseImmiscibleFluid * TwoPhaseImmiscibleFluidTest< true >::makeTwoPhaseImmiscibleFluid( TwoPhaseImmiscibleFluid & fluid )
{
  // 1D table with linear interpolation
  localIndex constexpr Naxis = 6;
  localIndex constexpr NaxisSingle = 1;

  real64_array densityCoordPhase0;
  fill( densityCoordPhase0, Feed< Naxis >{ 0.22, 0.3, 0.5, 0.6, 0.8, 1.0 } );
  real64_array densityValuesPhase0;
  fill( densityValuesPhase0, Feed< Naxis >{ 0.00603, 0.04224, 0.04224, 0.22423, 0.31311, 0.40203 } );

  real64_array densityCoordPhase1;
  fill( densityCoordPhase1, Feed< Naxis >{ 1.22, 1.3, 1.5, 1.6, 1.8, 2.0 } );
  real64_array densityValuesPhase1;
  fill( densityValuesPhase1, Feed< Naxis >{ 0.00603, 0.04224, 0.04224, 0.22423, 0.31311, 0.40203 } );


  real64_array viscosityCoordPhase0;
  fill( viscosityCoordPhase0, Feed< Naxis >{ 0.22, 0.3, 0.5, 0.6, 0.8, 1.0 } );
  real64_array viscosityValuesPhase0;
  fill( viscosityValuesPhase0, Feed< Naxis >{ 40203, 31311, 22423, 15011, 4224, 603 } );

  real64_array viscosityCoordPhase1;
  fill( viscosityCoordPhase1, Feed< NaxisSingle >{ 0.22 } );
  real64_array viscosityValuesPhase1;
  fill( viscosityValuesPhase1, Feed< NaxisSingle >{ 45 } );

  createTable( "densityTablePhase0", densityCoordPhase0, densityValuesPhase0 );
  createTable( "densityTablePhase1", densityCoordPhase1, densityValuesPhase1 );
  createTable( "viscosityTablePhase0", viscosityCoordPhase0, viscosityValuesPhase0 );
  createTable( "viscosityTablePhase1", viscosityCoordPhase1, viscosityValuesPhase1 );

  // 2) Set up the constitutive model

  string_array & phaseNames = fluid.getReference< string_array >( TwoPhaseImmiscibleFluid::viewKeyStruct::phaseNamesString() );
  phaseNames.emplace_back( "oil" );
  phaseNames.emplace_back( "water" );

  string_array & densityTableNames = fluid.getReference< string_array >( TwoPhaseImmiscibleFluid::viewKeyStruct::densityTableNamesString() );
  densityTableNames.emplace_back( "densityTablePhase0" );
  densityTableNames.emplace_back( "densityTablePhase1" );

  string_array & viscosityTableNames = fluid.getReference< string_array >( TwoPhaseImmiscibleFluid::viewKeyStruct::viscosityTableNamesString() );
  viscosityTableNames.emplace_back( "viscosityTablePhase0" );
  viscosityTableNames.emplace_back( "viscosityTablePhase1" );

  return &fluid;
}

template<>
TwoPhaseImmiscibleFluid * TwoPhaseImmiscibleFluidTest< false >::makeTwoPhaseImmiscibleFluid( TwoPhaseImmiscibleFluid & fluid )
{
  string_array & phaseNames = fluid.getReference< string_array >( TwoPhaseImmiscibleFluid::viewKeyStruct::phaseNamesString() );
  phaseNames.emplace_back( "oil" );
  phaseNames.emplace_back( "water" );

  path_array & tableNames = fluid.getReference< path_array >( TwoPhaseImmiscibleFluid::viewKeyStruct::tableFilesString() );
  tableNames.emplace_back( "phase0.txt" );
  tableNames.emplace_back( "phase1.txt" );

  return &fluid;
}

using TwoPhaseImmiscibleFluidTestFromFiles = TwoPhaseImmiscibleFluidTest< false >;
using TwoPhaseImmiscibleFluidTestFromTables = TwoPhaseImmiscibleFluidTest< true >;

TEST_F( TwoPhaseImmiscibleFluidTestFromTables, testNumericalDerivative_initFromTables )
{
  auto * fluid = getFluid( this->getFluidName() );
  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon());
  real64 constexpr relTol = 1.0e-8;
  real64 constexpr absTol = 1.0e-8;

  for( real64 const pressure : { 0.55, 1.0, 10.0 } )
  {
    testDerivatives( fluid, pressure, eps, relTol, absTol );
  }
}

TEST_F( TwoPhaseImmiscibleFluidTestFromFiles, testNumericalDerivative_initFromFiles )
{
  auto * fluid = getFluid( this->getFluidName() );
  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon());
  real64 constexpr relTol = 1.0e-8;
  real64 constexpr absTol = 1.0e-8;

  for( real64 const pressure : { 0.55, 1.0, 10.0 } )
  {
    testDerivatives( fluid, pressure, eps, relTol, absTol );
  }
}

} // testing
} // geos

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::setupEnvironment( argc, argv );

  int const result = RUN_ALL_TESTS();

  geos::cleanupEnvironment();

  return result;
}
