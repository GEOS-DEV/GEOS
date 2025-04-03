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

#include "constitutive/fluid/twophaseimmisciblefluid/TwoPhaseImmiscibleFluid.hpp"
#include "constitutive/fluid/twophaseimmisciblefluid/TwoPhaseImmiscibleFluidFields.hpp"

// Only for fill
#include "unitTests/constitutiveTests/MultiFluidTest.hpp"

// Only for initializeTable
#include "unitTests/constitutiveTests/constitutiveTestHelpers.hpp"


using namespace geos;
using namespace geos::testing;
using namespace geos::constitutive;
using namespace geos::dataRepository;   /// Only for group definition


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
class TwoPhaseImmiscibleFluidTest : public ConstitutiveTestBase< TwoPhaseImmiscibleFluid >
{
public:

  TwoPhaseImmiscibleFluidTest()
  {
    if constexpr (!FROM_TABLE)
    {
      writeTableToFile( "phase0.txt", tableContentPhase0 );
      writeTableToFile( "phase1.txt", tableContentPhase1 );
    }

    m_parent.resize( 1 );
    string const fluidName = GEOS_FMT( "fluid{}", (FROM_TABLE ? "Tables" : "Files"));
    m_model = makeTwoPhaseImmiscibleFluid( fluidName, m_parent );

    m_parent.initialize();
    m_parent.initializePostInitialConditions();
  }

  ~TwoPhaseImmiscibleFluidTest()
  {
    if constexpr (!FROM_TABLE)
    {
      removeFile( "phase0.txt" );
      removeFile( "phase1.txt" );
    }
  }

  constitutive::TwoPhaseImmiscibleFluid & getFluid() const { return *m_model; }

  dataRepository::Group & getParent() { return m_parent; }


  void testDerivatives( constitutive::TwoPhaseImmiscibleFluid & fluid,
                        dataRepository::Group * parent,
                        real64 const pressure,
                        real64 const perturbParameter,
                        real64 const relTol,
                        real64 const absTol = std::numeric_limits< real64 >::max() )
  {
    auto const & phaseNames = fluid.getReference< string_array >( TwoPhaseImmiscibleFluid::viewKeyStruct::phaseNamesString() );

    // create a clone of the fluid to run updates on
    string const fluidCopyName = fluid.getName() + "Copy";
    std::unique_ptr< constitutive::ConstitutiveBase > fluidCopyPtr = fluid.deliverClone( fluidCopyName, parent );
    constitutive::TwoPhaseImmiscibleFluid & fluidCopy = dynamicCast< constitutive::TwoPhaseImmiscibleFluid & >( *fluidCopyPtr );
    fluidCopy.initializePostSubGroups();

    fluid.allocateConstitutiveData( fluid.getParent(), 1 );
    fluidCopy.allocateConstitutiveData( fluid.getParent(), 1 );

    // extract data views from both fluids
    #define GET_FLUID_DATA( FLUID, TRAIT ) \
      FLUID.getReference< TRAIT::type >( TRAIT::key() )[0][0]

    constitutive::MultiFluidVarSlice< real64, 1, constitutive::multifluid::USD_PHASE - 2, constitutive::multifluid::USD_PHASE_DC - 2 > phaseVisc {
      GET_FLUID_DATA( fluid, fields::twophaseimmisciblefluid::phaseViscosity ),
      GET_FLUID_DATA( fluid, fields::twophaseimmisciblefluid::dPhaseViscosity )
    };

    constitutive::MultiFluidVarSlice< real64, 1, constitutive::multifluid::USD_PHASE - 2, constitutive::multifluid::USD_PHASE_DC - 2 > phaseDens {
      GET_FLUID_DATA( fluid, fields::twophaseimmisciblefluid::phaseDensity ),
      GET_FLUID_DATA( fluid, fields::twophaseimmisciblefluid::dPhaseDensity )
    };

    auto const & phaseDensCopy = GET_FLUID_DATA( fluidCopy, fields::twophaseimmisciblefluid::phaseDensity );
    auto const & phaseViscCopy = GET_FLUID_DATA( fluidCopy, fields::twophaseimmisciblefluid::phaseViscosity );
    #undef GET_FLUID_DATA


    // set the original fluid state to current
    constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
    {
      typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
      fluidWrapper.update( 0, 0, pressure );
    } );

    // now perturb variables and update the copied fluid's state
    constitutive::constitutiveUpdatePassThru( fluidCopy, [&] ( auto & castedFluid )
    {
      using Deriv = constitutive::multifluid::DerivativeOffset;

      typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

      // to be able to use the checkDerivative utility function, we have to invert the layout
      auto dPhaseDens     = invertLayout( phaseDens.derivs.toSliceConst(), 2, 1 );
      auto dPhaseVisc     = invertLayout( phaseVisc.derivs.toSliceConst(), 2, 1 );

      // update pressure and check derivatives
      real64 const dP = perturbParameter * (pressure + perturbParameter);
      fluidWrapper.update( 0, 0, pressure + dP );

      checkDerivative( phaseDensCopy.toSliceConst(), phaseDens.value.toSliceConst(), dPhaseDens[Deriv::dP].toSliceConst(),
                       dP, relTol, absTol, "phaseDens", "Pressure", phaseNames );
      checkDerivative( phaseViscCopy.toSliceConst(), phaseVisc.value.toSliceConst(), dPhaseVisc[Deriv::dP].toSliceConst(),
                       dP, relTol, absTol, "phaseVisc", "Pressure", phaseNames );
    } );
  }   // void testDerivatives


protected:
  static void writeTableToFile( string const & fileName, char const * content )
  {
    std::ofstream os( fileName );
    ASSERT_TRUE( os.is_open() );
    os << content;
    os.close();
  }

  static void removeFile( string const & fileName )
  {
    int const ret = std::remove( fileName.c_str() );
    ASSERT_TRUE( ret == 0 );
  }


private:
  static TwoPhaseImmiscibleFluid * makeTwoPhaseImmiscibleFluid( string const & name, Group & parent );

};  // class TwoPhaseImmiscibleFluidTest


template<>
TwoPhaseImmiscibleFluid * TwoPhaseImmiscibleFluidTest< true >::makeTwoPhaseImmiscibleFluid( string const & name, Group & parent )
{
  // 1D table with linear interpolation
  localIndex constexpr Naxis = 6;
  localIndex constexpr NaxisSingle = 1;

  array1d< real64_array > densityCoordPhase0( 1 );
  fill< Naxis >( densityCoordPhase0[0], { 0.22, 0.3, 0.5, 0.6, 0.8, 1.0 } );
  real64_array densityValuesPhase0;
  fill< Naxis >( densityValuesPhase0, { 0.00603, 0.04224, 0.04224, 0.22423, 0.31311, 0.40203 } );

  array1d< real64_array > densityCoordPhase1( 1 );
  fill< Naxis >( densityCoordPhase1[0], { 1.22, 1.3, 1.5, 1.6, 1.8, 2.0 } );
  real64_array densityValuesPhase1;
  fill< Naxis >( densityValuesPhase1, { 0.00603, 0.04224, 0.04224, 0.22423, 0.31311, 0.40203 } );


  array1d< real64_array > viscosityCoordPhase0( 1 );
  fill< Naxis >( viscosityCoordPhase0[0], { 0.22, 0.3, 0.5, 0.6, 0.8, 1.0 } );
  real64_array viscosityValuesPhase0;
  fill< Naxis >( viscosityValuesPhase0, { 40203, 31311, 22423, 15011, 4224, 603 } );

  array1d< real64_array > viscosityCoordPhase1( 1 );
  fill< NaxisSingle >( viscosityCoordPhase1[0], { 0.22 } );
  real64_array viscosityValuesPhase1;
  fill< NaxisSingle >( viscosityValuesPhase1, { 45 } );

  initializeTable( "densityTablePhase0", densityCoordPhase0, densityValuesPhase0 );
  initializeTable( "densityTablePhase1", densityCoordPhase1, densityValuesPhase1 );
  initializeTable( "viscosityTablePhase0", viscosityCoordPhase0, viscosityValuesPhase0 );
  initializeTable( "viscosityTablePhase1", viscosityCoordPhase1, viscosityValuesPhase1 );


  // 2) Set up the constitutive model
  TwoPhaseImmiscibleFluid & fluid = parent.registerGroup< TwoPhaseImmiscibleFluid >( name );

  string_array & phaseNames = fluid.getReference< string_array >( TwoPhaseImmiscibleFluid::viewKeyStruct::phaseNamesString() );
  phaseNames.emplace_back( "oil" );
  phaseNames.emplace_back( "water" );

  string_array & densityTableNames = fluid.getReference< string_array >( TwoPhaseImmiscibleFluid::viewKeyStruct::densityTableNamesString() );
  densityTableNames.emplace_back( "densityTablePhase0" );
  densityTableNames.emplace_back( "densityTablePhase1" );

  string_array & viscosityTableNames = fluid.getReference< string_array >( TwoPhaseImmiscibleFluid::viewKeyStruct::viscosityTableNamesString() );
  viscosityTableNames.emplace_back( "viscosityTablePhase0" );
  viscosityTableNames.emplace_back( "viscosityTablePhase1" );

  fluid.postInputInitializationRecursive();
  return &fluid;
}


template<>
TwoPhaseImmiscibleFluid * TwoPhaseImmiscibleFluidTest< false >::makeTwoPhaseImmiscibleFluid( string const & name, Group & parent )
{
  TwoPhaseImmiscibleFluid & fluid = parent.registerGroup< TwoPhaseImmiscibleFluid >( name );

  path_array & tableNames = fluid.getReference< path_array >( TwoPhaseImmiscibleFluid::viewKeyStruct::tableFilesString() );
  fill< 2 >( tableNames, {"phase0.txt", "phase1.txt"} );

  fluid.postInputInitializationRecursive();
  return &fluid;
}



using TwoPhaseImmiscibleFluidTestFromFiles = TwoPhaseImmiscibleFluidTest< false >;
using TwoPhaseImmiscibleFluidTestFromTables = TwoPhaseImmiscibleFluidTest< true >;


TEST_F( TwoPhaseImmiscibleFluidTestFromTables, testNumericalDerivative_initFromTables )
{
  auto & fluid = getFluid();
  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon());
  real64 constexpr relTol = 1.0e-8;
  real64 constexpr absTol = 1.0e-8;

  for( real64 const pressure : { 0.55, 1.0, 10.0 } )
  {
    testDerivatives( fluid, &getParent(), pressure, eps, relTol, absTol );
  }
}


TEST_F( TwoPhaseImmiscibleFluidTestFromFiles, testNumericalDerivative_initFromFiles )
{
  auto & fluid = getFluid();
  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon());
  real64 constexpr relTol = 1.0e-8;
  real64 constexpr absTol = 1.0e-8;

  for( real64 const pressure : { 0.55, 1.0, 10.0 } )
  {
    testDerivatives( fluid, &getParent(), pressure, eps, relTol, absTol );
  }
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  conduit::Node conduitNode;
  dataRepository::Group rootNode( "root", conduitNode );
  FunctionManager functionManager( "FunctionManager", &rootNode );

  int const result = RUN_ALL_TESTS();

  return result;
}
