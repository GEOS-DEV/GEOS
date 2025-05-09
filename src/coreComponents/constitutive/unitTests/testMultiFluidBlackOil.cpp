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
 * @file testMultiFluidBlackOil.cpp
 */

 #include "FluidModelTest.hpp"
 #include "constitutive/fluid/multifluid/blackOil/BlackOilFluid.hpp"
 #include "common/initializeEnvironment.hpp"

using namespace geos::constitutive;

namespace geos
{
namespace testing
{

static constexpr char const * pvdgTableContent = "# Pg(Pa) Bg(m3/sm3) Visc(Pa.s)\n"
                                                 "3000000  0.04234  0.00001344\n"
                                                 "6000000  0.02046  0.0000142\n"
                                                 "9000000  0.01328  0.00001526\n"
                                                 "12000000 0.00977  0.0000166\n"
                                                 "15000000 0.00773  0.00001818\n"
                                                 "18000000 0.006426 0.00001994\n"
                                                 "21000000 0.005541 0.00002181\n"
                                                 "24000000 0.004919 0.0000237\n"
                                                 "27000000 0.004471 0.00002559 -- this is a comment\n"
                                                 "29500000 0.004194 0.00002714\n"
                                                 "31000000 0.004031 0.00002806\n"
                                                 "33000000 0.00391  0.00002832\n"
                                                 "53000000 0.003868 0.00002935";

static const char * pvtoTableContent = "# Rs[sm3/sm3]\tPbub[Pa]\tBo[m3/sm3]\tVisc(Pa.s)\n"
                                       "\n"
                                       "  2\t            2000000\t    1.02\t    0.000975\n"
                                       "  5\t            5000000\t    1.03\t    0.00091\n"
                                       " 10\t            10000000\t1.04\t    0.00083\n"
                                       " 15\t            20000000\t1.05\t    0.000695\n"
                                       "                90000000\t1.03\t    0.000985  -- some line comment\n"
                                       " 30\t            30000000\t1.07\t    0.000594\n"
                                       " 40\t            40000000\t1.08\t    0.00051\n"
                                       "                50000000\t1.07\t    0.000549  -- another one\n"
                                       "                90000000\t1.06\t    0.00074\n"
                                       " 50\t            50000000.7\t1.09\t    0.000449\n"
                                       "                90000000.7\t1.08\t    0.000605";

static const char * pvtwTableContent = "#\tPref[Pa]\tBw[m3/sm3]\tCp[1/Pa]\t    Visc[Pa.s]\n"
                                       "\t30600000.1\t1.03\t\t0.00000000041\t0.0003";

template< typename USE_MASS_TYPE >
class MultiFluidBlackOilTestFixture : public FluidModelTest< BlackOilFluid, 3, 3 >
{
public:
  using Base = FluidModelTest< BlackOilFluid, 3, 3 >;

  static constexpr bool USE_MASS = USE_MASS_TYPE::value;
public:
  MultiFluidBlackOilTestFixture()
  {
    writeTableToFile( pvtoFileName, pvtoTableContent );
    writeTableToFile( pvdgFileName, pvdgTableContent );
    writeTableToFile( pvtwFileName, pvtwTableContent );

    Base::createFluid( getFluidName(), []( BlackOilFluid & fluid ){
      fillPhysicalProperties( fluid );
      fluid.setMassFlag( USE_MASS );
    } );
  }

  ~MultiFluidBlackOilTestFixture() override
  {
    removeFile( pvtoFileName );
    removeFile( pvdgFileName );
    removeFile( pvtwFileName );
  }

  static string getFluidName();

private:
  static void fillPhysicalProperties( BlackOilFluid & fluid );
  static constexpr const char * pvtoFileName = "pvto.txt";
  static constexpr const char * pvdgFileName = "pvdg.txt";
  static constexpr const char * pvtwFileName = "pvtw.txt";
};

template< typename USE_MASS_TYPE >
string MultiFluidBlackOilTestFixture< USE_MASS_TYPE >::getFluidName()
{
  return GEOS_FMT( "fluid{}", (USE_MASS ? "Mass" : "Molar"));
}

template< typename USE_MASS_TYPE >
void MultiFluidBlackOilTestFixture< USE_MASS_TYPE >::fillPhysicalProperties( BlackOilFluid & fluid )
{
  string_array & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames = {"oil", "water", "gas"};

  string_array & compNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  compNames = {"oil", "water", "gas"};

  array1d< real64 > & molarWgt = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
  fill( molarWgt, Feed< 3 >{114e-3, 18e-3, 16e-3} );

  array1d< real64 > & surfaceDens = fluid.getReference< array1d< real64 > >( BlackOilFluidBase::viewKeyStruct::surfacePhaseMassDensitiesString() );
  fill( surfaceDens, Feed< 3 >{800.0, 1022.0, 0.9907} );

  path_array & tableNames = fluid.getReference< path_array >( BlackOilFluidBase::viewKeyStruct::tableFilesString() );
  tableNames.emplace_back( Path( pvtoFileName ));
  tableNames.emplace_back( Path( pvtwFileName ));
  tableNames.emplace_back( Path( pvdgFileName ));
}

//using TestTypes = ::testing::Types< std::true_type, std::false_type >;
using TestTypes = ::testing::Types< std::true_type, std::false_type >;
class NameGenerator
{
public:
  template< typename T >
  static std::string GetName( int index )
  {
    return GEOS_FMT( "BlackOil_{}_{}", MultiFluidBlackOilTestFixture< T >::getFluidName(), index );
  }
};

TYPED_TEST_SUITE( MultiFluidBlackOilTestFixture, TestTypes, NameGenerator );

TYPED_TEST( MultiFluidBlackOilTestFixture, numericalDerivatives )
{
  BlackOilFluid * fluid = this->getFluid( this->getFluidName() );

  real64 constexpr eps = 1.0e-6;

  auto const samples = {
    Feed< 3 >{ 0.49500, 0.50000, 0.00500 },
    Feed< 3 >{ 0.60000, 0.20000, 0.20000 }
  };
  auto const pressures = { 2.000e5, 312.000e5, 800.000e5  };

  for( real64 const pressure : pressures )
  {
    for( auto const & sample : samples )
    {
      typename TestFixture::Base::TestPoint const data ( pressure, 297.15, sample );
      TestFixture::Base::testNumericalDerivatives( fluid, data, eps );
    }
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
