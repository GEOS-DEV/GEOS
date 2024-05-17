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
#include "TestFluid.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/RachfordRice.hpp"

using namespace geos;
using namespace geos::testing;
using namespace geos::constitutive;

template< int NC >
using RachfordRiceData = std::tuple<
  Feed< NC > const,         // k-values
  Feed< NC > const,         // composition
  real64 const              // expected vapour fraction
  >;

template< integer NUM_COMP >
class RachfordRiceTest : public ::testing::TestWithParam< RachfordRiceData< NUM_COMP > >
{
public:
  static constexpr real64 relTol = 1.0e-4;
  static constexpr integer numComps = NUM_COMP;
public:
  RachfordRiceTest() = default;
  ~RachfordRiceTest() = default;

  void testExpectedValue( RachfordRiceData< NUM_COMP > const & data )
  {
    stackArray1d< real64, numComps > kValues;
    stackArray1d< real64, numComps > composition;
    array1d< integer > presentComponents;
    TestFluid< NUM_COMP >::createArray( kValues, std::get< 0 >( data ));
    TestFluid< NUM_COMP >::createArray( composition, std::get< 1 >( data ));
    real64 const expectedVapourFraction = std::get< 2 >( data );

    calculatePresentComponents( composition.toSliceConst(), presentComponents );

    real64 const vapourFraction = RachfordRice::solve( kValues.toSliceConst(),
                                                       composition.toSliceConst(),
                                                       presentComponents.toSliceConst() );

    //for( integer ic = 0; ic < numComps; ++ic )
    //  std::cout << composition[ic] << " ";
    //std::cout << "| ";
    //for( integer ic : presentComponents )
    //  std::cout << ic << " ";
    //std::cout << "| " << std::scientific << vapourFraction << "\n";
    //GEOS_UNUSED_VAR( expectedVapourFraction, vapourFraction );

    checkRelativeError( vapourFraction, expectedVapourFraction, relTol );
  }

private:
  static integer calculatePresentComponents( arraySlice1d< real64 const > const & composition,
                                             array1d< integer > & presentComponents )
  {
    // Check for machine-zero feed values
    for( integer ic = 0; ic < numComps; ++ic )
    {
      if( MultiFluidConstants::epsilon < composition[ic] )
      {
        presentComponents.emplace_back( ic );
      }
    }
    return presentComponents.size();
  }
};

using RachfordRiceTest2 = RachfordRiceTest< 2 >;
using RachfordRiceTest4 = RachfordRiceTest< 4 >;
using RachfordRiceTest9 = RachfordRiceTest< 9 >;

TEST_P( RachfordRiceTest2, testExpectedValues )
{
  testExpectedValue( GetParam() );
}
TEST_P( RachfordRiceTest4, testExpectedValues )
{
  testExpectedValue( GetParam() );
}
TEST_P( RachfordRiceTest9, testExpectedValues )
{
  testExpectedValue( GetParam() );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------
/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P(
  RachfordRiceTest, RachfordRiceTest2,
  ::testing::Values(
    RachfordRiceData< 2 >{ {1.12223, 1.12223}, {0.500000, 0.500000}, 1.00000 },
    RachfordRiceData< 2 >{ {56.1091, 56.1091}, {0.500000, 0.500000}, 1.00000 },
    RachfordRiceData< 2 >{ {54.8660, 54.8660}, {0.100000, 0.900000}, 1.00000 },
    RachfordRiceData< 2 >{ {1.09733, 1.09733}, {0.100000, 0.900000}, 1.00000 },
    RachfordRiceData< 2 >{ {0.90000, 1.09733}, {1.e-10,  1.-1.e-10}, 1.00000 },
    RachfordRiceData< 2 >{ {1.09733, 0.90000}, {1.e-10,  1.-1.e-10}, 0.00000 },
    RachfordRiceData< 2 >{ {0.80000, 1.25000}, {0.500000, 0.500000}, 0.50000 }
  )
);

INSTANTIATE_TEST_SUITE_P(
  RachfordRiceTest, RachfordRiceTest4,
  ::testing::Values(
    RachfordRiceData< 4 >{ {17.4329, 0.000770753, 1.18694e-06, 0.00968376}, {0.1, 0.1, 0.1, 0.7}, 0.045999 },
    RachfordRiceData< 4 >{ {11.5506, 0.000210682, 1.57686e-08, 0.0442361}, {0.0984186, 0.297297, 0.593142, 0.0111427}, 0.0130259 }
  )
);

INSTANTIATE_TEST_SUITE_P(
  RachfordRiceTest, RachfordRiceTest9,
  ::testing::Values(
    RachfordRiceData< 9 >{
      {9.999990228775e-01, 7.866564714988e+00, 1.000001898362e+00, 3.465944574258e-01, 7.453388393000e-02, 2.040256983735e-02, 5.604465489973e-03, 1.519555065458e-03, 2.452615074973e-06},
      {5.000000000000e-01, 0.000000000000e+00, 5.000000000000e-01, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00},
      6.358674e+05
    },
    RachfordRiceData< 9 >{
      {9.999990228775e-01, 7.866564714988e+00, 1.000001898362e+00, 3.465944574258e-01, 7.453388393000e-02, 2.040256983735e-02, 5.604465489973e-03, 1.519555065458e-03, 2.452615074973e-06},
      {5.000050001000e-01, 0.000000000000e+00, 5.000000000000e-01, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00},
      2.483177e+05
    },
    RachfordRiceData< 9 >{
      {9.999370199314e-01, 7.866756794708e+00, 1.000002841977e+00, 3.466160938183e-01, 7.453979327738e-02, 2.040447363770e-02, 5.605068501236e-03, 1.519741440033e-03, 2.453098843789e-06},
      {5.000000000000e-01, 0.000000000000e+00, 5.000000000000e-01, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00},
      -2.599313e+05
    },
    RachfordRiceData< 9 >{
      {9.999370199314e-01, 7.866756794708e+00, 1.000002841977e+00, 3.466160938183e-01, 7.453979327738e-02, 2.040447363770e-02, 5.605068501236e-03, 1.519741440033e-03, 2.453098843789e-06},
      {5.000050001000e-01, 0.000000000000e+00, 5.000000000000e-01, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00},
      -1.679958e+05
    }
  )
);

/* UNCRUSTIFY-ON */
