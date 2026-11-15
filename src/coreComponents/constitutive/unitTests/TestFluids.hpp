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
 * @file TestFluids.hpp
 * @brief Library of named fluid mixtures shared by the constitutive unit tests.
 *
 * Each mixture is identified by a value of @c TestFluidName. A new fluid is created with
 * @code
 *   auto fluid = createTestFluid< TestFluidName::SourGas >();
 * @endcode
 * Names describe the kind of fluid each mixture represents. The components, in order, and whether
 * non-zero binary interaction coefficients (kij) are set are listed with each name. Each definition
 * also provides a list of representative compositions in @c testSamples.
 */

#ifndef GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUIDS_HPP_
#define GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUIDS_HPP_

#include "TestFluid.hpp"

namespace geos
{

namespace testing
{

/**
 * @brief Names of the predefined test fluids
 */
enum class TestFluidName : integer
{
  // 2 components
  CarbonStorage,        ///< CO2 dissolving in water: H2O, CO2 (with kij)
  CO2Solvent,           ///< CO2 as a miscible solvent for a light oil: CO2, C5H12 (with kij)
  RichGas,              ///< Natural gas with natural gas liquids: CH4, C3H8 (with kij)
  OilWater,             ///< Water with a heavy oil: H2O, C10H22

  // 3 components
  HydrogenStorage,      ///< Hydrogen stored in water-bearing rock, with sour gas: H2O, H2S, H2
  GasOilWater,          ///< Simple gas, oil and water system: CH4, C10H22, H2O (with kij)

  // 4 components
  EnhancedOilRecovery,  ///< CO2 injected into an oil with water and gas: CO2, H2O, CH4, C5H12
  AquiferGasStorage,    ///< Natural gas stored in an aquifer: N2, CH4, CO2, H2O (with kij)
  SourGas,              ///< Natural gas with H2S and water: CH4, CO2, H2S, H2O (with kij)

  // 6 components
  SourNaturalGas,       ///< Fuller sour natural gas with water: CH4, C2H6, CO2, H2S, N2, H2O (with kij)

  // 9 components
  VolatileOil,          ///< Volatile oil with water: H2O, CO2, N2, CH4, C2H6, C3H8, C4H10, C5H12, C10H22 (with kij)
};

/**
 * @brief Definition of a named test fluid
 * @tparam NAME the name of the fluid
 *
 * Each specialisation provides:
 *  - @c numComps the number of components
 *  - @c components the list of components (indices into @c Fluid)
 *  - @c binaryCoefficients the lower triangle of the binary interaction coefficients, row by row
 *    (empty if all coefficients are zero)
 *  - @c testSamples a list of representative compositions
 */
template< TestFluidName NAME >
struct TestFluidDefinition;

/**
 * @brief Number of components in a named test fluid
 */
template< TestFluidName NAME >
constexpr integer getNumberOfComponents = TestFluidDefinition< NAME >::numComps;

/**
 * @brief Create a named test fluid
 */
template< TestFluidName NAME >
std::unique_ptr< TestFluid< getNumberOfComponents< NAME > > > createTestFluid()
{
  using Definition = TestFluidDefinition< NAME >;
  constexpr integer numComps = Definition::numComps;
  constexpr integer numBics = static_cast< integer >( Definition::binaryCoefficients.size() );
  static_assert( numBics == 0 || numBics == numComps*(numComps-1)/2,
                 "Binary coefficients must be empty or contain the full lower triangle" );

  auto fluid = TestFluid< numComps >::create( Definition::components );
  if constexpr ( 0 < numBics )
  {
    fluid->setBinaryCoefficients( Definition::binaryCoefficients );
  }
  return fluid;
}

/**
 * @brief Get the test samples of a named test fluid
 */
template< TestFluidName NAME >
constexpr auto const & getTestSamples()
{
  return TestFluidDefinition< NAME >::testSamples;
}

/// Type for the samples
template< integer NC, std::size_t NS >
using TestSamples = std::array< Feed< NC >, NS >;

/// Type for no binary coefficients
using NoBinaryCoefficients = std::array< real64, 0 >;

/* UNCRUSTIFY-OFF */

//-------------------------------------------------------------------------------------------------
// 2 components
//-------------------------------------------------------------------------------------------------

template<>
struct TestFluidDefinition< TestFluidName::CarbonStorage >
{
  static constexpr integer numComps = 2;
  static constexpr std::array< integer, numComps > components{ Fluid::H2O, Fluid::CO2 };
  static constexpr std::array< real64, 1 > binaryCoefficients{ 0.1896 };
  static constexpr real64 e = 1.0e-7;
  static constexpr TestSamples< numComps, 4 > testSamples{{
    {  0.5,   0.5},
    {  0.0,   1.0},
    {    e, 1.0-e},
    {1.0-e,     e}
  }};
};

template<>
struct TestFluidDefinition< TestFluidName::CO2Solvent >
{
  static constexpr integer numComps = 2;
  static constexpr std::array< integer, numComps > components{ Fluid::CO2, Fluid::C5H12 };
  static constexpr std::array< real64, 1 > binaryCoefficients{ 0.025 };
  static constexpr real64 e = 1.0e-7;
  static constexpr TestSamples< numComps, 4 > testSamples{{
    {  0.5,   0.5},
    {  0.0,   1.0},
    {    e, 1.0-e},
    {1.0-e,     e}
  }};
};

template<>
struct TestFluidDefinition< TestFluidName::RichGas >
{
  static constexpr integer numComps = 2;
  static constexpr std::array< integer, numComps > components{ Fluid::CH4, Fluid::C3H8 };
  static constexpr std::array< real64, 1 > binaryCoefficients{ 0.1 };
  static constexpr real64 e = 1.0e-7;
  static constexpr TestSamples< numComps, 4 > testSamples{{
    {  0.5,   0.5},
    {  0.0,   1.0},
    {    e, 1.0-e},
    {1.0-e,     e}
  }};
};

template<>
struct TestFluidDefinition< TestFluidName::OilWater >
{
  static constexpr integer numComps = 2;
  static constexpr std::array< integer, numComps > components{ Fluid::H2O, Fluid::C10H22 };
  static constexpr NoBinaryCoefficients binaryCoefficients{};
  static constexpr TestSamples< numComps, 3 > testSamples{{
    {0.995, 0.005},
    {1.000, 0.000},
    {0.002, 0.998}
  }};
};

//-------------------------------------------------------------------------------------------------
// 3 components
//-------------------------------------------------------------------------------------------------

template<>
struct TestFluidDefinition< TestFluidName::HydrogenStorage >
{
  static constexpr integer numComps = 3;
  static constexpr std::array< integer, numComps > components{ Fluid::H2O, Fluid::H2S, Fluid::H2 };
  static constexpr NoBinaryCoefficients binaryCoefficients{};
  static constexpr real64 e = 1.0e-7;
  static constexpr TestSamples< numComps, 4 > testSamples{{
    {    0.4, 0.3, 0.3},
    {    1.0, 0.0, 0.0},
    {1.0-2*e,   e,   e},
    {    0.0, 0.5, 0.5}
  }};
};

template<>
struct TestFluidDefinition< TestFluidName::GasOilWater >
{
  static constexpr integer numComps = 3;
  static constexpr std::array< integer, numComps > components{ Fluid::CH4, Fluid::C10H22, Fluid::H2O };
  static constexpr std::array< real64, 3 > binaryCoefficients{ 0.25, 0.0, 0.0 };
  static constexpr real64 e = 1.0e-7;
  static constexpr TestSamples< numComps, 4 > testSamples{{
    {0.3, 0.3,     0.4},
    {0.0, 0.0,     1.0},
    {  e,   e, 1.0-2*e},
    {0.5, 0.5,     0.0}
  }};
};

//-------------------------------------------------------------------------------------------------
// 4 components
//-------------------------------------------------------------------------------------------------

template<>
struct TestFluidDefinition< TestFluidName::EnhancedOilRecovery >
{
  static constexpr integer numComps = 4;
  static constexpr std::array< integer, numComps > components{ Fluid::CO2, Fluid::H2O, Fluid::CH4, Fluid::C5H12 };
  static constexpr NoBinaryCoefficients binaryCoefficients{};
  static constexpr TestSamples< numComps, 4 > testSamples{{
    {0.250, 0.250, 0.250, 0.250},
    {0.100, 0.700, 0.100, 0.100},
    {0.300, 0.001, 0.400, 0.299},
    {0.050, 0.050, 0.100, 0.800}
  }};
};

template<>
struct TestFluidDefinition< TestFluidName::AquiferGasStorage >
{
  static constexpr integer numComps = 4;
  static constexpr std::array< integer, numComps > components{ Fluid::N2, Fluid::CH4, Fluid::CO2, Fluid::H2O };
  static constexpr std::array< real64, 6 > binaryCoefficients{ 0.0, 0.0, 0.0, 0.4778, 0.4850, 0.1896 };
  static constexpr TestSamples< numComps, 4 > testSamples{{
    {   0.099,    0.001,    0.300,    0.600},
    {   0.350,    0.100,    0.350,    0.200},
    {   0.050,    0.800,    0.050,    0.100},
    {0.030933, 0.319683, 0.637861, 0.011523}
  }};
};

template<>
struct TestFluidDefinition< TestFluidName::SourGas >
{
  static constexpr integer numComps = 4;
  static constexpr std::array< integer, numComps > components{ Fluid::CH4, Fluid::CO2, Fluid::H2S, Fluid::H2O };
  static constexpr std::array< real64, 6 > binaryCoefficients{ 0.0, 0.0, 0.0, 0.4850, 0.1896, 0.1353 };
  static constexpr real64 e = 1.0e-7;
  static constexpr TestSamples< numComps, 4 > testSamples{{
    {0.3, 0.2, 0.2,     0.3},
    {0.0, 0.0, 0.0,     1.0},
    {0.0, 0.5, 0.0,     0.5},
    {  e,   e,   e, 1.0-3*e}
  }};
};

//-------------------------------------------------------------------------------------------------
// 6 components
//-------------------------------------------------------------------------------------------------

template<>
struct TestFluidDefinition< TestFluidName::SourNaturalGas >
{
  static constexpr integer numComps = 6;
  static constexpr std::array< integer, numComps > components{ Fluid::CH4, Fluid::C2H6, Fluid::CO2, Fluid::H2S, Fluid::N2, Fluid::H2O };
  static constexpr std::array< real64, 15 > binaryCoefficients{
    0.0,
    0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
    0.4850, 0.4920, 0.1896, 0.1353, 0.4778
  };
  static constexpr real64 e = 1.0e-7;
  static constexpr real64 f = 0.16666;
  static constexpr TestSamples< numComps, 5 > testSamples{{
    {        f,           f,           f,               f,           f,    1.0-5*f},
    {        e,           e,           e,               e,           e,    1.0-5*e},
    {        e,     1.0-5*e,           e,               e,           e,          e},
    { 0.260871,   0.0663935,   0.0906286,       0.0524321,  0.52755293, 0.00212187},
    {6.233e-08, 2.99791e-09, 8.00657e-06, 8.187700789e-05, 5.10942e-08,    0.99991}
  }};
};

//-------------------------------------------------------------------------------------------------
// 9 components
//-------------------------------------------------------------------------------------------------

template<>
struct TestFluidDefinition< TestFluidName::VolatileOil >
{
  static constexpr integer numComps = 9;
  static constexpr std::array< integer, numComps > components{
    Fluid::H2O, Fluid::CO2, Fluid::N2, Fluid::CH4, Fluid::C2H6, Fluid::C3H8, Fluid::C4H10, Fluid::C5H12, Fluid::C10H22
  };
  // The H2O column (first entry of each row) uses the Soreide-Whitson (1992) gas-phase values
  // (Table 5) where available. No values are given there for C5H12 and C10H22.
  static constexpr std::array< real64, 36 > binaryCoefficients{
    0.1896,                                                     // CO2
    0.4778, 0.003732,                                           // N2
    0.4850, 0.01,  0.0,                                         // CH4
    0.4920, 0.01,  0.0,     0.0,                                // C2H6
    0.5525, 0.0,   0.0,     0.0,  0.0,                          // C3H8
    0.5091, 0.0,   0.0,     0.0,  0.0,  0.0,                    // C4H10
    0.01,   0.0,   0.028,   0.01, 0.01, 0.0, 0.0,               // C5H12
    0.01,   0.0,   0.04532, 0.01, 0.01, 0.0, 0.0, 0.0           // C10H22
  };
  static constexpr TestSamples< numComps, 3 > testSamples{{
    {0.839010, 0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120},
    {0.015100, 0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.169200},
    {0.001898, 0.007026, 0.006161, 0.827761, 0.091046, 0.045353, 0.015026, 0.004474, 0.001256}
  }};
};

/* UNCRUSTIFY-ON */

}// testing
}// geos

#endif //GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUIDS_HPP_
