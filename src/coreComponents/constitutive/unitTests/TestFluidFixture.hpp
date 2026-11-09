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
 * @file TestFluidFixture.hpp
 * @brief Common test data and test fixture for parameterised tests on a named test fluid
 */

#ifndef GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUIDFIXTURE_HPP_
#define GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUIDFIXTURE_HPP_

#include "TestFluids.hpp"
#include "TestFluidUtilities.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"

#include <gtest/gtest.h>

namespace geos
{

namespace testing
{

/**
 * @brief Common data for a single test point: the state at which a fluid property is evaluated
 * @tparam NC the number of components
 * @details This is an aggregate. Test-specific parameter types extend it by inheritance with their own
 *          inputs and expected values, and remain aggregates. They are initialised with the base fields
 *          in their own braces, e.g. for a type with two extra fields:
 * @code
 *   DensityData< 2 >{ { pressure, temperature, { z0, z1 } }, molarDensity, massDensity }
 * @endcode
 *          The inner braces around the base are required to avoid -Wmissing-braces with clang.
 */
template< integer NC >
struct TestDataBase
{
  real64 pressure;         ///< Pressure (Pa)
  real64 temperature;       ///< Temperature (K)
  Feed< NC > composition;   ///< Composition (mole fractions)
};

/**
 * @brief Writes the fields of a test data point as a comma separated list of numbers
 * @details Used to implement the stream operators of @c TestDataBase and the types derived from it,
 *          which gtest uses to print the parameters of a test. Each kind of value has a fixed format:
 *          pressure "%.2e", temperature "%.2f", composition "%.6f" and outputs "%.7e". Booleans are
 *          padded so that the fields that follow them line up in a column. Arrays are printed in braces.
 * @code
 *   template< integer NC >
 *   std::ostream & operator<<( std::ostream & os, DensityData< NC > const & data )
 *   {
 *     return TestDataPrinter( os, data ).output( data.expectedMolarDensity ).output( data.expectedMassDensity );
 *   }
 * @endcode
 */
class TestDataPrinter
{
public:
  /// Start printing with the common fields: pressure, temperature and composition
  template< integer NC >
  TestDataPrinter( std::ostream & os, TestDataBase< NC > const & data ):
    m_os( os )
  {
    pressure( data.pressure );
    temperature( data.temperature );
    composition( data.composition );
  }

  /// Print a pressure
  TestDataPrinter & pressure( real64 const value )
  {
    return field( GEOS_FMT( "{:.2e}", value ) );
  }

  /// Print a temperature
  TestDataPrinter & temperature( real64 const value )
  {
    return field( GEOS_FMT( "{:.2f}", value ) );
  }

  /**
   * @brief Print a composition (mole fractions)
   * @note The size of a @c std::array is a @c std::size_t, which cannot be deduced as an @c integer,
   *       so the array type is the parameter and its size is taken from @c std::tuple_size.
   *       This also restricts the template to arrays.
   */
  template< typename ARRAY, integer N = static_cast< integer >( std::tuple_size< ARRAY >::value ) >
  TestDataPrinter & composition( ARRAY const & values )
  {
    return field( formatArray< N >( values, []( real64 const v ){ return GEOS_FMT( "{:.6f}", v ); } ) );
  }

  /// Print an output (expected) value
  TestDataPrinter & output( real64 const value )
  {
    return field( formatOutput( value ) );
  }

  /// Print an array of output (expected) values
  template< typename ARRAY, integer N = static_cast< integer >( std::tuple_size< ARRAY >::value ) >
  TestDataPrinter & output( ARRAY const & values )
  {
    return field( formatArray< N >( values, &formatOutput ) );
  }

  /// Print an integer value, such as a phase index, an option or a status
  TestDataPrinter & value( integer const number )
  {
    return field( GEOS_FMT( "{}", number ) );
  }

  /// Print a boolean; the following field is aligned for both values ("true,  " and "false, ")
  TestDataPrinter & flag( bool const value )
  {
    field( value ? "true" : "false" );
    m_padding = value ? 1 : 0;
    return *this;
  }

  /// The stream being written to
  operator std::ostream &() const
  {
    return m_os;
  }

private:
  static string formatOutput( real64 const value )
  {
    return GEOS_FMT( "{:.7e}", value );
  }

  template< integer N, typename ARRAY, typename FORMAT >
  static string formatArray( ARRAY const & values, FORMAT && format )
  {
    string text = "{";
    for( integer i = 0; i < N; ++i )
    {
      text += ( 0 < i ? ", " : "" ) + format( values[i] );
    }
    return text + "}";
  }

  TestDataPrinter & field( string const & text )
  {
    if( m_started )
    {
      m_os << "," << string( 1 + m_padding, ' ' );
    }
    m_os << text;
    m_started = true;
    m_padding = 0;
    return *this;
  }

  std::ostream & m_os;
  bool m_started{ false };
  integer m_padding{ 0 };
};

/// Print the common fields of a test data point
template< integer NC >
std::ostream & operator<<( std::ostream & os, TestDataBase< NC > const & data )
{
  return TestDataPrinter( os, data );
}

/**
 * @brief Step sizes, scaling and tolerances for the numerical derivative checks
 * @details Default construct and change the fields that differ, e.g.
 * @code
 *   DerivativeCheckSettings settings;
 *   settings.compositionStep = 1.0e-6;
 * @endcode
 */
struct DerivativeCheckSettings
{
  real64 pressureStep{ 1.0e-4 };            ///< Pressure step, relative to the pressure
  real64 temperatureStep{ 1.0e-6 };         ///< Temperature step, relative to the temperature
  real64 compositionStep{ 1.0e-7 };         ///< Absolute composition step
  real64 pressureScale{ 1.0 };              ///< Scale applied to values and derivatives when checking pressure derivatives
  real64 temperatureScale{ 1.0 };           ///< Scale applied to values and derivatives when checking temperature derivatives
  real64 compositionScale{ 1.0 };           ///< Scale applied to values and derivatives when checking composition derivatives
  bool checkComposition{ true };            ///< Whether composition derivatives are checked
  real64 absTolerance{ internal::absTol };  ///< Absolute tolerance
  real64 relTolerance{ internal::relTol };  ///< Relative tolerance
};

/**
 * @brief Base fixture for parameterised tests on a named test fluid
 * @tparam FLUID the name of the test fluid
 * @tparam DATA the test parameter type, templated on the number of components
 *
 * @note Because this is a dependent base class, derived class templates must bring the members
 *       they use into scope, e.g. with @c using Base::numComps;
 */
template< TestFluidName FLUID, template< integer > class DATA >
class TestFluidFixture : public ::testing::TestWithParam< DATA< getNumberOfComponents< FLUID > > >
{
public:
  static constexpr integer numComps = getNumberOfComponents< FLUID >;
  static constexpr integer numDofs = numComps + 2;
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  using Deriv = constitutive::multifluid::DerivativeOffset;
  using DataType = DATA< numComps >;

  /**
   * @brief Check analytical derivatives with respect to pressure, temperature and every component
   *        against finite differences
   * @tparam NUM_VALUES the number of values whose derivatives are checked
   * @param data the state about which the derivatives are checked
   * @param derivatives the analytical derivatives, one row of size numDofs per value
   * @param evaluate a function with signature
   *        ( real64 pressure, real64 temperature, arraySlice1d< real64 const > composition, values )
   *        which evaluates the NUM_VALUES values at a perturbed state
   * @param settings the step sizes, scaling and tolerances
   * @param scales scaling factors applied to each value and its derivatives before comparison
   */
  template< integer NUM_VALUES, typename EVALUATE >
  static void checkNumericalDerivatives( TestDataBase< numComps > const & data,
                                         arraySlice2d< real64 const > const & derivatives,
                                         EVALUATE && evaluate,
                                         DerivativeCheckSettings const & settings = {},
                                         std::array< real64, NUM_VALUES > const & scales = unitScales< NUM_VALUES >() )
  {
    stackArray1d< real64, numComps > composition( numComps );
    for( integer ic = 0; ic < numComps; ++ic )
    {
      composition[ic] = data.composition[ic];
    }

    integer const numCheckedDofs = settings.checkComposition ? numDofs : Deriv::dC;
    stackArray1d< real64, NUM_VALUES > analytical( NUM_VALUES );
    for( integer idof = 0; idof < numCheckedDofs; ++idof )
    {
      real64 const dofScale = idof == Deriv::dP ? settings.pressureScale
                            : idof == Deriv::dT ? settings.temperatureScale
                            : settings.compositionScale;
      for( integer i = 0; i < NUM_VALUES; ++i )
      {
        analytical[i] = dofScale * scales[i] * derivatives( i, idof );
      }

      integer const ic = idof - Deriv::dC;
      real64 const step = idof == Deriv::dP ? settings.pressureStep * data.pressure
                        : idof == Deriv::dT ? settings.temperatureStep * data.temperature
                        : settings.compositionStep;

      internal::testNumericalDerivative< NUM_VALUES >( 0.0, step, analytical.toSliceConst(),
                                                       [&]( real64 const delta, auto & values )
      {
        real64 const pressure = data.pressure + ( idof == Deriv::dP ? delta : 0.0 );
        real64 const temperature = data.temperature + ( idof == Deriv::dT ? delta : 0.0 );
        if( 0 <= ic )
        {
          composition[ic] += delta;
        }
        evaluate( pressure, temperature, composition.toSliceConst(), values );
        if( 0 <= ic )
        {
          composition[ic] = data.composition[ic];
        }
        for( integer i = 0; i < NUM_VALUES; ++i )
        {
          values[i] *= dofScale * scales[i];
        }
      }, settings.absTolerance, settings.relTolerance );
    }
  }

  /// Unit scaling factors
  template< integer NUM_VALUES >
  static constexpr std::array< real64, NUM_VALUES > unitScales()
  {
    std::array< real64, NUM_VALUES > scales{};
    for( auto & scale : scales )
    {
      scale = 1.0;
    }
    return scales;
  }

protected:
  TestFluidFixture():
    m_fluid( createTestFluid< FLUID >() )
  {}

  std::unique_ptr< TestFluid< numComps > > m_fluid{};
};

} // namespace testing

} // namespace geos

#endif //GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUIDFIXTURE_HPP_
