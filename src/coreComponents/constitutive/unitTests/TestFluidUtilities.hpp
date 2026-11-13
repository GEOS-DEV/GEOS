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
 * @file TestFluidUtilities.hpp
 */

#ifndef GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUIDUTILITIES_HPP_
#define GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUIDUTILITIES_HPP_

#include "codingUtilities/UnitTestUtilities.hpp"
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "common/format/table/TableLayout.hpp"

#include <algorithm>
#include <initializer_list>

namespace geos
{

namespace testing
{

namespace internal
{

static constexpr real64 relTol = 1.0e-4;
static constexpr real64 absTol = 1.0e-7;

/// @return the scale used for the relative error between two derivatives (as in @c checkRelativeError)
inline real64 derivativeErrorScale( real64 const analytical, real64 const numerical )
{
  return std::max( LvArray::math::abs( analytical ), LvArray::math::abs( numerical ) ) + 1.0;
}

/// @return the relative error between an analytical and a numerical derivative
inline real64 derivativeRelativeError( real64 const analytical, real64 const numerical )
{
  return LvArray::math::abs( analytical - numerical ) / derivativeErrorScale( analytical, numerical );
}

/// @return true if a numerical derivative matches the analytical derivative (same criterion as @c checkRelativeError)
inline bool derivativeWithinTolerance( real64 const analytical, real64 const numerical,
                                       real64 const relTolerance, real64 const absTolerance )
{
  real64 const delta = LvArray::math::abs( analytical - numerical );
  return delta <= absTolerance || delta <= relTolerance * derivativeErrorScale( analytical, numerical );
}

/**
 * @brief Add a row to a derivative check table
 * @param table the table
 * @param index the index of the value
 * @param values the function values and derivatives to show
 * @param error the relative error
 * @param passed whether the check passed
 * @return 1 if the check failed, 0 otherwise
 */
inline integer addDerivativeRow( TableData & table,
                                 integer const index,
                                 std::initializer_list< real64 > const values,
                                 real64 const error,
                                 bool const passed )
{
  stdVector< TableData::CellData > row;
  row.push_back( { CellType::Value, GEOS_FMT( "{}", index ) } );
  for( real64 const value : values )
  {
    row.push_back( { CellType::Value, GEOS_FMT( "{:.5e}", value ) } );
  }
  row.push_back( { CellType::Value, GEOS_FMT( "{:.5e}", error ) } );
  row.push_back( { CellType::Value, passed ? "ok" : "FAIL" } );
  table.addRow( row );
  return passed ? 0 : 1;
}

/**
 * @brief Issue a single test failure showing the derivative check table if any check failed
 * @param valueColumns the names of the value columns
 * @param table the table of results
 * @param numFailures the number of failed checks
 * @param numValues the number of values checked
 * @param x the value at which the derivatives are evaluated
 * @param dx the perturbation used for the numerical derivatives
 * @param relTolerance the relative tolerance
 * @param absTolerance the absolute tolerance
 */
inline void reportDerivativeFailures( std::initializer_list< string > const valueColumns,
                                      TableData const & table,
                                      integer const numFailures,
                                      integer const numValues,
                                      real64 const x,
                                      real64 const dx,
                                      real64 const relTolerance,
                                      real64 const absTolerance )
{
  if( numFailures == 0 )
  {
    return;
  }
  stdVector< string > columns{ "index" };
  columns.insert( columns.end(), valueColumns.begin(), valueColumns.end() );
  columns.emplace_back( "rel. error" );
  columns.emplace_back( "status" );
  TableLayout const layout( "", columns );
  TableTextFormatter const formatter( layout );
  ADD_FAILURE() << GEOS_FMT( "Numerical derivative check failed for {} of {} value(s) at x = {:.5e} with dx = {:.5e} "
                             "(relTol = {:.5e}, absTol = {:.5e})\n",
                             numFailures, numValues, x, dx, relTolerance, absTolerance )
                << formatter.toString( table );
}

/**
 * @brief Tests a multi-valued function against a derivative
 * @details Will calculate the left-sided and the right-sided numerical derivatives of a function
 *          and compare this against a analytically calculated values provided.
 * @tparam numValues the number of values that the function returns
 * @tparam FUNCTION the type of function (typically a lambda)
 * @param x The value at which the function should be evaluated
 * @param dx The value to use to perturb @c x in the calculation of the numerical derivatives
 * @param derivatives The values of the analytically calculated derivatives to use for comparison
 * @param function The function which is being tested. This should be a function that takes 2 parameters.
 *       The first is the value at which the function is being evaluated (x) and the second is an array
 *       of size @c numValues which is the result of the execution of the function.
 * @param absTolerance The absolute tolerance to use for the comparison
 * @param relTolerance The relative tolerance to use for the comparison
 */
template< integer numValues, typename FUNCTION >
void testNumericalDerivative( real64 const x,
                              real64 const dx,
                              arraySlice1d< real64 const > const & derivatives,
                              FUNCTION && function,
                              real64 const absTolerance = absTol,
                              real64 const relTolerance = relTol )
{
  stackArray1d< real64, numValues > leftValues( numValues );
  stackArray1d< real64, numValues > centreValues( numValues );
  stackArray1d< real64, numValues > rightValues( numValues );
  function( x-dx, leftValues );
  function( x, centreValues );
  function( x+dx, rightValues );

  TableData table;
  integer numFailures = 0;
  for( integer i = 0; i < numValues; ++i )
  {
    real64 const leftDerivative = (centreValues[i] - leftValues[i]) / dx;
    real64 const rightDerivative = (rightValues[i] - centreValues[i]) / dx;
    real64 const centralDerivative = 0.5*(rightValues[i] - leftValues[i]) / dx;

    // Choose from the left, central and right derivatives, the one that's nearest the analytical value
    real64 selectedDerivative = leftDerivative;
    for( real64 const deriv : {rightDerivative, centralDerivative} )
    {
      if( LvArray::math::abs( deriv - derivatives[i] ) < LvArray::math::abs( selectedDerivative - derivatives[i] ) )
      {
        selectedDerivative = deriv;
      }
    }
    numFailures += addDerivativeRow( table, i,
                                     {leftValues[i], centreValues[i], rightValues[i], derivatives[i],
                                      leftDerivative, rightDerivative, centralDerivative},
                                     derivativeRelativeError( derivatives[i], selectedDerivative ),
                                     derivativeWithinTolerance( derivatives[i], selectedDerivative, relTolerance, absTolerance ) );
  }
  reportDerivativeFailures( {"f(x-dx)", "f(x)", "f(x+dx)", "analytical", "left", "right", "central"},
                            table, numFailures, numValues, x, dx, relTolerance, absTolerance );
}

/**
 * @brief Tests a function against a derivative
 * @details Single-valued version of the multi-valued check below: the test passes if the left, right
 *          or central numerical derivative matches the analytical value.
 * @param x The value at which the function should be evaluated
 * @param dx The value to use to perturb @c x in the calculation of the numerical derivatives
 * @param derivative The value of the analytically calculated derivative to use for comparison
 * @param function The function which is being tested. This should be a function that returns
 *        a single real value.
 * @param absTolerance The absolute tolerance to use for the comparison
 * @param relTolerance The relative tolerance to use for the comparison
 */
template< typename FUNCTION >
void testNumericalDerivative( real64 const x, real64 const dx, real64 const derivative,
                              FUNCTION && function,
                              real64 const absTolerance = absTol, real64 const relTolerance = relTol )
{
  stackArray1d< real64, 1 > derivatives( 1 );
  derivatives[0] = derivative;
  testNumericalDerivative< 1 >( x, dx, derivatives.toSliceConst(),
                                [&]( real64 const xi, auto & values ){ values[0] = function( xi ); },
                                absTolerance, relTolerance );
}

/**
 * @brief Tests a multi-valued function against a second derivative
 * @details Will calculate the left-sided and the right-sided numerical second derivatives of a function
 *          and compare this against a analytically calculated values provided.
 * @tparam numValues the number of values that the function returns
 * @tparam FUNCTION the type of function (typically a lambda)
 * @param x The value at which the function should be evaluated
 * @param dx The value to use to perturb @c x in the calculation of the numerical derivatives
 * @param derivatives The values of the analytically calculated derivatives to use for comparison
 * @param function The function which is being tested. This should be a function that takes 2 parameters.
 *       The first is the value at which the function is being evaluated (x) and the second is an array
 *       of size @c numValues which is the result of the execution of the function.
 * @param absTolerance The absolute tolerance to use for the comparison
 * @param relTolerance The relative tolerance to use for the comparison
 */
template< integer numValues, typename FUNCTION >
void testNumericalSecondDerivative( real64 const x,
                                    real64 const dx,
                                    arraySlice1d< real64 const > const & derivatives,
                                    FUNCTION && function,
                                    real64 const absTolerance = absTol,
                                    real64 const relTolerance = relTol )
{
  stackArray2d< real64, 5*numValues > values( 5, numValues );
  function( x-2.0*dx, values[0] );
  function( x-dx, values[1] );
  function( x, values[2] );
  function( x+dx, values[3] );
  function( x+2.0*dx, values[4] );

  real64 constexpr stencils[6][5] = {
    {1.0, -2.0, 1.0, 0.0, 0.0},
    {-1.0, 4.0, -5.0, 2.0, 0.0},
    {-1.0/12.0, 16.0/12.0, -30.0/12.0, 16.0/12.0, -1.0/12.0},
    {0.0, 1.0, -2.0, 1.0, 0.0},
    {0.0, -1.0, 4.0, -5.0, 2.0},
    {0.0, 0.0, 1.0, -2.0, 1.0},
  };
  real64 const invdx2 = 1.0 / (dx*dx);
  TableData table;
  integer numFailures = 0;
  for( integer i = 0; i < numValues; ++i )
  {
    // Choose the stencil whose second derivative is nearest the analytical value
    real64 minError = LvArray::NumericLimits< real64 >::max;
    real64 selectedDerivative = 0.0;
    for( integer si = 0; si < 6; si++ )
    {
      real64 deriv = 0.0;
      for( integer ci = 0; ci < 5; ci++ )
      {
        deriv += stencils[si][ci]*values( ci, i );
      }
      deriv *= invdx2;
      real64 const error = LvArray::math::abs( deriv - derivatives[i] );
      if( error < minError )
      {
        minError = error;
        selectedDerivative = deriv;
      }
    }
    numFailures += addDerivativeRow( table, i,
                                     {values( 0, i ), values( 1, i ), values( 2, i ), values( 3, i ), values( 4, i ),
                                      derivatives[i], selectedDerivative},
                                     derivativeRelativeError( derivatives[i], selectedDerivative ),
                                     derivativeWithinTolerance( derivatives[i], selectedDerivative, relTolerance, absTolerance ) );
  }
  reportDerivativeFailures( {"f(x-2dx)", "f(x-dx)", "f(x)", "f(x+dx)", "f(x+2dx)", "analytical", "numerical"},
                            table, numFailures, numValues, x, dx, relTolerance, absTolerance );
}

}// namespace internal

}// namespace testing

}// namespace geos

#endif //GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUIDUTILITIES_HPP_
