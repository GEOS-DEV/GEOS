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
 * @file FluidModelTest.hpp
 */

#ifndef GEOS_CORECOMPONENTS_CONSTITUTIVE_UNITTESTS_FLUIDMODELTEST_HPP_
#define GEOS_CORECOMPONENTS_CONSTITUTIVE_UNITTESTS_FLUIDMODELTEST_HPP_

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "TestFluid.hpp"

#include <conduit.hpp>
#include <gtest/gtest.h>

namespace geos
{

class TableFunction;

namespace testing
{

/**
 * @brief A generic test fixture for a multiphase fluid
 * @details Contains the tools and features to create a fluid model and test it with particular inputs.
 * @tparam FLUID_TYPE - the type of fluid model to be tested. This should be a class derived from @c MultiFluidBase
 * @tparam NUM_COMP - the number of components in the fluid model
 * @tparam NUM_PHASE - the number of phases in the fluid model
 */
template< typename FLUID_TYPE, integer NUM_COMP, integer NUM_PHASE = 2 >
class FluidModelTest : public ::testing::Test
{
public:
  using FluidModel = FLUID_TYPE;
  using FluidWrapper = typename FLUID_TYPE::KernelWrapper;
  static constexpr integer numComp = NUM_COMP;
  static constexpr integer numDof = numComp + 2;
  static constexpr integer numPhase = NUM_PHASE;
  static constexpr real64 relTolerance = 1.0e-5;
  static constexpr real64 absTolerance = 1.0e-4;

public:
  /**
   * A point at which the fluid model will be tests. The inputs to the fluid model.
   */
  using TestPoint = std::tuple<
    real64 const,           // pressure
    real64 const,           // temperature
    Feed< NUM_COMP > const  // composition
    >;

  /**
   * The result of a fluid calculation at a point
   */
  using TestResult = std::tuple<
    Feed< NUM_PHASE > const,  // phase fraction
    Feed< NUM_PHASE > const,  // phase density,
    Feed< NUM_PHASE > const,  // phase mass density
    Feed< NUM_PHASE > const,  // phase viscosity
    Feed< NUM_PHASE > const,  // phase enthalpy
    Feed< NUM_PHASE > const,  // phase internal energy
    real64 const              // total density,
    >;

public:
  FluidModelTest();
  ~FluidModelTest() override = default;

  /**
   * @brief Retrieve a fluid model
   * @param name - the name of the fluid model to be retrieved
   * @return A pointer to the named fluid model
   */
  FluidModel * getFluid( string const & name );

  /**
   * @brief Create a fluid model
   * @tparam LAMBDA - a function that takes a reference to the
   * @param name - a name for the fluid model to be created
   * @param function - a function to process the fluid model. This will be of the type @c std::function<void(FluidModel&)>.
   *         It will be called to populate the parameters of the model to make sure it has a valid definition.
   * @return A pointer to the newly created fluid model
   */
  template< typename LAMBDA >
  FluidModel * createFluid( string const & name, LAMBDA && function );

  /**
   * @brief Test implementation against known (expected) valuse
   * @details Will test the values returned by the fluid compute against values that are known or expected
   * @param fluid - a pointer to the fluid model. Must be derived from @c MultiFluidBase
   * @param testPoint - the test point input data (pressure, temperature and composition)
   * @param expectedValues - the expected values
   * @param relTol - the relative tolerance to use in the check
   * @param absTol - the absolute tolerance to use in the check
   */
  void testValuesAgainstPreviousImplementation( FluidModel * fluid,
                                                TestPoint const & testPoint,
                                                TestResult const & expectedValues,
                                                real64 const relTol = relTolerance,
                                                real64 const absTol = absTolerance );

  /**
   * @brief Tests the derivatives of the fluid model
   * @details Will test the derivatives returned by the fluid model at the test point specified. The input
   *          values will be perturbed to calculate left, central and right numerical derivatives at the
   *          selected point. These will then be compared with the calculated analytical derivative. The if
   *          the smallest difference between the analytical derivative and the 3 numerical derivatives is
   *          larger than the tolerance, the test will fail.
   * @param fluid - a pointer to the fluid model. Must be derived from @c MultiFluidBase
   * @param data - the test point input data (pressure, temperature and composition)
   * @param perturbationLevel - fraction of parameter to use to perturb when calculating numerical derivatives
   * @param relTol - the relative tolerance to use in the check
   * @param absTol - the absolute tolerance to use in the check
   */
  void testNumericalDerivatives( FluidModel * fluid,
                                 TestPoint const & data,
                                 real64 const perturbationLevel = 1.0e-4,
                                 real64 const relTol = relTolerance,
                                 real64 const absTol = absTolerance );

protected:
  /**
   * @brief Writes content to a file
   * @param fileName - the name of the file. Will be overwritten in it already exists
   * @param content - the content to be written to the file
   */
  static void writeTableToFile( string const & fileName, char const * content );

  /**
   * @brief Removes a file
   * @param fileName - the name of the file to remove
   */
  static void removeFile( string const & fileName );

  /**
   * @brief Creates a table function
   * @details Creates a linearly interpolated 1D table function
   * @param tableName - the name of the table function
   * @param coordinates - an array of coordinates for the table
   * @param values - an array of values of the function at the corresponding coordinates
   * @return a pointer to the newly created table function
   */
  static TableFunction const * createTable( string const & tableName, real64_array const & coordinates, real64_array const & values );

  /**
   * @brief Creates a table function
   * @details Creates a linearly interpolated potentially multi-dimensional table function
   * @param tableName - the name of the table function
   * @param coordinates - an array of coordinates for the table
   * @param values - an array of values of the function at the corresponding coordinates
   * @return a pointer to the newly created table function
   */
  static TableFunction const * createTable( string const & tableName, array1d< real64_array > const & coordinates, real64_array const & values );

  /**
   * @brief Fill in an array
   * @details Will fill an array possibly extending the array with values from a list
   * @tparam ARRAY the type of output array. Should support @c emplace_back
   * @tparam LIST the type of list for the data
   * @param array - the output array. This should ideally be unallocated and will be extended with values from the input list.
   * @param data - the input array.
   */
  template< typename ARRAY, typename LIST >
  static void fill( ARRAY & array, LIST const & data );

  /**
   * @brief Populate an array
   * @details Will populate an array with values from a list. The array should be pre-allocated prioir to calling this method.
   * @tparam ARRAY the type of output array.
   * @tparam LIST the type of list for the data
   * @param array - the output array. The values in this will be replaced by values from the input list.
   * @param data - the input array.
   */
  template< typename ARRAY, typename LIST >
  static void populate( ARRAY & array, LIST const & data );

  /**
   * @brief Convert an array slice into a string
   * @details Will print out an array slice into a string
   * @tparam USD The leading dimension stride
   * @param array - The array slice.
   * @return A string printout of the array
   */
  template< integer USD >
  static string toString( arraySlice1d< real64 const, USD > const & array );

  /**
   * @brief Populate an array with linearly spaced data
   * @details Will allocated and populate the array with linearly spaced values in the interval [x0,x1] including the end
   *          points using @c n intervals. Note that @n is the number of intervals so the list will be of length @c n+1 on exit.
   * @param array - the array to be populated. This should not be allocated
   * @param x0 - the start value for the list
   * @param x1 - the end value for the list
   * @param n - the number of intervals (array will have n+1 points)
   */
  static void populateLinearScale( array1d< real64 > & array, real64 const x0, real64 const x1, integer const n );

  /**
   * @brief Populate an array with logarithmically spaced data
   * @details Will allocated and populate the array with logarithmically spaced values in the interval [x0,x1] including the end
   *          points using @c n intervals. Note that @n is the number of intervals so the list will be of length @c n+1 on exit.
   * @param array - the array to be populated. This should not be allocated
   * @param x0 - the start value for the list (must be positive)
   * @param x1 - the end value for the list (must be positive)
   * @param n - the number of intervals (array will have n+1 points)
   */
  static void populateLogScale( array1d< real64 > & array, real64 const x0, real64 const x1, integer const n );

  /**
   * @brief Checks the derivatives of a property against numerical values
   * @details Tests the analytical derivatives calculated for a property in an @c NDIM dimensonal array. The
   *          values are provided in @c valueArray which includes all the values calculated for pertubations of
   *          the input primary variables. The analytical derivatives are provided in the @c derivArray array.
   * @tparam NDIM - the dimension of the value array
   * @param propName - the name of the property whose derivatives are being tests
   * @param testValues - a description of the test point (pressure, temperature and composition)
   * @param valueArray - array of calculated values including values at perturbed
   * @param derivArray - array of calculated analytical derivatives
   * @param displacements - values of displacements at each of the test points
   * @param valueScale - a scale to use on the values before calculating the finite differences to reduce
   *                     round off error.
   * @param dofNames - names for the degrees of freedom variables (for debug)
   * @param relTol - the relative tolerance to use in the check
   * @param absTol - the absolute tolerance to use in the check
   * @param indices - An indicator of which value in the array is being tested. This are the trailing 2
   *            indices in @c valueArray
   */
  template< integer NDIM, typename ... INDICES, integer USD1, integer USD2, integer USD3,
            typename=std::enable_if_t< sizeof ... ( INDICES ) == NDIM-2 > >
  static void testDerivatives( string const propName,
                               string const testValues,
                               ArrayView< real64 const, NDIM, USD1 > const & valueArray,
                               ArrayView< real64 const, NDIM+1, USD2 > const & derivArray,
                               ArraySlice< real64 const, 1, USD3 > const & displacements,
                               real64 const valueScale,
                               string_array const & dofNames,
                               real64 const relTol,
                               real64 const absTol,
                               INDICES const ... indices );

private:
  void createFunctionManager();

protected:
  conduit::Node m_node;
  dataRepository::Group m_parent;
};

} // namespace testing
} // namespace geos

// Implementation details in separate file
#include "FluidModelTest_impl.hpp"

#endif // GEOS_CORECOMPONENTS_CONSTITUTIVE_UNITTESTS_FLUIDMODELTEST_HPP_
