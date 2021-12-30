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

/**
 * @file MultivariableTableFunction.hpp
 */

#ifndef GEOSX_FUNCTIONS_MULTIVARIABLETABLEFUNCTION_HPP_
#define GEOSX_FUNCTIONS_MULTIVARIABLETABLEFUNCTION_HPP_

#include "FunctionBase.hpp"

#include "codingUtilities/EnumStrings.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{

/**
 * @class MultivariableTableFunction
 *
 * An interface class for multivariable table function (function with multiple inputs and outputs)
 * Prepares input data for MultivariableStaticInterpolatorKernel, which performes actual interpolation
 */

class MultivariableTableFunction : public FunctionBase
{
public:

  /**
   * @brief The constructor
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  MultivariableTableFunction( const string & name,
                              Group * const parent );

  /**
   * @brief The catalog name interface
   * @return name of the MultivariableTableFunction in the FunctionBase catalog
   */
  static string catalogName() { return "MultivariableTableFunction"; }

  /**
   * @brief Initialize the table function
   */
  virtual void initializeFunction() override;

  void evaluate( dataRepository::Group const & group,
                 real64 const time,
                 SortedArrayView< localIndex const > const & set,
                 real64_array & result ) const
  {};

  /**
   * @brief Method to evaluate a function
   * @param input a scalar input
   * @return the function evaluation
   */
  virtual real64 evaluate( real64 const * const input ) const
  {
    return 0;
  };

  /**
   * @brief Get the table axes definitions
   * @return a reference to an array of arrays that define each table axis
   */
  real64_array & getAxisMinimums() { return m_axisMinimums; }

  /**
   * @brief Get the table axes definitions
   * @return a reference to an array of arrays that define each table axis
   */
  real64_array & getAxisMaximums() { return m_axisMaximums; }

  /**
   * @brief Get the table axes definitions
   * @return a reference to an array of arrays that define each table axis
   */
  integer_array & getAxisPoints() { return m_axisPoints; }


  /**
   * @brief Get the table axes definitions
   * @return a reference to an array of arrays that define each table axis
   */
  real64_array & getAxisSteps() { return m_axisSteps; }

  /**
   * @brief Get the table axes definitions
   * @return a reference to an array of arrays that define each table axis
   */
  real64_array & getAxisStepInvs() { return m_axisStepInvs; }


  globalIndex_array & getAxisHypercubeMults() { return m_axisHypercubeMults; }
  /**
   * @brief Get the Values object
   *
   * @return real64_array& a reference to an array where table values are stored in fortran order
   */
  real64_array & getHypercubeData() { return m_hypercubeData; }

  /**
   * @brief Set the table coordinates
   */
  void setTableCoordinates( integer const numDims,
                            integer const numOps,
                            real64_array const & axisMinimums,
                            real64_array const & m_axisMaximums,
                            integer_array const & m_axisPoints );

  /**
   * @brief Set the table values
   * @param values An array of table values in fortran order
   */
  void setTableValues( real64_array values );

private:

  /**
   * @brief Parse a table file.
   * @tparam T The type for table or axis values.
   * @param[in] target The place to store values.
   * @param[in] filename The name of the file to read.
   * @param[in] delimiter The delimiter used for file entries.
   */
  template< typename T >
  void parseFile( string const & filename, array1d< T > & target );


  void getHypercubePoints( globalIndex hypercubeIndex, globalIndex_array & hypercubePoints );

  /// Number of table dimensions
  integer m_numDims;

  /// Number of operators - functions to interpolate
  integer m_numOps;

  /// Number of vertices in each hypercube
  integer m_numVerts;

  /// Array [numDims] of axis minimum values
  real64_array m_axisMinimums;

  /// Array [numDims] of axis maximum values
  real64_array m_axisMaximums;

  /// Array [numDims] of axis discretization points
  integer_array m_axisPoints;

  // inputs : service data derived from table discretization data

  ///  Array [numDims] of axis interval lengths (axes are discretized uniformly)
  real64_array m_axisSteps;

  ///  Array [numDims] of inversions of axis interval lengths (axes are discretized uniformly)
  real64_array m_axisStepInvs;

  ///  Array [numDims] of point index mult factors for each axis
  globalIndex_array m_axisPointMults;

  ///  Array [numDims] of hypercube index mult factors for each axis
  globalIndex_array m_axisHypercubeMults;

  // inputs: operator sample data

  ///  Main table data stored per point: all operator values for each table coordinate
  real64_array m_pointData;

  ///  Main table data stored per hypercube: all values required for interpolation withing give hypercube are stored contiguously
  real64_array m_hypercubeData;

};


} /* namespace geosx */

#endif /* GEOSX_FUNCTIONS_MULTIVARIABLETABLEFUNCTION_HPP_ */
