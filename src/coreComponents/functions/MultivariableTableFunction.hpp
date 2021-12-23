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

  /// maximum dimensions for the coordinates in the table
  static constexpr integer maxDimensions = 5;

  /**
   * @brief The constructor
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  MultivariableTableFunction( const string & name,
                              dataRepository::Group * const parent );

  /**
   * @brief The catalog name interface
   * @return name of the MultivariableTableFunction in the FunctionBase catalog
   */
  static string catalogName() { return "MultivariableTableFunction"; }

  /**
   * @brief Initialize the table function
   */
  virtual void initializeFunction() override;

  /**
   * @brief Build the maps used to evaluate the table function
   */
  void reInitializeFunction();

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
  ArrayOfArrays< real64 > & getCoordinates() { return m_coordinates; }

  /**
   * @brief Get the Values object
   * 
   * @return array1d< real64 >& a reference to an array where table values are stored in fortran order
   */
  array1d< real64 > & getValues() { return m_values; }

  /**
   * @brief Set the table coordinates
   * @param coordinates An array of arrays containing table coordinate definitions
   */
  void setTableCoordinates( array1d< real64_array > const & coordinates );

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

  /// An array of table axes
  ArrayOfArrays< real64 > m_coordinates;

  /// Table values (in fortran order)
  array1d< real64 > m_values;

};


} /* namespace geosx */

#endif /* GEOSX_FUNCTIONS_MULTIVARIABLETABLEFUNCTION_HPP_ */
