/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TableFunction.hpp
 */

#ifndef GEOSX_MANAGERS_FUNCTIONS_TABLEFUNCTION_HPP_
#define GEOSX_MANAGERS_FUNCTIONS_TABLEFUNCTION_HPP_

#include "codingUtilities/EnumStrings.hpp"
#include "managers/Functions/FunctionBase.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{

/**
 * @class TableFunction
 *
 * An interface for a dense table-based function
 */
class TableFunction : public FunctionBase
{
public:

  /// Enumerator of available interpolation types
  enum class InterpolationType : integer
  {
    Linear,
    Nearest,
    Upper,
    Lower
  };

  /// maximum dimensions for the coordinates in the table
  static constexpr localIndex maxDimensions = 4;

  /// maximum number of corners
  static constexpr localIndex maxNumCorners = 16;

  /**
   * @class KernelWrapper
   *
   * A nested class encapsulating the kernel function doing the interpolation in the table
   */
  class KernelWrapper
  {
public:

    /**
     * @brief The constructor of the kernel wrapper
     * @param[in] interpolationMethod table interpolation method
     * @param[in] coordinates array of table axes
     * @param[in] values table values (in fortran order)
     * @param[in] dimensions number of active table dimensions
     * @param[in] size size of the table
     * @param[in] indexIncrement array used to locate values within ND tables
     * @param[in] corners corners of the box that surround the value in N dimensions
     * @param[in] numCorners number of active table corners
     */
    KernelWrapper( TableFunction::InterpolationType interpolationMethod,
                   ArrayOfArraysView< real64 const > const & coordinates,
                   arrayView1d< real64 const > const & values,
                   localIndex dimensions,
                   arrayView1d< localIndex const > const & size,
                   arrayView1d< localIndex const > const & indexIncrement,
                   localIndex const (&corners)[TableFunction::maxDimensions][maxNumCorners],
                   localIndex const numCorners );

    /// Default constructor for the kernel wrapper
    KernelWrapper();

    /// Default copy constructor
    KernelWrapper( KernelWrapper const & ) = default;

    /// Default move constructor
    KernelWrapper( KernelWrapper && ) = default;

    /// Copy assignment operator
    KernelWrapper & operator=( KernelWrapper const & ) = delete;

    /// Deleted move assignment operator
    KernelWrapper & operator=( KernelWrapper && ) = delete;

    /**
     * @brief The function that populates the wrapper
     * @param[in] interpolationMethod table interpolation method
     * @param[in] coordinates array of table axes
     * @param[in] values table values (in fortran order)
     * @param[in] dimensions number of active table dimensions
     * @param[in] size size of the table
     * @param[in] indexIncrement array used to locate values within ND tables
     * @param[in] corners corners of the box that surround the value in N dimensions
     * @param[in] numCorners number of active table corners
     */
    void create( TableFunction::InterpolationType interpolationMethod,
                 ArrayOfArraysView< real64 const > const & coordinates,
                 arrayView1d< real64 const > const & values,
                 localIndex dimensions,
                 arrayView1d< localIndex const > const & size,
                 arrayView1d< localIndex const > const & indexIncrement,
                 localIndex const (&corners)[TableFunction::maxDimensions][maxNumCorners],
                 localIndex const numCorners );

    /**
     * @brief Main compute function to interpolate in the table and return the derivatives
     * @param[in] input vector of input value
     * @param[out] value interpolated value
     * @param[out] derivatives vector of derivatives of interpolated value wrt the variables present in input
     */
    template< typename IN_ARRAY, typename OUT_ARRAY >
    GEOSX_HOST_DEVICE
    void
    compute( IN_ARRAY const & input, real64 & value, OUT_ARRAY && derivatives ) const;

private:

    /// Table interpolation method
    TableFunction::InterpolationType m_interpolationMethod;

    /// An array of table axes
    ArrayOfArraysView< real64 const > m_coordinates;

    /// Table values (in fortran order)
    arrayView1d< real64 const > m_values;

    /// Number of active table dimensions
    localIndex m_dimensions;

    /// Size of the table
    arrayView1d< localIndex const > m_size;

    /// Array used to locate values within ND tables
    arrayView1d< localIndex const > m_indexIncrement;

    /**
     * @brief The corners of the box that surround the value in N dimensions
     * m_corners should be of size m_maxDimensions x (2^m_maxDimensions)
     */
    localIndex m_corners[TableFunction::maxDimensions][maxNumCorners];

    /// The number of active table corners
    localIndex m_numCorners;

  };

  /**
   * @brief The constructor
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  TableFunction( const string & name,
                 dataRepository::Group * const parent );

  /**
   * @brief The destructor
   */
  virtual ~TableFunction() override;

  /**
   * @brief The catalog name interface
   * @return name of the TableFunction in the FunctionBase catalog
   */
  static string catalogName() { return "TableFunction"; }

  /**
   * @brief Parse a table file.
   *
   * @tparam T The type for table or axis values.
   * @param[in] target The place to store values.
   * @param[in] filename The name of the file to read.
   * @param[in] delimiter The delimiter used for file entries.
   */
  template< typename T >
  void parseFile( array1d< T > & target, string const & filename, char delimiter );

  /**
   * @brief Initialize the table function
   */
  virtual void initializeFunction() override;

  /**
   * @brief Build the maps used to evaluate the table function
   */
  void reInitializeFunction();

  /**
   * @brief Method to evaluate a function on a target object
   * @param group a pointer to the object holding the function arguments
   * @param time current time
   * @param set the subset of nodes to apply the function to
   * @param result an array to hold the results of the function
   */
  virtual void evaluate( dataRepository::Group const & group,
                         real64 const time,
                         SortedArrayView< localIndex const > const & set,
                         real64_array & result ) const override final
  {
    FunctionBase::evaluateT< TableFunction >( group, time, set, result );
  }

  /**
   * @brief Method to evaluate a function
   * @param input a scalar input
   * @return the function result
   */
  virtual real64 evaluate( real64 const * const input ) const override final;

  /**
   * @brief Get the table axes definitions
   * @return a reference to an array of arrays that define each table axis
   */
  ArrayOfArraysView< real64 const > getCoordinates() const { return m_coordinates.toViewConst(); }

  /**
   * @copydoc getCoordinates() const
   */
  ArrayOfArraysView< real64 > getCoordinates() { return m_coordinates.toView(); }

  /**
   * @brief Get the table values
   * @return a reference to the 1d array of table values.  For ND arrays, values are stored in Fortran order.
   */
  array1d< real64 > const & getValues() const { return m_values; }

  /**
   * @copydoc getValues() const
   */
  array1d< real64 > & getValues()       { return m_values; }

  /**
   * @brief Get the interpolation method
   * @return The interpolation method
   */
  InterpolationType getInterpolationMethod() const { return m_interpolationMethod; }

  /**
   * @brief Set the interpolation method
   * @param method The interpolation method
   */
  void setInterpolationMethod( InterpolationType const method );

  /**
   * @brief Set the table coordinates
   * @param coordinates An array of arrays containing table coordinate definitions
   */
  void setTableCoordinates( array1d< real64_array > coordinates );

  /**
   * @brief Set the table values
   * @param values An array of table values in fortran order
   */
  void setTableValues( real64_array values );

  /**
   * @brief Create an instance of the kernel wrapper
   * @return the kernel wrapper
   */
  KernelWrapper createKernelWrapper() const;

private:
  /// Coordinates for 1D table
  real64_array m_tableCoordinates1D;

  /// List of table coordinate file names
  path_array m_coordinateFiles;

  /// Table voxel file names
  Path m_voxelFile;

  /// Table interpolation method
  InterpolationType m_interpolationMethod;

  /// An array of table axes
  ArrayOfArrays< real64 > m_coordinates;

  /// Table values (in fortran order)
  real64_array m_values;

  /// Number of active table dimensions
  localIndex m_dimensions;

  /// Size of the table
  localIndex_array m_size;

  /// Array used to locate values within ND tables
  localIndex_array m_indexIncrement;

  /**
   * @brief The corners of the box that surround the value in N dimensions
   * m_corners should be of size maxDimensions x (2^maxDimensions)
   */
  localIndex m_corners[maxDimensions][maxNumCorners];

  /// The number of active table corners
  localIndex m_numCorners;

  /// Kernel wrapper to interpolate in table and return derivatives
  KernelWrapper m_kernelWrapper;

};

template< typename IN_ARRAY, typename OUT_ARRAY >
GEOSX_HOST_DEVICE
void
TableFunction::KernelWrapper::compute( IN_ARRAY const & input, real64 & value, OUT_ARRAY && derivatives ) const
{
  value = 0.0;
  for( localIndex i = 0; i < m_dimensions; ++i )
  {
    derivatives[i] = 0.0;
  }

  // Linear interpolation
  if( m_interpolationMethod == TableFunction::InterpolationType::Linear )
  {
    localIndex bounds[TableFunction::maxDimensions][2]{};
    real64 weights[TableFunction::maxDimensions][2]{};
    real64 dWeights_dInput[TableFunction::maxDimensions][2]{};
    real64 dCornerValue_dInput[TableFunction::maxDimensions]{};

    // Determine position, weights
    for( localIndex ii=0; ii<m_dimensions; ++ii )
    {
      if( input[ii] <= m_coordinates[ii][0] )
      {
        // Coordinate is to the left of this axis
        bounds[ii][0] = 0;
        bounds[ii][1] = 0;
        weights[ii][0] = 0;
        weights[ii][1] = 1;
        dWeights_dInput[ii][0] = 0;
        dWeights_dInput[ii][1] = 0;
      }
      else if( input[ii] >= m_coordinates[ii][m_size[ii] - 1] )
      {
        // Coordinate is to the right of this axis
        bounds[ii][0] = m_size[ii] - 1;
        bounds[ii][1] = bounds[ii][0];
        weights[ii][0] = 1;
        weights[ii][1] = 0;
        dWeights_dInput[ii][0] = 0;
        dWeights_dInput[ii][1] = 0;
      }
      else
      {
        // Find the coordinate index
        ///TODO make this fast
        // Note: find uses a binary search...  If we assume coordinates are
        // evenly spaced, we can speed things up considerably
        auto lower = LvArray::sortedArrayManipulation::find( m_coordinates[ii].begin(), m_coordinates.sizeOfArray( ii ), input[ii] );
        bounds[ii][1] = LvArray::integerConversion< localIndex >( lower );
        bounds[ii][0] = bounds[ii][1] - 1;

        real64 dx = m_coordinates[ii][bounds[ii][1]] - m_coordinates[ii][bounds[ii][0]];
        weights[ii][0] = 1.0 - (input[ii] - m_coordinates[ii][bounds[ii][0]]) / dx;
        weights[ii][1] = 1.0 - weights[ii][0];
        dWeights_dInput[ii][0] = -1.0 / dx;
        dWeights_dInput[ii][1] = -dWeights_dInput[ii][0];
      }
    }

    // Calculate the result
    for( localIndex ii=0; ii<m_numCorners; ++ii )
    {
      // Find array index
      localIndex tableIndex = 0;
      for( localIndex jj=0; jj<m_dimensions; ++jj )
      {
        tableIndex += bounds[jj][m_corners[jj][ii]] * m_indexIncrement[jj];
      }

      // Determine weighted value
      real64 cornerValue = m_values[tableIndex];
      for( localIndex jj=0; jj<m_dimensions; ++jj )
      {
        dCornerValue_dInput[jj] = cornerValue;
      }

      for( localIndex jj=0; jj<m_dimensions; ++jj )
      {
        cornerValue *= weights[jj][m_corners[jj][ii]];
        for( localIndex kk = 0; kk<m_dimensions; ++kk )
        {
          dCornerValue_dInput[kk] *= ( jj == kk )
                             ? dWeights_dInput[jj][m_corners[jj][ii]]
                       : weights[jj][m_corners[jj][ii]];
        }
      }

      for( localIndex jj=0; jj<m_dimensions; ++jj )
      {
        derivatives[jj] += dCornerValue_dInput[jj];
      }
      value += cornerValue;
    }
  }
  // Nearest, Upper, Lower interpolation methods
  else
  {
    // Determine the index to the nearest table entry
    localIndex tableIndex = 0;
    for( localIndex ii=0; ii<m_dimensions; ++ii )
    {
      // Determine the index along each table axis
      localIndex subIndex = 0;

      if( input[ii] <= m_coordinates[ii][0] )
      {
        // Coordinate is to the left of the table axis
        subIndex = 0;
      }
      else if( input[ii] >= m_coordinates[ii][m_size[ii] - 1] )
      {
        // Coordinate is to the right of the table axis
        subIndex = m_size[ii] - 1;
      }
      else
      {
        // Coordinate is within the table axis
        // Note: std::distance will return the index of the upper table vertex
        auto lower = LvArray::sortedArrayManipulation::find( m_coordinates[ii].begin(), m_coordinates.sizeOfArray( ii ), input[ii] );
        subIndex = LvArray::integerConversion< localIndex >( lower );

        // Interpolation types:
        //   - Nearest returns the value of the closest table vertex
        //   - Upper returns the value of the next table vertex
        //   - Lower returns the value of the previous table vertex
        if( m_interpolationMethod == TableFunction::InterpolationType::Nearest )
        {
          if((input[ii] - m_coordinates[ii][subIndex - 1]) <= (m_coordinates[ii][subIndex] - input[ii]))
          {
            --subIndex;
          }
        }
        else if( m_interpolationMethod == TableFunction::InterpolationType::Lower )
        {
          if( subIndex > 0 )
          {
            --subIndex;
          }
        }
      }

      // Increment the global table index
      tableIndex += subIndex * m_indexIncrement[ii];
    }

    // Retrieve the nearest value
    value = m_values[tableIndex];
  }
}

ENUM_STRINGS( TableFunction::InterpolationType, "linear", "nearest", "upper", "lower" )

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_FUNCTIONS_TABLEFUNCTION_HPP_ */
