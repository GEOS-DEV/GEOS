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
 * @file TableFunction.hpp
 */

#ifndef GEOSX_FUNCTIONS_TABLEFUNCTION_HPP_
#define GEOSX_FUNCTIONS_TABLEFUNCTION_HPP_

#include "FunctionBase.hpp"

#include "codingUtilities/EnumStrings.hpp"
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
  static constexpr integer maxDimensions = 4;

  /**
   * @class KernelWrapper
   *
   * A nested class encapsulating the kernel function doing the interpolation in the table
   */
  class KernelWrapper
  {
public:

    /// @cond DO_NOT_DOCUMENT
    /// We need these SMFs to enable array1d< KernelWrapper > and avoid
    /// host-device errors with CUDA. Otherwise rule of 0 would be fine.

    KernelWrapper() = default;
    KernelWrapper( KernelWrapper const & ) = default;
    KernelWrapper( KernelWrapper && ) = default;
    KernelWrapper & operator=( KernelWrapper const & ) = default;

    /// Note: move assignment not deleted, not defaulted on purpose!
    /// This is needed to avoid a compilation warning with CUDA
    KernelWrapper & operator=( KernelWrapper && other )
    {
      m_coordinates = std::move( other.m_coordinates );
      m_values = std::move( other.m_values );
      m_interpolationMethod = other.m_interpolationMethod;
      return *this;
    }

    /// @endcond

    /**
     * @brief Interpolate in the table.
     * @tparam IN_ARRAY type of input value array
     * @param[in] input vector of input value
     * @return interpolated value
     */
    template< typename IN_ARRAY >
    GEOSX_HOST_DEVICE
    real64 compute( IN_ARRAY const & input ) const;

    /**
     * @brief Interpolate in the table with derivatives.
     * @param[in] input vector of input value
     * @param[out] derivatives vector of derivatives of interpolated value wrt the variables present in input
     * @return interpolated value
     */
    template< typename IN_ARRAY, typename OUT_ARRAY >
    GEOSX_HOST_DEVICE
    real64 compute( IN_ARRAY const & input, OUT_ARRAY && derivatives ) const;

    /**
     * @brief Move the KernelWrapper to the given execution space, optionally touching it.
     * @param space the space to move the KernelWrapper to
     * @param touch whether the KernelWrapper should be touched in the new space or not
     * @note This function exists to enable holding KernelWrapper objects in an ArrayView
     *       and have their contents properly moved between memory spaces.
     */
    void move( LvArray::MemorySpace const space, bool const touch )
    {
      m_coordinates.move( space, touch );
      m_values.move( space, touch );
    }

private:

    friend class TableFunction; // Allow only parent class to construct the wrapper

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
    KernelWrapper( InterpolationType interpolationMethod,
                   ArrayOfArraysView< real64 const > const & coordinates,
                   arrayView1d< real64 const > const & values );

    /**
     * @brief Interpolate in the table using linear method.
     * @param[in] input vector of input value
     * @return interpolated value
     */
    template< typename IN_ARRAY >
    GEOSX_HOST_DEVICE
    real64
    interpolateLinear( IN_ARRAY const & input ) const;

    /**
     * @brief Interpolate in the table with derivatives using linear method.
     * @param[in] input vector of input value
     * @param[out] derivatives vector of derivatives of interpolated value wrt the variables present in input
     * @return interpolated value
     */
    template< typename IN_ARRAY, typename OUT_ARRAY >
    GEOSX_HOST_DEVICE
    real64
    interpolateLinear( IN_ARRAY const & input, OUT_ARRAY && derivatives ) const;

    /**
     * @brief Interpolate in the table by rounding to an exact point
     * @param[in] input vector of input value
     * @return interpolated value
     */
    template< typename IN_ARRAY >
    GEOSX_HOST_DEVICE
    real64
    interpolateRound( IN_ARRAY const & input ) const;

    /**
     * @brief Interpolate in the table with derivatives using linear method.
     * @param[in] input vector of input value
     * @param[out] derivatives vector of derivatives of interpolated value wrt the variables present in input
     * @return interpolated value
     */
    template< typename IN_ARRAY, typename OUT_ARRAY >
    GEOSX_HOST_DEVICE
    real64
    interpolateRound( IN_ARRAY const & input, OUT_ARRAY && derivatives ) const;

    /// Table interpolation method
    TableFunction::InterpolationType m_interpolationMethod = InterpolationType::Linear;

    /// An array of table axes
    ArrayOfArraysView< real64 const > m_coordinates;

    /// Table values (in fortran order)
    arrayView1d< real64 const > m_values;
  };

  /**
   * @brief The constructor
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  TableFunction( const string & name,
                 dataRepository::Group * const parent );

  /**
   * @brief The catalog name interface
   * @return name of the TableFunction in the FunctionBase catalog
   */
  static string catalogName() { return "TableFunction"; }

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
   * @brief @return Number of table dimensions
   */
  integer numDimensions() const { return LvArray::integerConversion< integer >( m_coordinates.size() ); }

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
  arrayView1d< real64 const > getValues() const { return m_values.toViewConst(); }

  /**
   * @copydoc getValues() const
   */
  array1d< real64 > & getValues() { return m_values; }

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
  void setTableCoordinates( array1d< real64_array > const & coordinates );

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

  /// Struct containing lookup keys for data repository wrappers
  struct viewKeyStruct
  {
    /// @return Key for coordinate arrays
    static constexpr char const * coordinatesString() { return "coordinates"; }
    /// @return Key for value array
    static constexpr char const * valuesString() { return "values"; }
    /// @return Key for interpolation type
    static constexpr char const * interpolationString() { return "interpolation"; }
    /// @return Key for list of files containing table coordinates
    static constexpr char const * coordinateFilesString() { return "coordinateFiles"; }
    /// @return Key for name of file containing table values
    static constexpr char const * voxelFileString() { return "voxelFile"; }
  };

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

  /// Coordinates for 1D table
  array1d< real64 > m_tableCoordinates1D;

  /// List of table coordinate file names
  path_array m_coordinateFiles;

  /// Table voxel file names
  Path m_voxelFile;

  /// Table interpolation method
  InterpolationType m_interpolationMethod;

  /// An array of table axes
  ArrayOfArrays< real64 > m_coordinates;

  /// Table values (in fortran order)
  array1d< real64 > m_values;

  /// Kernel wrapper object used in evaluate() interface
  KernelWrapper m_kernelWrapper;

};

template< typename IN_ARRAY >
GEOSX_HOST_DEVICE
real64
TableFunction::KernelWrapper::compute( IN_ARRAY const & input ) const
{
  if( m_interpolationMethod == TableFunction::InterpolationType::Linear )
  {
    return interpolateLinear( input );
  }
  else // Nearest, Upper, Lower interpolation methods
  {
    return interpolateRound( input );
  }
}

template< typename IN_ARRAY >
GEOSX_HOST_DEVICE
real64
TableFunction::KernelWrapper::interpolateLinear( IN_ARRAY const & input ) const
{
  integer const numDimensions = LvArray::integerConversion< integer >( m_coordinates.size() );
  localIndex bounds[maxDimensions][2]{};
  real64 weights[maxDimensions][2]{};

  // Determine position, weights
  for( localIndex dim = 0; dim < numDimensions; ++dim )
  {
    arraySlice1d< real64 const > const coords = m_coordinates[dim];
    if( input[dim] <= coords[0] )
    {
      // Coordinate is to the left of this axis
      bounds[dim][0] = 0;
      bounds[dim][1] = 0;
      weights[dim][0] = 0.0;
      weights[dim][1] = 1.0;
    }
    else if( input[dim] >= coords[coords.size() - 1] )
    {
      // Coordinate is to the right of this axis
      bounds[dim][0] = coords.size() - 1;
      bounds[dim][1] = bounds[dim][0];
      weights[dim][0] = 1.0;
      weights[dim][1] = 0.0;
    }
    else
    {
      // Find the coordinate index
      ///TODO make this fast
      // Note: find uses a binary search...  If we assume coordinates are
      // evenly spaced, we can speed things up considerably
      auto const lower = LvArray::sortedArrayManipulation::find( coords.begin(), coords.size(), input[dim] );
      bounds[dim][1] = LvArray::integerConversion< localIndex >( lower );
      bounds[dim][0] = bounds[dim][1] - 1;

      real64 const dx = coords[bounds[dim][1]] - coords[bounds[dim][0]];
      weights[dim][0] = 1.0 - ( input[dim] - coords[bounds[dim][0]]) / dx;
      weights[dim][1] = 1.0 - weights[dim][0];
    }
  }

  // Calculate the result
  real64 value = 0.0;
  integer const numCorners = 1 << numDimensions;
  for( integer point = 0; point < numCorners; ++point )
  {
    // Find array index
    localIndex tableIndex = 0;
    localIndex stride = 1;
    for( integer dim = 0; dim < numDimensions; ++dim )
    {
      integer const corner = (point >> dim) & 1;
      tableIndex += bounds[dim][corner] * stride;
      stride *= m_coordinates.sizeOfArray( dim );
    }

    // Determine weighted value
    real64 cornerValue = m_values[tableIndex];
    for( integer dim = 0; dim < numDimensions; ++dim )
    {
      integer const corner = (point >> dim) & 1;
      cornerValue *= weights[dim][corner];
    }
    value += cornerValue;
  }
  return value;
}

template< typename IN_ARRAY >
GEOSX_HOST_DEVICE
real64
TableFunction::KernelWrapper::interpolateRound( IN_ARRAY const & input ) const
{
  integer const numDimensions = LvArray::integerConversion< integer >( m_coordinates.size() );

  // Determine the index to the nearest table entry
  localIndex tableIndex = 0;
  localIndex stride = 1;
  for( integer dim = 0; dim < numDimensions; ++dim )
  {
    arraySlice1d< real64 const > const coords = m_coordinates[dim];
    // Determine the index along each table axis
    localIndex subIndex;
    if( input[dim] <= coords[0] )
    {
      // Coordinate is to the left of the table axis
      subIndex = 0;
    }
    else if( input[dim] >= coords[coords.size() - 1] )
    {
      // Coordinate is to the right of the table axis
      subIndex = coords.size() - 1;
    }
    else
    {
      // Coordinate is within the table axis
      // Note: find() will return the index of the upper table vertex
      auto const lower = LvArray::sortedArrayManipulation::find( coords.begin(), coords.size(), input[dim] );
      subIndex = LvArray::integerConversion< localIndex >( lower );

      // Interpolation types:
      //   - Nearest returns the value of the closest table vertex
      //   - Upper returns the value of the next table vertex
      //   - Lower returns the value of the previous table vertex
      if( m_interpolationMethod == TableFunction::InterpolationType::Nearest )
      {
        if( ( input[dim] - coords[subIndex - 1]) <= ( coords[subIndex] - input[dim]) )
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
    tableIndex += subIndex * stride;
    stride *= coords.size();
  }

  // Retrieve the nearest value
  return m_values[tableIndex];
}

template< typename IN_ARRAY, typename OUT_ARRAY >
GEOSX_HOST_DEVICE
real64
TableFunction::KernelWrapper::compute( IN_ARRAY const & input, OUT_ARRAY && derivatives ) const
{
  // Linear interpolation
  if( m_interpolationMethod == TableFunction::InterpolationType::Linear )
  {
    return interpolateLinear( input, derivatives );
  }
  // Nearest, Upper, Lower interpolation methods
  else
  {
    return interpolateRound( input, derivatives );
  }
}

template< typename IN_ARRAY, typename OUT_ARRAY >
GEOSX_HOST_DEVICE
real64
TableFunction::KernelWrapper::interpolateLinear( IN_ARRAY const & input, OUT_ARRAY && derivatives ) const
{
  integer const numDimensions = LvArray::integerConversion< integer >( m_coordinates.size() );

  localIndex bounds[maxDimensions][2]{};
  real64 weights[maxDimensions][2]{};
  real64 dWeights_dInput[maxDimensions][2]{};

  // Determine position, weights
  for( integer dim = 0; dim < numDimensions; ++dim )
  {
    arraySlice1d< real64 const > const coords = m_coordinates[dim];
    if( input[dim] <= coords[0] )
    {
      // Coordinate is to the left of this axis
      bounds[dim][0] = 0;
      bounds[dim][1] = 0;
      weights[dim][0] = 0;
      weights[dim][1] = 1;
      dWeights_dInput[dim][0] = 0;
      dWeights_dInput[dim][1] = 0;
    }
    else if( input[dim] >= coords[coords.size() - 1] )
    {
      // Coordinate is to the right of this axis
      bounds[dim][0] = coords.size() - 1;
      bounds[dim][1] = bounds[dim][0];
      weights[dim][0] = 1;
      weights[dim][1] = 0;
      dWeights_dInput[dim][0] = 0;
      dWeights_dInput[dim][1] = 0;
    }
    else
    {
      // Find the coordinate index
      ///TODO make this fast
      // Note: find uses a binary search...  If we assume coordinates are
      // evenly spaced, we can speed things up considerably
      auto lower = LvArray::sortedArrayManipulation::find( coords.begin(), coords.size(), input[dim] );
      bounds[dim][1] = LvArray::integerConversion< localIndex >( lower );
      bounds[dim][0] = bounds[dim][1] - 1;

      real64 const dx = coords[bounds[dim][1]] - coords[bounds[dim][0]];
      weights[dim][0] = 1.0 - ( input[dim] - coords[bounds[dim][0]]) / dx;
      weights[dim][1] = 1.0 - weights[dim][0];
      dWeights_dInput[dim][0] = -1.0 / dx;
      dWeights_dInput[dim][1] = -dWeights_dInput[dim][0];
    }
  }

  // Calculate the result
  real64 value = 0.0;
  for( integer dim = 0; dim < numDimensions; ++dim )
  {
    derivatives[dim] = 0.0;
  }

  integer const numCorners = 1 << numDimensions;
  for( integer point = 0; point < numCorners; ++point )
  {
    // Find array index
    localIndex tableIndex = 0;
    localIndex stride = 1;
    for( integer dim = 0; dim < numDimensions; ++dim )
    {
      integer const corner = (point >> dim) & 1;
      tableIndex += bounds[dim][corner] * stride;
      stride *= m_coordinates.sizeOfArray( dim );
    }

    // Determine weighted value
    real64 cornerValue = m_values[tableIndex];
    real64 dCornerValue_dInput[maxDimensions];
    for( integer dim = 0; dim < numDimensions; ++dim )
    {
      dCornerValue_dInput[dim] = cornerValue;
    }

    for( integer dim = 0; dim < numDimensions; ++dim )
    {
      integer const corner = (point >> dim) & 1;
      cornerValue *= weights[dim][corner];
      for( integer kk = 0; kk < numDimensions; ++kk )
      {
        dCornerValue_dInput[kk] *= ( dim == kk ) ? dWeights_dInput[dim][corner] : weights[dim][corner];
      }
    }

    for( integer dim = 0; dim < numDimensions; ++dim )
    {
      derivatives[dim] += dCornerValue_dInput[dim];
    }
    value += cornerValue;
  }
  return value;
}

template< typename IN_ARRAY, typename OUT_ARRAY >
GEOSX_HOST_DEVICE
real64
TableFunction::KernelWrapper::interpolateRound( IN_ARRAY const & input, OUT_ARRAY && derivatives ) const
{
  GEOSX_UNUSED_VAR( input, derivatives );
  GEOSX_ERROR( "Rounding interpolation with derivatives not implemented" );
  return 0.0;
}

/// Declare strings associated with enumeration values.
ENUM_STRINGS( TableFunction::InterpolationType,
              "linear",
              "nearest",
              "upper",
              "lower" );

} /* namespace geosx */

#endif /* GEOSX_FUNCTIONS_TABLEFUNCTION_HPP_ */
