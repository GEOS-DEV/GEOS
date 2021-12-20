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
 * An interface for a dense table-based multivariable function (function with multiple inputs and outputs)
 */
class MultivariableTableFunction : public FunctionBase
{
public:

  /// maximum dimensions for the coordinates in the table
  static constexpr integer maxDimensions = 4;



  /**
   * @class KernelWrapperBase
   *
   * A nested base class for the kernel function doing the interpolation in the table
   */
  class KernelWrapperBase
  {
public:

    /// @cond DO_NOT_DOCUMENT
    /// We need these SMFs to enable array1d< KernelWrapper > and avoid
    /// host-device errors with CUDA. Otherwise rule of 0 would be fine.

    KernelWrapperBase() = default;
    /// @endcond

  };


  /**
   * @class KernelWrapper
   *
   * A nested class encapsulating the kernel function doing the interpolation in the table
   */
  template< integer NUM_DIMS, integer NUM_OPS >
  class KernelWrapper : public KernelWrapperBase
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
      return *this;
    }

    /// @endcond
    /**
     * @brief Interpolate in the table.
     * @tparam IN_ARRAY type of input value array
     * @param[in] input vector of input value
     * @return interpolated value
     */
    // template< typename IN_ARRAY >
    // GEOSX_HOST_DEVICE
    // real64 compute( IN_ARRAY const & input ) const
    // {
    //   return interpolateLinear( input );
    // }

    GEOSX_HOST_DEVICE
    void compute( localIndex const ei ) const
    {
      interpolateLinear( m_input );
    }

    /**
     * @brief Interpolate in the table.
     * @tparam IN_ARRAY type of input value array
     * @param[in] input vector of input value
     * @param[out] values vector of interpolated values
     *
     * @return interpolated value
     */
    template< typename IN_ARRAY, typename OUT_ARRAY >
    GEOSX_HOST_DEVICE
    void compute( IN_ARRAY const & input, OUT_ARRAY && values ) const
    {
      for( integer op = 0; op < NUM_OPS; ++op )
      {
        values[op] = interpolateLinear( input );
      }
    }

    /**
     * @brief Interpolate in the table with derivatives.
     * @param[in] input vector of input value
     * @param[out] values vector of interpolated values
     * @param[out] derivatives vector of derivatives of interpolated values wrt the variables present in input
     * @return interpolated value
     */
    template< typename IN_ARRAY, typename OUT_ARRAY >
    GEOSX_HOST_DEVICE
    real64 compute( IN_ARRAY const & input, OUT_ARRAY && values, OUT_ARRAY && derivatives ) const;

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

    template< typename POLICY, typename KERNEL_TYPE >
    static void
    launch( localIndex const numElems,
            KERNEL_TYPE const & kernelComponent )
    {
      forAll< POLICY >( numElems, [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        kernelComponent.compute( ei );
      } );
    }


private:

    friend class MultivariableTableFunction; // Allow only parent class to construct the wrapper

    /**
     * @brief The constructor of the kernel wrapper
     * @param[in] coordinates array of table axes
     * @param[in] values table values (in fortran order)
     */
    KernelWrapper( arrayView1d< real64 const > const & input,
                   ArrayOfArraysView< real64 const > const & coordinates,
                   arrayView1d< real64 > const & values );
    /**
     * @brief Interpolate in the table using linear method.
     * @param[in] input vector of input value
     * @return interpolated value
     */
    template< typename IN_ARRAY >
    GEOSX_HOST_DEVICE
    real64
    interpolateLinear( IN_ARRAY const & input ) const
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

    /// Input array view
    arrayView1d< real64 const > m_input;

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
    FunctionBase::evaluateT< MultivariableTableFunction >( group, time, set, result );
  }

  /**
   * @brief Method to evaluate a function
   * @param input a scalar input
   * @return the function result
   */
  virtual real64 evaluate( real64 const * const input ) const override final;


  real64 evaluate( arrayView1d< real64 const > const & input ) const;

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
   * @brief Set the table coordinates
   * @param coordinates An array of arrays containing table coordinate definitions
   */
  void setTableCoordinates( array1d< real64_array > const & coordinates );

  /**
   * @brief Set the table values
   * @param values An array of table values in fortran order
   */
  void setTableValues( real64_array values );

  // /**
  //  * @brief Create an instance of the kernel wrapper
  //  * @return the kernel wrapper
  //  */
  // KernelWrapper< 2, 3 > createKernelWrapper() const;

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

  template< typename POLICY >
  static void
  createAndLaunch( arrayView1d< real64 const > const & input,
                   ArrayOfArraysView< real64 const > const & coordinates,
                   arrayView1d< real64 > const & values );

  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent );



  /// Coordinates for 1D table
  array1d< real64 > m_tableCoordinates1D;

  /// List of table coordinate file names
  path_array m_coordinateFiles;

  /// Table voxel file names
  Path m_voxelFile;

  /// An array of table axes
  ArrayOfArrays< real64 > m_coordinates;

  /// Table values (in fortran order)
  array1d< real64 > m_values;

  /// Kernel wrapper object used in evaluate() interface
  //KernelWrapper< 2, 3 > m_kernelWrapper;

};

template< integer NUM_DIMS, integer NUM_OPS >
MultivariableTableFunction::KernelWrapper< NUM_DIMS, NUM_OPS >::KernelWrapper( arrayView1d< real64 const > const & input,
                                                                               ArrayOfArraysView< real64 const > const & coordinates,
                                                                               arrayView1d< real64 > const & values )
  :
  m_input( input ),
  m_coordinates( coordinates ),
  m_values( values )
{}

// template< typename IN_ARRAY, integer NUM_DIMS, integer NUM_OPS >
// GEOSX_HOST_DEVICE
// real64
// MultivariableTableFunction::KernelWrapper< NUM_DIMS, NUM_OPS >::compute( IN_ARRAY const & input ) const
// {
//   return interpolateLinear( input );
// }

// template< typename IN_ARRAY, integer NUM_DIMS, integer NUM_OPS >
// GEOSX_HOST_DEVICE
// real64
// MultivariableTableFunction::KernelWrapper< NUM_DIMS, NUM_OPS >::interpolateLinear( IN_ARRAY const & input ) const
// {}

// template< typename IN_ARRAY, typename OUT_ARRAY, integer NUM_DIMS, integer NUM_OPS >
// GEOSX_HOST_DEVICE
// void
// MultivariableTableFunction::KernelWrapper< NUM_DIMS, NUM_OPS >::compute( IN_ARRAY const & input, OUT_ARRAY && derivatives ) const
// {}


} /* namespace geosx */

#endif /* GEOSX_FUNCTIONS_MULTIVARIABLETABLEFUNCTION_HPP_ */
