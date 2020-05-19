/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FieldSpecificationOps.hpp
 */

#ifndef GEOSX_MANAGERS_FIELDSPECIFICATION_FIELDSPECIFICATIONOPS_HPP
#define GEOSX_MANAGERS_FIELDSPECIFICATION_FIELDSPECIFICATIONOPS_HPP

#include "common/DataTypes.hpp"

namespace geosx
{

struct OpEqual
{
  template< typename T, typename U >
  GEOSX_HOST_DEVICE static inline
  void apply( T & lhs, U const & rhs )
  {
    lhs = static_cast< T >( rhs );
  }
};

struct OpAdd
{
  template< typename T, typename U >
  GEOSX_HOST_DEVICE static inline
  void apply( T & lhs, U const & rhs )
  {
    lhs += static_cast< T >( rhs );
  }
};

template< typename OP >
struct FieldSpecificationOp
{
  using OpType = OP;

  /**
   * @brief Pointwise application of a value to a field
   * @tparam T The type of the array1d field variable specified in @p field.
   * @param[in] field The array1d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component not used.
   * @param[in] value The value to apply to @p field.
   *
   * This function performs field[index] (+)= value.
   */
  template< typename T >
  GEOSX_HOST_DEVICE
  static inline typename std::enable_if< !traits::is_tensorT< T >, void >::type
  SpecifyFieldValue( arrayView1d< T > const & field,
                     localIndex const index,
                     integer const GEOSX_UNUSED_PARAM( component ),
                     real64 const value )
  {
    OP::template apply( field( index ), value );
  }

  /**
   * @brief Pointwise application of value to a field variable.
   * @tparam T The type of the array1d field variable specified in @p field.
   * @param[in] field The array1d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component The component of @p field to apply @p value to. If @p T is a scalar type,
   *                      this will not be used.
   * @param[in] value The value to apply to @p field.
   *
   * This function performs field[index][component] (+)= value.
   */
  template< typename T >
  GEOSX_HOST_DEVICE
  static inline typename std::enable_if< traits::is_tensorT< T >, void >::type
  SpecifyFieldValue( arrayView1d< T > const & field,
                     localIndex const index,
                     integer const component,
                     real64 const value )
  {
    OP::template apply( field( index ).Data()[component], value );
  }

  /**
   * @brief Read a field value
   * @tparam T The type of the array1d field variable specified in @p field.
   * @param[in] field The array1d field variable to read @p value from.
   * @param[in] index The index in field to read @p value from.
   * @param[in] component not used.
   * @param[out] value The value read from @p field.
   *
   * This function performs value (+)= field[index].
   */
  template< typename T >
  GEOSX_HOST_DEVICE
  static inline typename std::enable_if< !traits::is_tensorT< T >, void >::type
  ReadFieldValue( arrayView1d< T const > const & field,
                  localIndex const index,
                  integer const GEOSX_UNUSED_PARAM( component ),
                  real64 & value )
  {
    OP::template apply( value, field( index ) );
  }

  /**
   * @brief Read a field value
   * @tparam T The type of the array1d field variable specified in @p field.
   * @param[in] field The array1d field variable to read @p value from.
   * @param[in] index The index in field to read @p value from.
   * @param[in] component not used.
   * @param[out] value The value read from @p field.
   *
   * This function performs value (+)= field[index][component].
   */
  template< typename T >
  GEOSX_HOST_DEVICE
  static inline typename std::enable_if< traits::is_tensorT< T >, void >::type
  ReadFieldValue( arrayView1d< T const > const & field,
                  localIndex const index,
                  integer const component,
                  real64 & value )
  {
    OP::template apply( value, field( index ).Data()[component] );
  }

  /**
   * @brief Pointwise application of a value to a field variable.
   * @tparam T The type of the array2d field variable specified in @p field.
   * @tparam USD the unit stride dimension of the array @p field.
   * @param[in] field The array2d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component The index along second dimension of 2d array.
   * @param[in] value The value to apply to @p field.
   *
   * This function performs field[index][component] (+)= value.
   */
  template< typename T, int USD >
  GEOSX_HOST_DEVICE
  static inline typename std::enable_if< !traits::is_tensorT< T >, void >::type
  SpecifyFieldValue( arrayView2d< T, USD > const & field,
                     localIndex const index,
                     integer const component,
                     real64 const value )
  {
    if( component >= 0 )
    {
      OP::template apply( field( index, component ), value );
    }
    else
    {
      for( localIndex a = 0; a < field.size( 1 ); ++a )
      {
        OP::template apply( field( index, a ), value );
      }
    }
  }

  /**
   * @brief Pointwise application of a value to a field variable.
   * @tparam T The type of the array2d field variable specified in @p field.
   * @tparam USD the unit stride dimension of the array @p field.
   * @param[in] field The array2d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component The component of @p field to apply @p value to. If @p T is a scalar type,
   *                      this will not be used.
   * @param[in] value The value to apply to @p field.
   *
   * This function performs field[index][component] (+)= value for all values of field[index].
   */
  template< typename T, int USD >
  GEOSX_HOST_DEVICE
  static inline typename std::enable_if< traits::is_tensorT< T >, void >::type
  SpecifyFieldValue( arrayView2d< T, USD > const & field,
                     localIndex const index,
                     integer const component,
                     real64 const value )
  {
    if( component >= 0 )
    {
      for( localIndex a = 0; a < field.size( 1 ); ++a )
      {
        OP::template apply( field( index, a ).Data()[component], value );
      }
    }
    else
    {
      for( localIndex a = 0; a < field.size( 1 ); ++a )
      {
        for( localIndex c = 0; c < T::Length(); ++c )
        {
          OP::template apply( field( index, a ).Data()[c], value );
        }
      }
    }
  }

  /**
   * @brief Read value from a 2d field.
   * @tparam T The type of the array2d field variable specified in @p field.
   * @tparam USD the unit stride dimension of the array @p field.
   * @param[in] field The array2d field variable to read @p value from.
   * @param[in] index The index in field to read @p value from.
   * @param[in] component The index along second dimension of 2d array.
   * @param[out] value The value that is read from @p field.
   *
   * This function performs value (+)= field[index][component].
   */
  template< typename T, int USD >
  GEOSX_HOST_DEVICE
  static inline typename std::enable_if< !traits::is_tensorT< T >, void >::type
  ReadFieldValue( arrayView2d< T const, USD > const & field,
                  localIndex const index,
                  integer const component,
                  real64 & value )
  {
    GEOSX_ASSERT( component >= 0 );
    OP::template apply( value, field( index, component ) );
  }

  /**
   * @brief This function is not meaningful. It exists for generic purposes, but will result in an error if called.
   */
  template< typename T, int USD >
  GEOSX_HOST_DEVICE
  static inline typename std::enable_if< traits::is_tensorT< T >, void >::type
  ReadFieldValue( arrayView2d< T const, USD > const & GEOSX_UNUSED_PARAM( field ),
                  localIndex const GEOSX_UNUSED_PARAM( index ),
                  integer const GEOSX_UNUSED_PARAM( component ),
                  real64 & GEOSX_UNUSED_PARAM( value ) )
  {
    GEOSX_ERROR( "ReadFieldValue: unsupported operation" );
  }

  /**
   * @brief Pointwise application of a value to a field variable.
   * @tparam T The type of the array2d field variable specified in @p field.
   * @tparam USD the unit stride dimension of the array @p field.
   * @param[in] field The array2d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component The index along third dimension of 3d array.
   * @param[in] value The value to apply to @p field.
   *
   * This function performs field[index] (+)= value for all values of field[index].
   */
  template< typename T, int USD >
  GEOSX_HOST_DEVICE
  static inline typename std::enable_if< !traits::is_tensorT< T >, void >::type
  SpecifyFieldValue( arrayView3d< T, USD > const & field,
                     localIndex const index,
                     integer const component,
                     real64 const value )
  {
    if( component >= 0 )
    {
      for( localIndex a = 0; a < field.size( 1 ); ++a )
      {
        OP::template apply( field( index, a, component ), value );
      }
    }
    else
    {
      for( localIndex a = 0; a < field.size( 1 ); ++a )
      {
        for( localIndex b = 0; b < field.size( 2 ); ++b )
        {
          OP::template apply( field( index, a, b ), value );
        }
      }
    }
  }

  /**
   * @brief Pointwise application of a value to a field variable.
   * @tparam T The type of the array2d field variable specified in @p field.
   * @tparam USD the unit stride dimension of the array @p field.
   * @param[in] field The array2d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component The component of @p field to apply @p value to. If @p T is a scalar type,
   *                      this will not be used.
   * @param[in] value The value to apply to @p field.
   *
   * This function performs field[index][component] (+)= value for all values of field[index].
   */
  template< typename T, int USD >
  GEOSX_HOST_DEVICE
  static inline typename std::enable_if< traits::is_tensorT< T >, void >::type
  SpecifyFieldValue( arrayView3d< T, USD > const & field,
                     localIndex const index,
                     integer const component,
                     real64 const value )
  {
    if( component >= 0 )
    {
      for( localIndex a = 0; a < field.size( 1 ); ++a )
      {
        for( localIndex b = 0; b < field.size( 2 ); ++b )
        {
          OP::template apply( field( index, a, b ).Data()[component], value );
        }
      }
    }
    else
    {
      for( localIndex a = 0; a < field.size( 1 ); ++a )
      {
        for( localIndex b = 0; b < field.size( 2 ); ++b )
        {
          for( localIndex c = 0; c < T::Length(); ++c )
          {
            OP::template apply( field( index, a, b ).Data()[c], value );
          }
        }
      }
    }
  }

  /**
   * @brief This function is not meaningful. It exists for generic purposes, but will result in an error if called.
   */
  template< typename T, int USD >
  GEOSX_HOST_DEVICE
  static inline void
  ReadFieldValue( arrayView3d< T const, USD > const & GEOSX_UNUSED_PARAM( field ),
                  localIndex const GEOSX_UNUSED_PARAM( index ),
                  integer const GEOSX_UNUSED_PARAM( component ),
                  real64 & GEOSX_UNUSED_PARAM( value ) )
  {
    GEOSX_ERROR( "ReadFieldValue: unsupported operation" );
  }

};

/**
 * @struct FieldSpecificationEqual
 * this struct a collection of static functions which adhere to an assumed interface for overwriting
 * a value for a field.
 */
struct FieldSpecificationEqual : public FieldSpecificationOp< OpEqual >
{
  using base_type = FieldSpecificationOp< OpEqual >;
  using base_type::SpecifyFieldValue;

  /**
   * @brief Function to apply a Dirichlet like boundary condition to a single dof in a system of
   *        equations.
   * @param[in] dof The degree of freedom that is to be set.
   * @param[inout] matrix A ParallelMatrix object: the system matrix.
   * @param[out] rhs The rhs contribution resulting from the application of the BC.
   * @param[in] bcValue The target value of the Boundary Condition
   * @param[in] fieldValue The current value of the variable to be set.
   *
   * This function clears the matrix row for the specified \p dof, sets the diagonal to some
   * appropriately scaled value, and sets \p rhs to the negative product of the scaled value
   * of the diagonal and the difference between \p bcValue and \p fieldValue.
   *
   * @note This function assumes the user is doing a Newton-type nonlinear solve and will
   * negate the rhs vector upon assembly. Thus, it sets the value to negative of the desired
   * update for the field. For a linear problem, this may lead to unexpected results.
   */
  template< typename LAI >
  static inline void SpecifyFieldValue( globalIndex const dof,
                                        typename LAI::ParallelMatrix & matrix,
                                        real64 & rhs,
                                        real64 const bcValue,
                                        real64 const fieldValue )
  {
    if( matrix.getLocalRowID( dof ) >= 0 )
    {
      real64 const diag = matrix.clearRow( dof, true );
      rhs = -diag * (bcValue - fieldValue);
    }
    else
    {
      rhs = 0.0;
    }
  }

  static inline void GEOSX_HOST_DEVICE
  SpecifyFieldValue( globalIndex const dof,
                     globalIndex const dofRankOffset,
                     LvArray::CRSMatrixView< real64, globalIndex const, localIndex const > const & matrix,
                     real64 & rhs,
                     real64 const bcValue,
                     real64 const fieldValue )
  {
    globalIndex const localRow = dof - dofRankOffset;
    if( localRow >= 0 && localRow < matrix.numRows() )
    {
      arraySlice1d< globalIndex const > const columns = matrix.getColumns( localRow );
      arraySlice1d< real64 > const entries = matrix.getEntries( localRow );
      localIndex const numEntries = matrix.numNonZeros( localRow );

      real64 diagonal = 0;
      for( localIndex j = 0; j < numEntries; ++j )
      {
        if( columns[ j ] == dof )
        { diagonal = entries[ j ]; }
        else
        { entries[ j ] = 0; }
      }

      rhs = -diagonal * (bcValue - fieldValue);
    }
    else
    {
      rhs = 0.0;
    }
  }

  /**
   * @brief Function to add some values of a vector.
   * @param rhs A ParallelVector object.
   * @param num The number of values in \p rhs to replace
   * @param dof A pointer to the global DOF to be replaced
   * @param values A pointer to the values corresponding to \p dof that will be added to \p rhs.
   */
  template< typename LAI >
  static inline void PrescribeRhsValues( typename LAI::ParallelVector & rhs,
                                         localIndex const num,
                                         globalIndex * const dof,
                                         real64 * const values )
  {
    for( localIndex a = 0; a < num; ++a )
    {
      if( rhs.getLocalRowID( dof[a] ) >= 0 )
      {
        rhs.set( dof[a], values[a] );
      }
    }
  }

  template< typename POLICY >
  static inline void PrescribeRhsValues( arrayView1d< real64 > const & rhs,
                                         arrayView1d< globalIndex const > const & dof,
                                         globalIndex const dofRankOffset,
                                         arrayView1d< real64 const > const & values )
  {
    GEOSX_ASSERT_EQ( dof.size(), values.size() );
    forAll< POLICY >( dof.size(), [rhs, dof, dofRankOffset, values] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      globalIndex const localRow = dof[ a ] - dofRankOffset;
      if( localRow >= 0 && localRow < rhs.size() )
      { rhs[ localRow ] = values[ a ]; }
    } );
  }
};

/**
 * @struct FieldSpecificationAdd
 * this struct a collection of static functions which adhere to an assumed interface for adding
 * a value for a field.
 */
struct FieldSpecificationAdd : public FieldSpecificationOp< OpAdd >
{
  using base_type = FieldSpecificationOp< OpAdd >;
  using base_type::SpecifyFieldValue;

  /**
   * @brief Function to apply a value to a vector field for a single dof.
   * @param[in] dof The degree of freedom that is to be modified.
   * @param[in] matrix A ParalleMatrix object: the system matrix.
   * @param[out] rhs The rhs contribution to be modified
   * @param[in] bcValue The value to add to rhs
   * @param[in] fieldValue unused.
   *
   */
  template< typename LAI >
  static inline void SpecifyFieldValue( globalIndex const GEOSX_UNUSED_PARAM( dof ),
                                        typename LAI::ParallelMatrix & GEOSX_UNUSED_PARAM( matrix ),
                                        real64 & rhs,
                                        real64 const bcValue,
                                        real64 const GEOSX_UNUSED_PARAM( fieldValue ) )
  {
    rhs += bcValue;
  }

  /**
   * @brief Function to apply a value to a vector field for a single dof.
   * @param[in] dof The degree of freedom that is to be modified.
   * @param[in] matrix A ParalleMatrix object: the system matrix.
   * @param[out] rhs The rhs contribution to be modified
   * @param[in] bcValue The value to add to rhs
   * @param[in] fieldValue unused.
   *
   */
  GEOSX_HOST_DEVICE
  static inline void SpecifyFieldValue( globalIndex const GEOSX_UNUSED_PARAM( dof ),
                                        globalIndex const GEOSX_UNUSED_PARAM( dofRankOffset ),
                                        LvArray::CRSMatrixView< real64, globalIndex const, localIndex const > const & GEOSX_UNUSED_PARAM( matrix ),
                                        real64 & rhs,
                                        real64 const bcValue,
                                        real64 const GEOSX_UNUSED_PARAM( fieldValue ) )
  {
    rhs += bcValue;
  }

  /**
   * @brief Function to add some values of a vector.
   * @param rhs A ParallelVector object.
   * @param num The number of values in \p rhs to replace
   * @param dof A pointer to the global DOF to be replaced
   * @param values A pointer to the values corresponding to \p dof that will be added to \p rhs.
   */
  template< typename LAI >
  static inline void PrescribeRhsValues( typename LAI::ParallelVector & rhs,
                                         localIndex const num,
                                         globalIndex * const dof,
                                         real64 * const values )
  {
    rhs.add( dof, values, num );
  }

  template< typename POLICY >
  static inline void PrescribeRhsValues( arrayView1d< real64 > const & rhs,
                                         arrayView1d< globalIndex const > const & dof,
                                         globalIndex const dofRankOffset,
                                         arrayView1d< real64 const > const & values )
  {
    GEOSX_ASSERT_EQ( dof.size(), values.size() );
    forAll< POLICY >( dof.size(), [rhs, dof, dofRankOffset, values] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      globalIndex const localRow = dof[ a ] - dofRankOffset;
      if( localRow >= 0 && localRow < rhs.size() )
      { rhs[ localRow ] += values[ a ]; }
    } );
  }

};

} //namespace geosx

#endif //GEOSX_MANAGERS_FIELDSPECIFICATION_FIELDSPECIFICATIONOPS_HPP
