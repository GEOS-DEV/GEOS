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
 * @brief This file contains the definition of the R1TensorT class
 * @file R1TensorT.h
 */

#ifndef R1_TENSOR_T_H_
#define R1_TENSOR_T_H_

#include "TensorBaseT.h"
#include "LvArray/src/ArraySlice.hpp"

#include <cstdlib>

/**
 * @brief R1TensorT is a rank-1 tensor object type
 * @tparam T_dim length of tensor index
 *
 * R1TensorT derives from TensorBaseT, and defines basic operations that can be
 * done on a rank-1 tensor, as well as operations that result in a rank-1
 * tensor.
 */
template< int T_dim >
class R1TensorT : public TensorBaseT< T_dim >
{

public:
  //**** CONSTRUCTORS AND DESTRUCTORS ******************************************

  R1TensorT() = default;
  R1TensorT( R1TensorT const & ) = default;
  ~R1TensorT() = default;

  /**
   * @param[in] data use for initialization of t_data
   */
  explicit R1TensorT( const realT data ): TensorBaseT< T_dim >( data ) {}

  /**
   * @param[in] data naked array used for initialization of t_data
   */
  explicit R1TensorT( realT * const data ): TensorBaseT< T_dim >( data ) {}

  /**
   * @param[in] data use for initialization of t_data
   */
  explicit R1TensorT( const int data ): TensorBaseT< T_dim >( realT( data ) ) {}

  template< int USD >
  R1TensorT( LvArray::ArraySlice< realT const, 1, USD, std::ptrdiff_t > const & src ):
    TensorBaseT< T_dim >()
  { *this = src; }

  //**** CONSTRUCTORS AND DESTRUCTORS
  // *******************************************

  /**
   * Explicit constructors - will throw compile-time errors if not called with
   * the correct dimension
   */
  R1TensorT( realT x, realT y );  //2D only

  GEOSX_HOST_DEVICE
  R1TensorT( realT x, realT y, realT z ); //3D only

  //***** ASSIGNMENT OPERATORS *************************************************
  /// assignment of all data to an integer
  R1TensorT & operator=( const int & rhs );

  /// assignment to all data to a realT
  GEOSX_HOST_DEVICE
  R1TensorT & operator=( const realT & rhs );

  /// assignment to another R1TensorT
  R1TensorT & operator=( const R1TensorT & rhs ) = default;


  template< int USD >
  GEOSX_HOST_DEVICE inline
  R1TensorT & operator=( LvArray::ArraySlice< realT const, 1, USD, std::ptrdiff_t > const & src )
  {
    GEOSX_ASSERT_EQ( src.size(), T_dim );

    for( int i = 0; i < T_dim; ++i )
    {
      this->t_data[ i ] = src[ i ];
    }

    return *this;
  }

  //***** ACCESS OPERATORS ****************************************************
  /// const access to data
  GEOSX_HOST_DEVICE inline realT operator()( const int i ) const { return this->t_data[i]; }

  /// non-const access to data
  GEOSX_HOST_DEVICE inline realT & operator()( const int i )       { return this->t_data[i]; }

  /// const access to data
  GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE realT operator[]( const int i ) const { return this->t_data[i]; }

  /// non-const access to data
  GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE realT & operator[]( const int i )       { return this->t_data[i]; }

private:
};


//*****************************************************************************
//***** END DECLARATION *******************************************************
//*****************************************************************************

/// Explicit 2D constructor
///
/// Template specialisation - if templated on another dimension constructor will
// throw a compile time error.
template<>
inline R1TensorT< 2 >::R1TensorT( realT x, realT y ):
  TensorBaseT< 2 >()
{
  this->t_data[0] = x;
  this->t_data[1] = y;
}


/// Explicit 3D constructor
///
/// Template specialisation - if templated on another dimension constructor will
// throw a compile time error.
template<>
inline R1TensorT< 3 >::R1TensorT( realT x, realT y, realT z ):
  TensorBaseT< 3 >()
{
  this->t_data[0] = x;
  this->t_data[1] = y;
  this->t_data[2] = z;
}



//*****************************************************************************
//***** R1TensorT Member Function Definition **********************************
//*****************************************************************************



//***** ASSIGNMENT OPERATORS **************************************************

/**
 * @param[in] rhs value to set each member of t_data to
 * @return reference to *this
 */
template< int T_dim >
inline R1TensorT< T_dim > & R1TensorT< T_dim >::operator=( const int & rhs )
{
  TensorBaseT< T_dim >::operator=( rhs );
  return *this;
}

/**
 * @param[in] rhs value to set each member of t_data to
 * @return reference to *this
 */
template< int T_dim >
GEOSX_HOST_DEVICE
inline R1TensorT< T_dim > & R1TensorT< T_dim >::operator=( const realT & rhs )
{
  TensorBaseT< T_dim >::operator=( rhs );
  return *this;
}

#endif
