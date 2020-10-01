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
 * @brief This file contains the definition of the TensorBaseT class
 * @file TensorBaseT.h
 */

#ifndef TENSOR_BASE_T_H_
#define TENSOR_BASE_T_H_
#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <exception>
#include <limits>
#include "common/Logger.hpp"
#include "common/GeosxMacros.hpp"

using realT = double;

/**
 * @brief TensorBaseT is the base class for the tensor library.
 * @tparam T_length length
 *
 * TensorBaseT defines basic operations on the data for use by the derived
 * class.
 */
template< int T_length >
class TensorBaseT
{

  friend std ::istream & operator>>( std::istream & in, TensorBaseT< T_length > & t )
  {
    realT *tp = t.Data();
    for( int ii = 0; ii < T_length; ++ii )
    {
      while( in.peek() == ',' || in.peek() == ' ' )
      {
        in.ignore();
      }

      in >> tp[ii];
    }
    return in;
  }

  friend std ::ostream & operator<<( std::ostream & out, const TensorBaseT< T_length > & t )
  {
    const realT *tp = t.Data();
    for( int ii = 0; ii < T_length; ++ii )
    {
      if( ii > 0 ) out << " ";
      out << tp[ii];
    }
    return out;
  }


public:
  //**** CONSTRUCTORS AND DESTRUCTORS ******************************************

  /// default constructor
  GEOSX_HOST_DEVICE
  TensorBaseT( void );

  /// constructor initialized by single value
  GEOSX_HOST_DEVICE
  explicit TensorBaseT( const realT data );

  /// constructor initialized by raw data
  GEOSX_HOST_DEVICE
  explicit TensorBaseT( const realT data[T_length] );

  /// constructor initialized by another TensorBaseT object
  TensorBaseT( const TensorBaseT< T_length > & rhs ) = default;

  /// non-virtual destructor. This means that this class in NOT intended to be
  /// used as a polymorphically.
  ~TensorBaseT() = default;

  //***** ASSIGNMENT OPERATORS *************************************************
  /// assignment of all data to an integer
  GEOSX_HOST_DEVICE
  TensorBaseT & operator=( const int & rhs );

  /// assignment to all data to a realT
  GEOSX_HOST_DEVICE
  TensorBaseT & operator=( const realT & rhs );

  /// assignment to another TensorBaseT
  TensorBaseT & operator=( const TensorBaseT & rhs ) = default;

  bool operator<( const TensorBaseT & rhs ) const
  {
    bool rval = true;
    for( int i = 0; i < T_length; ++i )
    {
      if( t_data[i] >= rhs.t_data[i] )
      {
        rval = false;
      }
    }
    return rval;
  }

  bool operator<=( const TensorBaseT< T_length > & rhs ) const
  {
    bool rval = true;
    for( int i = 0; i < T_length; ++i )
    {
      if( t_data[i] > rhs.t_data[i] )
      {
        rval = false;
      }
    }
    return rval;
  }

  bool operator>( const TensorBaseT< T_length > & rhs ) const
  {
    bool rval = true;
    for( int i = 0; i < T_length; ++i )
    {
      if( t_data[i] <= rhs.t_data[i] )
      {
        rval = false;
      }
    }
    return rval;
  }

  bool operator>=( const TensorBaseT< T_length > & rhs ) const
  {
    bool rval = true;
    for( int i = 0; i < T_length; ++i )
    {
      if( t_data[i] < rhs.t_data[i] )
      {
        rval = false;
      }
    }
    return rval;
  }

  bool operator==( const TensorBaseT< T_length > & rhs ) const
  {
    for( int i = 0; i < T_length; ++i )
    {
      if( (t_data[i] > rhs.t_data[i]) || (t_data[i] < rhs.t_data[i]) )
      {
        return false;
      }
    }
    return true;
  }

//***** DATA MEMBERS ********************************************************
protected:
  /// Tensor data array
  realT t_data[T_length];

  //***** MEMBER ACCESS *******************************************************
public:
  /**
   * @return gives a non-const realT* which points to t_data
   * @brief returns a non-const realT* which points to t_data
   */
  GEOSX_HOST_DEVICE inline
  realT * Data( void )
  {
    return t_data;
  }

  /**
   * @return gives a const realT* which points to t_data
   * @brief gives a const realT* which points to t_data
   */
  GEOSX_HOST_DEVICE inline constexpr
  const realT * Data( void ) const
  {
    return t_data;
  }

  /**
   * @return the number of data entries (Length) for the tensor
   * @brief gives the number of data entries (Length) for the tensor
   */
  GEOSX_HOST_DEVICE constexpr
  static int Length( void )
  {
    return T_length;
  }

private:
  //  TensorBaseT(TensorBaseT<T_length>&);

};
//*****************************************************************************
//***** END DECLARATION *******************************************************
//*****************************************************************************



//*****************************************************************************
//***** TensorBaseT Member Function Definition ********************************
//*****************************************************************************

//**** CONSTRUCTORS AND DESTRUCTORS *******************************************

/**
 * @return none
 */
template< int T_length >
GEOSX_HOST_DEVICE
TensorBaseT< T_length >::TensorBaseT( void )//:
{
  *this = 0.0;
}

/**
 * @param[in] data naked array used for initialization of t_data
 * @return none
 */
template< int T_length >
TensorBaseT< T_length >::TensorBaseT( const realT data )
{
  for( int i = 0; i < T_length; ++i )
    t_data[i] = data;
}

/**
 * @param[in] data naked array used for initialization of t_data
 * @return none
 */
template< int T_length >
TensorBaseT< T_length >::TensorBaseT( const realT data[T_length] )
{
  for( int i = 0; i < T_length; ++i )
    t_data[i] = data[i];
}

//***** ASSIGNMENT OPERATORS **************************************************

/**
 * @param[in] rhs value to set each member of t_data to
 * @return none
 */
template< int T_length >
GEOSX_FORCE_INLINE
TensorBaseT< T_length > &
TensorBaseT< T_length >::operator=( const int & rhs )
{
  operator=( static_cast< realT >( rhs ) );
  return *this;
}

/**
 * @param[in] rhs value to set each member of t_data to
 * @return none
 */
template< int T_length >
GEOSX_FORCE_INLINE
GEOSX_HOST_DEVICE
TensorBaseT< T_length > &
TensorBaseT< T_length >::operator=( const realT & rhs )
{
  for( int i = 0; i < T_length; ++i )
    t_data[i] = rhs;
  return *this;
}

#endif
