/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @brief This file contains the definition of the TensorBaseT class
 * @file TensorBaseT.h
 * @author Randolph Settgast
 */

#ifndef TENSOR_BASE_T_H_
#define TENSOR_BASE_T_H_
#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <exception>
#include <limits>
#include "TensorOps.h"
#include "Logger.hpp"

/**
 * @brief TensorBaseT is the base class for the tensor library.
 * @author Randolph Settgast
 * @tparam T_length length
 *
 * TensorBaseT defines basic operations on the data for use by the derived
 * class.
 */
template<int T_length>
class TensorBaseT
{

  friend std ::istream& operator>>(std::istream& in, TensorBaseT<T_length>& t){
    realT *tp = t.Data();
    for(int ii = 0 ; ii < T_length ; ++ii)
    {
      while(in.peek() == ',' || in.peek() == ' ')
      {
        in.ignore();
      }

      in >> tp[ii];
    }
    return in;
  }

  friend std ::ostream& operator<<(std::ostream& out, const TensorBaseT<T_length>& t){
    const realT *tp = t.Data();
    for(int ii = 0 ; ii < T_length ; ++ii)
    {
      if(ii > 0) out << " ";
      out << tp[ii];
    }
    return out;
  }


public:
  //**** CONSTRUCTORS AND DESTRUCTORS ******************************************

  /// default constructor
  TensorBaseT( void );

  /// constructor initialized by single value
  explicit TensorBaseT( const realT data );

  /// constructor initialized by raw data
  explicit TensorBaseT( const realT data[T_length] );

  /// constructor initialized by another TensorBaseT object
  TensorBaseT( const TensorBaseT< T_length >& rhs );

  /// non-virtual destructor. This means that this class in NOT intended to be
  /// used as a polymorphically.
  ~TensorBaseT( void );

  //***** ASSIGNMENT OPERATORS *************************************************
  /// assignment of all data to an integer
  TensorBaseT< T_length >& operator=( const int& rhs );

  /// assignment to all data to a realT
  TensorBaseT< T_length >& operator=( const realT& rhs );

  /// assignment to another TensorBaseT
  TensorBaseT< T_length >& operator=( const TensorBaseT< T_length >& rhs );

  /// add a realT to data
  TensorBaseT< T_length >& operator+=( const realT& rhs );

  /// subtract a realT from data
  TensorBaseT< T_length >& operator-=( const realT& rhs );

  /// multiply each entry in t_data by a realT
  TensorBaseT< T_length >& operator*=( const realT& rhs );

  /// divide each entry in t_data by a realT
  TensorBaseT< T_length >& operator/=( const realT& rhs );

  /// add another tensor
  TensorBaseT< T_length >& operator+=( const TensorBaseT< T_length >& rhs );

  /// subtract a tensor
  TensorBaseT< T_length >& operator-=( const TensorBaseT< T_length >& rhs );

  /// multiply by a tensor (data component by component)
  TensorBaseT< T_length >& operator*=( const TensorBaseT< T_length >& rhs );

  /// divide by a tensor (data component by component)
  TensorBaseT< T_length >& operator/=( const TensorBaseT< T_length >& rhs );

  bool operator<( const TensorBaseT< T_length >& rhs ) const
  {
    bool rval = true;
    for (int i = 0 ; i < T_length ; ++i)
    {
      if( t_data[i] >= rhs.t_data[i] )
      {
        rval = false;
      }
    }
    return rval;
  }
  bool operator<=( const TensorBaseT< T_length >& rhs ) const
  {
    bool rval = true;
    for (int i = 0 ; i < T_length ; ++i)
    {
      if( t_data[i] > rhs.t_data[i] )
      {
        rval = false;
      }
    }
    return rval;
  }
  bool operator>( const TensorBaseT< T_length >& rhs ) const
  {
    bool rval = true;
    for (int i = 0 ; i < T_length ; ++i)
    {
      if( t_data[i] <= rhs.t_data[i] )
      {
        rval = false;
      }
    }
    return rval;
  }

  bool operator>=( const TensorBaseT< T_length >& rhs ) const
  {
    bool rval = true;
    for (int i = 0 ; i < T_length ; ++i)
    {
      if( t_data[i] < rhs.t_data[i] )
      {
        rval = false;
      }
    }
    return rval;
  }


  bool operator==( const TensorBaseT< T_length >& rhs ) const
  {
    for (int i = 0; i < T_length; ++i)
    {
      if( (t_data[i] > rhs.t_data[i]) || (t_data[i] < rhs.t_data[i]) )
      {
        return false;
      }
    }
    return true;
  }


  /// function to add the product of a scalar and tensor
  inline void plus_cA( const realT& c, const TensorBaseT< T_length >& A )
  {
    for (int i = 0 ; i < T_length ; ++i)
      t_data[i] += c * A.t_data[i];

  }

  /// function to take the product of a scalar and tensor
  inline void cA( const realT& c, const TensorBaseT< T_length >& A )
  {
    for (int i = 0 ; i < T_length ; ++i)
      t_data[i] = c * A.t_data[i];

  }

  /// function to take the quotient of a tensor by a scalar
  inline void Adivc( const realT& c, const TensorBaseT< T_length >& A )
  {
    for (int i = 0 ; i < T_length ; ++i)
      t_data[i] = A.t_data[i] / c;

  }

  //***** OUTPUT **************************************************************

  /// function to cast data array to float
  void CastDataToFloat( float rval[T_length] ) const
  {
    for( int a=0 ; a<T_length ; ++a )
    {
      rval[a] = static_cast<float>(t_data[a]);
    }
  }

  inline void StrVal(const std::string& str)
  {
    std::istringstream iss(str, std::istringstream::in);
    for ( int i = 0 ; i < T_length ; i++ )
      GEOS_ERROR_IF(!(iss >> t_data[i]), "Error");
  }

/*
   /// ouput function
   virtual void print( ostream& os ) const = 0;

   /// stream function
   friend ostream &operator<<( ostream &os, const TensorBaseT< T_length >& A )
   {
    A.print( os );
    return os;
   }
 */
//***** DATA MEMBERS ********************************************************
protected:
  /// Tensor data array
  realT t_data[T_length];

  //***** MEMBER ACCESS *******************************************************
public:
  /**
   * @author Randolph R. Settgast
   * @return gives a non-const realT* which points to t_data
   * @brief returns a non-const realT* which points to t_data
   */
  realT* Data( void )
  {
    return t_data;
  }

  /**
   * @author Randolph R. Settgast
   * @return gives a const realT* which points to t_data
   * @brief gives a const realT* which points to t_data
   */
  const realT* Data( void ) const
  {
    return t_data;
  }

  /**
   * @author walsh24
   * @return gives a non-const realT* which points to t_data
   * @brief returns a non-const realT* which points to t_data
   */
  realT* begin( void )
  {
    return t_data;
  }

  /**
   * @author walsh24
   * @return gives a const realT* which points to t_data
   * @brief gives a const realT* which points to t_data
   */
  const realT* begin( void ) const
  {
    return t_data;
  }

  /**
   * @author walsh24
   * @return gives a non-const realT* which points to the past-the-end element
   * of t_data
   * @brief returns a non-const realT* which points to the past-the-end element
   * of t_data
   */
  realT* end( void )
  {
    return t_data+T_length;
  }

  /**
   * @author walsh24
   * @return gives a const realT* which points to the past-the-end element of
   * t_data
   * @brief gives a const realT* which points to the past-the-end element of
   * t_data
   */
  const realT* end( void ) const
  {
    return t_data+T_length;
  }



  /**
   * @author Randolph R. Settgast
   * @return the number of data entries (Length) for the tensor
   * @brief gives the number of data entries (Length) for the tensor
   */
  static int Length( void )
  {
    return T_length;
  }

  /**
   * @author Randolph R. Settgast
   * @return the maximum single value in the t_data
   * @brief gives the maximum single value in the t_data
   */
  realT MaxVal( void ) const
  {
    realT rval = 0;
    for (int i = 0 ; i < T_length ; ++i)
      if (fabs( t_data[i] ) > rval)
        rval = fabs( t_data[i] );
    return rval;
  }

  /**
   * @author Scott Johnson
   * @return the minimum single value in the t_data
   * @brief gives the minimum single value in the t_data
   */
  realT MinVal( void ) const
  {
    realT rval = std::numeric_limits<realT>::max();
    for (int i = 0 ; i < T_length ; ++i)
      if (fabs( t_data[i] ) < rval)
        rval = fabs( t_data[i] );
    return rval;
  }

  void SetMax( const TensorBaseT<T_length>& newval )
  {
    for (int i = 0 ; i < T_length ; ++i)
      t_data[i] = (t_data[i] < newval.t_data[i]) ? newval.t_data[i] : t_data[i];
  }

  void SetMin( const TensorBaseT<T_length>& newval )
  {
    for (int i = 0 ; i < T_length ; ++i)
      t_data[i] = (t_data[i] > newval.t_data[i]) ? newval.t_data[i] : t_data[i];
  }

  friend inline
  realT Dot( const TensorBaseT<T_length>& A,  const TensorBaseT<T_length>& B )
  {
    realT rval = 0;
    for( int i=0 ; i<T_length ; ++i )
    {
      rval += A.t_data[i] * B.t_data[i];
    }
    return rval;
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
 * @author Randolph Settgast
 * @return none
 */
template<int T_length>
TensorBaseT< T_length >::TensorBaseT( void )//:
{
  *this = 0.0;
}


/**
 * @author Randolph Settgast
 * @param[in] rhs reference to TensorBaseT object to use in initialization
 * @return none
 */
template<int T_length>
TensorBaseT< T_length >::TensorBaseT( const TensorBaseT< T_length >& rhs )
{
  TensorBaseT< T_length >::operator=( rhs );
}


/**
 * @author Randolph Settgast
 * @param[in] data naked array used for initialization of t_data
 * @return none
 */
template<int T_length>
TensorBaseT< T_length >::TensorBaseT( const realT data )
{
  for (int i = 0 ; i < T_length ; ++i)
    t_data[i] = data;
}

/**
 * @author Randolph Settgast
 * @param[in] data naked array used for initialization of t_data
 * @return none
 */
template<int T_length>
TensorBaseT< T_length >::TensorBaseT( const realT data[T_length] )
{
  for (int i = 0 ; i < T_length ; ++i)
    t_data[i] = data[i];
}

/**
 * @author Randolph Settgast
 * @return none
 */
template<int T_length>
TensorBaseT< T_length >::~TensorBaseT( void )
{}

//***** ASSIGNMENT OPERATORS **************************************************

/**
 * @author Randolph Settgast
 * @param[in] rhs value to set each member of t_data to
 * @return none
 */
template<int T_length>
inline TensorBaseT< T_length >&
TensorBaseT< T_length >::operator=( const int& rhs )
{
  operator=( static_cast< realT > ( rhs ) );
  return *this;
}

/**
 * @author Randolph Settgast
 * @param[in] rhs value to set each member of t_data to
 * @return none
 */
template<int T_length>
inline TensorBaseT< T_length >&
TensorBaseT< T_length >::operator=( const realT& rhs )
{
  for (int i = 0 ; i < T_length ; ++i)
    t_data[i] = rhs;
  return *this;
}

/**
 * @author Randolph Settgast
 * @param[in] rhs tensor to copy
 * @return none
 */
template<int T_length>
inline TensorBaseT< T_length >&
TensorBaseT< T_length >::operator=( const TensorBaseT< T_length >& rhs )
{
  for (int i = 0 ; i < T_length ; ++i)
    t_data[i] = rhs.t_data[i];

//  memcpy(t_data,rhs.t_data,sizeof(realT)*T_length);

  return *this;
}


// intel compiler doesn't seem to be unrolling these loops
template<>
inline TensorBaseT< 3 >&
TensorBaseT<3>::operator=( const TensorBaseT< 3 >& rhs )
{
  t_data[0] = rhs.t_data[0];
  t_data[1] = rhs.t_data[1];
  t_data[2] = rhs.t_data[2];
  return *this;
}


/**
 * @author Randolph Settgast
 * @param[in] rhs value to add to t_data
 * @return none
 */
template<int T_length>
inline TensorBaseT< T_length >&
TensorBaseT< T_length >::operator+=( const realT& rhs )
{
  for (int i = 0 ; i < T_length ; ++i)
    t_data[i] += rhs;
  return *this;
}

template<>
inline TensorBaseT< 3 >&
TensorBaseT<3>::operator+=( const realT& rhs )
{
  t_data[0] += rhs;
  t_data[1] += rhs;
  t_data[2] += rhs;
  return *this;
}

/**
 * @author Randolph Settgast
 * @param[in] rhs value to subtract from t_data
 * @return none
 */
template<int T_length>
inline TensorBaseT< T_length >&
TensorBaseT< T_length >::operator-=( const realT& rhs )
{
  for (int i = 0 ; i < T_length ; ++i)
    t_data[i] -= rhs;
  return *this;
}

/**
 * @author Randolph Settgast
 * @param[in] rhs value to multiply t_data with
 * @return none
 */
template<int T_length>
inline TensorBaseT< T_length >&
TensorBaseT< T_length >::operator*=( const realT& rhs )
{
  for (int i = 0 ; i < T_length ; ++i)
    t_data[i] *= rhs;
  return *this;
}

template<>
inline TensorBaseT< 3 >&
TensorBaseT<3>::operator*=( const realT& rhs )
{
  t_data[0] *= rhs;
  t_data[1] *= rhs;
  t_data[2] *= rhs;
  return *this;
}

/**
 * @author Randolph Settgast
 * @param[in] rhs value to divide t_data with
 * @return none
 */
template<int T_length>
inline TensorBaseT< T_length >&
TensorBaseT< T_length >::operator/=( const realT& rhs )
{
  const realT irhs = 1 / rhs;
  operator*=( irhs );
  return *this;
}

/**
 * @author Randolph Settgast
 * @param[in] rhs tensor to add
 * @return none
 */
template<int T_length>
inline TensorBaseT< T_length >&
TensorBaseT< T_length >::operator+=( const TensorBaseT< T_length >& rhs )
{
  for (int i = 0 ; i < T_length ; ++i)
    t_data[i] += rhs.t_data[i];
  return *this;
}


template<>
inline TensorBaseT< 3 >&
TensorBaseT<3>::operator+=( const TensorBaseT< 3 >& rhs )
{
  t_data[0] += rhs.t_data[0];
  t_data[1] += rhs.t_data[1];
  t_data[2] += rhs.t_data[2];
  return *this;
}

/**
 * @author Randolph Settgast
 * @param[in] rhs tensor to subract
 * @return none
 */
template<int T_length>
inline TensorBaseT< T_length >&
TensorBaseT< T_length >::operator-=( const TensorBaseT< T_length >& rhs )
{
  for (int i = 0 ; i < T_length ; ++i)
    t_data[i] -= rhs.t_data[i];
  return *this;
}

/**
 * @author Randolph Settgast
 * @param[in] rhs tensor to multiply by
 * @return none
 */
template<int T_length>
inline TensorBaseT< T_length >&
TensorBaseT< T_length >::operator*=( const TensorBaseT< T_length >& rhs )
{
  for (int i = 0 ; i < T_length ; ++i)
    t_data[i] *= rhs.t_data[i];
  return *this;
}

/**
 * @author Randolph Settgast
 * @param[in] rhs tensor to divide by
 * @return none
 */
template<int T_length>
inline TensorBaseT< T_length >&
TensorBaseT< T_length >::operator/=( const TensorBaseT< T_length >& rhs )
{
  for (int i = 0 ; i < T_length ; ++i)
    t_data[i] /= rhs.t_data[i];
  return *this;
}

/*
   template<int N>
   TensorBaseT<N>&& operator+(TensorBaseT<N> &&src1, const TensorBaseT<N> &src2)
   {
   for (int i=0; i<N; ++i)
    src1[i] += src2[i];

   return src1;
   }*/


#endif
