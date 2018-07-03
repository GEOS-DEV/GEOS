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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @brief This file contains the definition of the Array1dT class
 * @file Array1dT.h
 * @author Randolph Settgast
 */

#ifndef ARRAY_1D_T_H_
#define ARRAY_1D_T_H_

//#include "VectorT.h"
#include "VectorT_derived.h"

template<typename TYPE> class Array2dT;

/**
 * @brief Array1dT is a derived class of VectorT.
 * @author Randolph Settgast
 * @tparam TYPE type of data that is contained.
 *
 * Array1dT is one-dimensional array that is a derived class of VectorT.
 * Array1dT provides some
 * legacy function wrappers and accessors.
 *
 */
template <typename TYPE>
class Array1dT : public VectorT<TYPE>
{
public:

  typedef TYPE* pointer;
  typedef long int size_type;

  //***** Constructors & Destructors ******************************************
  /// default constructor
  Array1dT(void);

  /// constructor that initializes the size
  Array1dT( const size_type num_elem );

  /// repedative sequence constructor - creates a vector with num_elem copies of
  // value
  Array1dT(const size_type num_elem, const TYPE value);

  /// copy constructor
  Array1dT(const Array1dT<TYPE>& source);

  /// default destructor
  virtual ~Array1dT(void);

  Array1dT& operator=( const Array1dT& rhs );

  //***** Assignment Operators ************************************************
  /// equals operator for other arrays of any type
  template <class rhsTYPE> Array1dT<TYPE>& operator=( const Array1dT<rhsTYPE>& rhs );

  template <class rhsTYPE> Array1dT<TYPE>& operator=( const std::vector<rhsTYPE>& rhs );

  /// += operator for other arrays of any type
  template <class rhsTYPE> Array1dT<TYPE>& operator+=( const Array1dT<rhsTYPE>& rhs );

  /// -= operator for other arrays of any type
  template <class rhsTYPE> Array1dT<TYPE>& operator-=( const Array1dT<rhsTYPE>& rhs );

  /// *= operator for other arrays of any type
  template <class rhsTYPE> Array1dT<TYPE>& operator*=( const Array1dT<rhsTYPE>& rhs );

  /// /= operator for other arrays of any type
  template <class rhsTYPE> Array1dT<TYPE>& operator/=( const Array1dT<rhsTYPE>& rhs );

  /// equals operator for individual value of any type
  template <class rhsTYPE> Array1dT<TYPE>& operator=( const rhsTYPE& rhs );

  /// += operator for individual value of any type
  template <class rhsTYPE> Array1dT<TYPE>& operator+=( const rhsTYPE& rhs );

  /// -= operator for individual value of any type
  template <class rhsTYPE> Array1dT<TYPE>& operator-=( const rhsTYPE& rhs );

  /// *= operator for individual value of any type
  template <class rhsTYPE> Array1dT<TYPE>& operator*=( const rhsTYPE& rhs );

  /// /= operator for individual value of any type
  template <class rhsTYPE> Array1dT<TYPE>& operator/=( const rhsTYPE& rhs );

  //***** Accessors **********************************************************
  /// const square bracket accessor
  const TYPE& operator[](const size_type dim1) const;
  /// square bracket accessor
  TYPE& operator[](const size_type dim1);

  /// const parentheses accessor
  const TYPE& operator()(const size_type dim1) const;
  /// parentheses accessor
  TYPE& operator()(const size_type dim1);

  const Array1dT& Slice( Array2dT<TYPE>& array2d, const size_type index ) const;


protected:

public:
  //***** Class Information Functions ****************************************


};


// *****************************************************************************
// **** Class Implementation ***************************************************
// *****************************************************************************

#include "Array2dT.h"


template <class TYPE>
const Array1dT<TYPE>& Array1dT<TYPE>::Slice( Array2dT<TYPE>& array2d, const size_type index ) const
{
  this->VectorT<TYPE>::SetValue( array2d.FirstIndexPointer(index), array2d.dimension(1) );
  return *this;
}

//*****************************************************************************
//***** CONSTRUCTOR/DESTRUCTOR ************************************************
//*****************************************************************************
template <class TYPE>
Array1dT<TYPE>::Array1dT(void)
{}

template <class TYPE>
Array1dT<TYPE>::Array1dT(const size_type num_elem):
  VectorT<TYPE>(num_elem)
{}

template <class TYPE>
Array1dT<TYPE>::Array1dT(const size_type num_elem, TYPE value):
  VectorT<TYPE>(num_elem,value)
{}

template <class TYPE>
Array1dT<TYPE>::Array1dT(const Array1dT<TYPE>& source):
  VectorT<TYPE>(source)
{
  operator=(source);
}

template <class TYPE>
Array1dT<TYPE>::~Array1dT(void)
{}

//*****************************************************************************
//***** ASSIGNMENT OPERATORS **************************************************
//*****************************************************************************
template <class TYPE>
template <class rhsTYPE>
Array1dT<TYPE>& Array1dT<TYPE>::operator=( const Array1dT<rhsTYPE>& rhs )
{
  VectorT<TYPE>::operator=( static_cast<const VectorT<rhsTYPE>& >(rhs) );
  return (*this);
}


template <class TYPE>
inline Array1dT<TYPE>& Array1dT<TYPE>::operator=( const Array1dT<TYPE>& rhs )
{
  VectorT<TYPE>::operator=( static_cast<const VectorT<TYPE>& >(rhs) );
  return (*this);
}


template <class TYPE>
template <class rhsTYPE>
Array1dT<TYPE>& Array1dT<TYPE>::operator=( const std::vector<rhsTYPE>& rhs )
{
  VectorT<TYPE>::operator=( rhs );
  return (*this);
}


template <class TYPE>
template <class rhsTYPE>
Array1dT<TYPE>& Array1dT<TYPE>::operator+=( const Array1dT<rhsTYPE>& rhs )
{
  VectorT<TYPE>::operator+=( static_cast<const VectorT<rhsTYPE>& >(rhs) );
  return (*this);
}

template <class TYPE>
template <class rhsTYPE>
Array1dT<TYPE>& Array1dT<TYPE>::operator-=( const Array1dT<rhsTYPE>& rhs )
{
  VectorT<TYPE>::operator-=( static_cast<const VectorT<rhsTYPE>& >(rhs) );
  return (*this);
}

template <class TYPE>
template <class rhsTYPE>
Array1dT<TYPE>& Array1dT<TYPE>::operator*=( const Array1dT<rhsTYPE>& rhs )
{
  VectorT<TYPE>::operator*=( static_cast<const VectorT<rhsTYPE>& >(rhs) );
  return (*this);
}

template <class TYPE>
template <class rhsTYPE>
Array1dT<TYPE>& Array1dT<TYPE>::operator/=( const Array1dT<rhsTYPE>& rhs )
{
  VectorT<TYPE>::operator/=( static_cast<const VectorT<rhsTYPE>& >(rhs) );
  return (*this);
}



template <class TYPE>
template <class rhsTYPE>
Array1dT<TYPE>& Array1dT<TYPE>::operator=( const rhsTYPE& rhs )
{
  VectorT<TYPE>::operator=(rhs);
  return (*this);
}

template <class TYPE>
template <class rhsTYPE>
Array1dT<TYPE>& Array1dT<TYPE>::operator+=( const rhsTYPE& rhs )
{
  VectorT<TYPE>::operator+=(rhs);
  return (*this);
}

template <class TYPE>
template <class rhsTYPE>
Array1dT<TYPE>& Array1dT<TYPE>::operator-=( const rhsTYPE& rhs )
{
  VectorT<TYPE>::operator-=(rhs);
  return (*this);
}

template <class TYPE>
template <class rhsTYPE>
Array1dT<TYPE>& Array1dT<TYPE>::operator*=( const rhsTYPE& rhs )
{
  VectorT<TYPE>::operator*=(rhs);
  return (*this);
}

template <class TYPE>
template <class rhsTYPE>
Array1dT<TYPE>& Array1dT<TYPE>::operator/=( const rhsTYPE& rhs )
{
  VectorT<TYPE>::operator/=(rhs);
  return (*this);
}



//*****************************************************************************
//***** ELEMENT ACCESS ********************************************************
//*****************************************************************************

template <class TYPE>
inline const TYPE& Array1dT<TYPE>::operator[](const size_type index) const
{
  return VectorT<TYPE>::operator[](index );
}

template <class TYPE>
inline TYPE& Array1dT<TYPE>::operator[](const size_type index)
{
  return VectorT<TYPE>::operator[](index );
}

template <class TYPE>
inline const TYPE& Array1dT<TYPE>::operator()(const size_type index) const
{
  return VectorT<TYPE>::operator[](index );
}

template <class TYPE>
inline TYPE& Array1dT<TYPE>::operator()(const size_type index)
{
  return VectorT<TYPE>::operator[](index );
}



#endif
