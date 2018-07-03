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

#ifndef ARRAY_3D_T_H_
#define ARRAY_3D_T_H_

//#include "VectorT.h"
#include "VectorT_derived.h"
#include <memory>

// *****************************************************************************
// **** Class Declaration ******************************************************
// *****************************************************************************
template<class TYPE>
class Array3dT : public VectorT<TYPE>
{
public:
  //***** Constructors & Destructors ******************************************
  Array3dT(void);
  Array3dT(const int dim1, const int dim2, const int dim3);
  Array3dT(const Array3dT<TYPE>& source);
  ~Array3dT(void);

  //***** Assignment Operators ************************************************
  inline Array3dT<TYPE>& operator=(const TYPE& rhs);
  inline Array3dT<TYPE>& operator+=(const TYPE& rhs);
  inline Array3dT<TYPE>& operator-=(const TYPE& rhs);
  inline Array3dT<TYPE>& operator*=(const TYPE& rhs);
  inline Array3dT<TYPE>& operator/=(const TYPE& rhs);

  inline Array3dT<TYPE>& operator=(const Array3dT<TYPE>& rhs);
  inline Array3dT<TYPE>& operator+=(const Array3dT<TYPE>& rhs);
  inline Array3dT<TYPE>& operator-=(const Array3dT<TYPE>& rhs);
  inline Array3dT<TYPE>& operator*=(const Array3dT<TYPE>& rhs);
  inline Array3dT<TYPE>& operator/=(const Array3dT<TYPE>& rhs);

  //***** Memory Allocation and Release ***************************************
  void Allocate(const int dim1, const int dim2, const int dim3);

  //***** Accessors **********************************************************

//  inline const TYPE& operator[](const int dim1) const ;
//  inline TYPE& operator[](const int dim1) ;

  inline const TYPE& operator()(const int dim1, const int dim2, const int dim3) const;
  inline TYPE& operator()(const int dim1, const int dim2, const int dim3);

protected:
  int dimension[3];

public:
  //***** Class Information Functions ****************************************
  inline int Dimension(const int dimnum) const;

};

// *****************************************************************************
// **** Class Implementation ***************************************************
// *****************************************************************************

//*****************************************************************************
//***** CONSTRUCTOR/DESTRUCTOR ************************************************
//*****************************************************************************
template<class TYPE>
Array3dT<TYPE>::Array3dT(void)
{
  dimension[0] = 0;
  dimension[1] = 0;
  dimension[2] = 0;
}

template<class TYPE>
Array3dT<TYPE>::Array3dT(const int dim1, const int dim2, const int dim3):
  VectorT<TYPE>::VectorT()
{
  Allocate(dim1, dim2, dim3);
}

template<class TYPE>
Array3dT<TYPE>::Array3dT(const Array3dT<TYPE>& source)
{
  operator=(source);
}

template<class TYPE>
Array3dT<TYPE>& Array3dT<TYPE>::operator=(const Array3dT<TYPE>& rhs)
{
  dimension[0] = rhs.dimension[0];
  dimension[1] = rhs.dimension[1];
  dimension[2] = rhs.dimension[2];
  VectorT<TYPE>::operator=(rhs);
  return (*this);
}

template<class TYPE>
Array3dT<TYPE>::~Array3dT(void)
{}

// *****************************************************************************
// ***** MEMORY ALLOCATION *****************************************************
// *****************************************************************************
template<class TYPE>
void Array3dT<TYPE>::Allocate(const int dim1, const int dim2, const int dim3)
{
  dimension[0] = dim1;
  dimension[1] = dim2;
  dimension[2] = dim3;
  VectorT<TYPE>::resize(dim1 * dim2 * dim3);
}

//*****************************************************************************
//***** ASSIGNMENT OPERATORS **************************************************
//*****************************************************************************
template<class TYPE>
Array3dT<TYPE>& Array3dT<TYPE>::operator=(const TYPE& rhs)
{
  VectorT<TYPE>::operator=(rhs);
  return (*this);
}

template<class TYPE>
Array3dT<TYPE>& Array3dT<TYPE>::operator+=(const TYPE& rhs)
{
  VectorT<TYPE>::operator+=(rhs);
  return (*this);
}

template<class TYPE>
Array3dT<TYPE>& Array3dT<TYPE>::operator-=(const TYPE& rhs)
{
  VectorT<TYPE>::operator-=(rhs);
  return (*this);
}

template<class TYPE>
Array3dT<TYPE>& Array3dT<TYPE>::operator*=(const TYPE& rhs)
{
  VectorT<TYPE>::operator*=(rhs);
  return (*this);
}

template<class TYPE>
Array3dT<TYPE>& Array3dT<TYPE>::operator/=(const TYPE& rhs)
{
  VectorT<TYPE>::operator/=(rhs);
  return (*this);
}

template<class TYPE>
Array3dT<TYPE>& Array3dT<TYPE>::operator+=(const Array3dT<TYPE>& rhs)
{
#if RANGE_CHECKING==1
  if ( dimension[0] != rhs.dimension[0] ||
       dimension[1] != rhs.dimension[1] ||
       dimension[2] != rhs.dimension[2] )
    throw 1;
#endif
  VectorT<TYPE>::operator+=(rhs);
  return (*this);
}

template<class TYPE>
Array3dT<TYPE>& Array3dT<TYPE>::operator-=(const Array3dT<TYPE>& rhs)
{
#if RANGE_CHECKING==1
  if ( dimension[0] != rhs.dimension[0] ||
       dimension[1] != rhs.dimension[1] ||
       dimension[2] != rhs.dimension[2] )
    throw 1;
#endif
  VectorT<TYPE>::operator-=(rhs);
  return (*this);
}

template<class TYPE>
Array3dT<TYPE>& Array3dT<TYPE>::operator*=(const Array3dT<TYPE>& rhs)
{
#if RANGE_CHECKING==1
  if ( dimension[0] != rhs.dimension[0] ||
       dimension[1] != rhs.dimension[1] ||
       dimension[2] != rhs.dimension[2] )
    throw 1;
#endif
  VectorT<TYPE>::operator*=(rhs);
  return (*this);
}

template<class TYPE>
Array3dT<TYPE>& Array3dT<TYPE>::operator/=(const Array3dT<TYPE>& rhs)
{
#if RANGE_CHECKING==1
  if ( dimension[0] != rhs.dimension[0] ||
       dimension[1] != rhs.dimension[1] ||
       dimension[2] != rhs.dimension[2] )
    throw 1;
#endif
  VectorT<TYPE>::operator/=(rhs);
  return (*this);
}

//*****************************************************************************
//***** ELEMENT ACCESS ********************************************************
//*****************************************************************************
/*
   template <class TYPE>
   const TYPE& Array3dT<TYPE>::operator[](const int index) const
   {
   throw 1;
   return VectorT<TYPE>::operator[](index );
   }

   template <class TYPE>
   TYPE& Array3dT<TYPE>::operator[](const int index)
   {
   throw 1;
   return VectorT<TYPE>::operator[](index );
   }
 */
template<class TYPE>
const TYPE& Array3dT<TYPE>::operator()(const int dim1, const int dim2, const int dim3) const
{
#if RANGE_CHECKING==1
  if ((dim1) < 0 || (dim1 ) >= dimension[0] ||
      (dim2) < 0 || (dim2 ) >= dimension[1] ||
      (dim3) < 0 || (dim3 ) >= dimension[2] )
    throw 1;
#endif

  return VectorT<TYPE>::operator[](
    (dim1) * dimension[1] * dimension[2] + (dim2) * dimension[2] + (dim3));
}

template<class TYPE>
TYPE& Array3dT<TYPE>::operator()(const int dim1, const int dim2, const int dim3)
{
#if RANGE_CHECKING==1
  if ((dim1) < 0 || (dim1 ) >= dimension[0] ||
      (dim2) < 0 || (dim2 ) >= dimension[1] ||
      (dim3) < 0 || (dim3 ) >= dimension[2] )
    throw 1;
#endif

  return VectorT<TYPE>::operator[](
    (dim1) * dimension[1] * dimension[2] + (dim2) * dimension[2] + (dim3));
}

template<class TYPE>
int Array3dT<TYPE>::Dimension(const int dimnum) const
{
#if RANGE_CHECKING==1
  if ( dimnum>2 )
    throw 1;
#endif
  return dimension[dimnum];
}

#endif
