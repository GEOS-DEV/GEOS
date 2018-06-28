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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
// SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S.
// Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National
// Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy,
// completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe
// privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or
// services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or
// favoring by the United States Government or Lawrence Livermore National
// Security,
//     LLC. The views and opinions of authors expressed herein do not
// necessarily state or reflect those of the United States Government or
// Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising
// or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The
// BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
