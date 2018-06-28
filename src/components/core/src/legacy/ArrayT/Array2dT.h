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

#ifndef ARRAY_2D_T_H_
#define ARRAY_2D_T_H_

//#include "VectorT.h"
#include "VectorT_derived.h"

#include <memory>


// *****************************************************************************
// **** Class Declaration ******************************************************
// *****************************************************************************
template <typename TYPE>
class Array2dT : public VectorT<TYPE>
{
public:

  typedef long size_type;

  //***** Constructors & Destructors ******************************************
  Array2dT(void);
  Array2dT( const long dim1, const long dim2 );
  Array2dT(const Array2dT<TYPE>& source);


  //***** Assignment Operators ************************************************
  template <class rhsTYPE> Array2dT<TYPE>& operator=( const Array2dT<rhsTYPE>& rhs );
  template <class rhsTYPE> Array2dT<TYPE>& operator+=( const Array2dT<rhsTYPE>& rhs );
  template <class rhsTYPE> Array2dT<TYPE>& operator-=( const Array2dT<rhsTYPE>& rhs );
  template <class rhsTYPE> Array2dT<TYPE>& operator*=( const Array2dT<rhsTYPE>& rhs );
  template <class rhsTYPE> Array2dT<TYPE>& operator/=( const Array2dT<rhsTYPE>& rhs );


  template <class rhsTYPE> Array2dT<TYPE>& operator=( const rhsTYPE& rhs );
  template <class rhsTYPE> Array2dT<TYPE>& operator+=( const rhsTYPE& rhs );
  template <class rhsTYPE> Array2dT<TYPE>& operator-=( const rhsTYPE& rhs );
  template <class rhsTYPE> Array2dT<TYPE>& operator*=( const rhsTYPE& rhs );
  template <class rhsTYPE> Array2dT<TYPE>& operator/=( const rhsTYPE& rhs );



  //***** Memory Allocation and Release ***************************************
	void setDimensions(const size_type numDims, const size_type * dims);
  void resize2( const size_type dim1 , const size_type dim2 );
  void resize(const size_type num_elem)
  { resize2( num_elem, dimension[1] );  }

  void Insert( const size_type dim1, const TYPE& t);
  void Erase( const size_type dim1 );


  //***** Accessors **********************************************************



  inline const TYPE& operator()(const size_type dim1, const size_type dim2) const;
  inline TYPE& operator()(const size_type dim1, const size_type dim2);



private:
  size_type dimension[2];

public:
  //***** Class Information Functions ****************************************

  inline int numDimensions() const
  { return 2; }

  inline size_type Dimension( const size_type dimnum ) const ;

  inline const TYPE* operator[]( const size_type index ) const
  {
    return ( &(VectorT<TYPE>::operator[](index*dimension[1])) );
  }

  inline TYPE* operator[]( const size_type index )
  {
    return ( &(VectorT<TYPE>::operator[](index*dimension[1])) );
  }

};


// *****************************************************************************
// **** Class Implementation ***************************************************
// *****************************************************************************


//*****************************************************************************
//***** CONSTRUCTOR/DESTRUCTOR ************************************************
//*****************************************************************************
template <class TYPE>
Array2dT<TYPE>::Array2dT(void)
{
  dimension[0] = 0;
  dimension[1] = 0;
}

template <class TYPE>
Array2dT<TYPE>::Array2dT( const size_type dim1, const size_type dim2):
  VectorT<TYPE>::VectorT()
{
  resize2(dim1,dim2);
}

template <class TYPE>
Array2dT<TYPE>::Array2dT(const Array2dT<TYPE>& source):
  VectorT<TYPE>(source)
{
  operator=(source);
}



// *****************************************************************************
// ***** MEMORY ALLOCATION *****************************************************
// *****************************************************************************
template <class TYPE>
void Array2dT<TYPE>::setDimensions( long numDims, const long * dims)
{
  if (numDims != 2) throw 1;
  dimension[0] = dims[0];
  dimension[1] = dims[1];

  std::vector<TYPE>::resize(dimension[0] * dimension[1]);
}

template <class TYPE> 
void Array2dT<TYPE>::resize2( const size_type dim1 , const size_type dim2)
{
  dimension[0] = dim1;
  dimension[1] = dim2;

  std::vector<TYPE>::resize(dimension[0] * dimension[1]);
}

template <class TYPE>
void Array2dT<TYPE>::Insert( const size_type dim1, const TYPE& t)
{
  const size_type index = (dim1)*dimension[1];
  for(size_type tt = 0 ; tt < dimension[1] ; tt++)
    this->insert(this->begin() + index, t);
  ++dimension[0];
}

template <class TYPE>
void Array2dT<TYPE>::Erase( const size_type dim1 )
{
  const size_type index = (dim1)*dimension[1];
  for(size_type tt = 0 ; tt < dimension[1] ; tt++)
    this->erase(this->begin() + index );
  --dimension[0];
}

//*****************************************************************************
//***** ASSIGNMENT OPERATORS **************************************************
//*****************************************************************************

template <class TYPE>
template <class rhsTYPE>
Array2dT<TYPE>& Array2dT<TYPE>::operator=( const Array2dT<rhsTYPE>& rhs )
{
  dimension[0] = rhs.dimension[0];
  dimension[1] = rhs.dimension[1];
  VectorT<TYPE>::operator=( static_cast<VectorT<rhsTYPE> >(rhs) );
  return (*this);
}

template <class TYPE>
template <class rhsTYPE>
Array2dT<TYPE>& Array2dT<TYPE>::operator+=( const Array2dT<rhsTYPE>& rhs )
{
#if RANGE_CHECKING==1
  if ( dimension[0] != rhs.dimension[0] ||
       dimension[1] != rhs.dimension[1])
    throw 1;
#endif
  VectorT<TYPE>::operator+=( static_cast<VectorT<rhsTYPE> >(rhs) );
  return (*this);
}

template <class TYPE>
template <class rhsTYPE>
Array2dT<TYPE>& Array2dT<TYPE>::operator-=( const Array2dT<rhsTYPE>& rhs )
{
#if RANGE_CHECKING==1
  if ( dimension[0] != rhs.dimension[0] ||
       dimension[1] != rhs.dimension[1])
    throw 1;
#endif
  VectorT<TYPE>::operator-=( static_cast<VectorT<rhsTYPE> >(rhs) );
  return (*this);
}

template <class TYPE>
template <class rhsTYPE>
Array2dT<TYPE>& Array2dT<TYPE>::operator*=( const Array2dT<rhsTYPE>& rhs )
{
#if RANGE_CHECKING==1
  if ( dimension[0] != rhs.dimension[0] ||
       dimension[1] != rhs.dimension[1])
    throw 1;
#endif
  VectorT<TYPE>::operator*=( static_cast<VectorT<rhsTYPE> >(rhs) );
  return (*this);
}

template <class TYPE>
template <class rhsTYPE>
Array2dT<TYPE>& Array2dT<TYPE>::operator/=( const Array2dT<rhsTYPE>& rhs )
{
#if RANGE_CHECKING==1
  if ( dimension[0] != rhs.dimension[0] ||
       dimension[1] != rhs.dimension[1])
    throw 1;
#endif
  VectorT<TYPE>::operator/=( static_cast<VectorT<rhsTYPE> >(rhs) );
  return (*this);
}



template <class TYPE>
template <class rhsTYPE>
Array2dT<TYPE>& Array2dT<TYPE>::operator=( const rhsTYPE& rhs )
{
  VectorT<TYPE>::operator=(rhs);
  return (*this);
}


template <class TYPE>
template <class rhsTYPE>
Array2dT<TYPE>& Array2dT<TYPE>::operator+=( const rhsTYPE& rhs )
{
  VectorT<TYPE>::operator+=(rhs);
  return (*this);
}

template <class TYPE>
template <class rhsTYPE>
Array2dT<TYPE>& Array2dT<TYPE>::operator-=( const rhsTYPE& rhs )
{
  VectorT<TYPE>::operator-=(rhs);
  return (*this);
}

template <class TYPE>
template <class rhsTYPE>
Array2dT<TYPE>& Array2dT<TYPE>::operator*=( const rhsTYPE& rhs )
{
  VectorT<TYPE>::operator*=(rhs);
  return (*this);
}

template <class TYPE>
template <class rhsTYPE>
Array2dT<TYPE>& Array2dT<TYPE>::operator/=( const rhsTYPE& rhs )
{
  VectorT<TYPE>::operator/=(rhs);
  return (*this);
}



//*****************************************************************************
//***** ELEMENT ACCESS ********************************************************
//*****************************************************************************

template <class TYPE>
const TYPE& Array2dT<TYPE>::operator()( const size_type dim1, const size_type dim2 ) const
{
#if RANGE_CHECKING==1
  if (/*(dim1) < 0 ||*/ (dim1 ) >= dimension[0] ||
                        /*(dim2) < 0 ||*/ (dim2 ) >= dimension[1] )
    throw 1;
#endif

  return VectorT<TYPE>::operator[]((dim1)*dimension[1] + (dim2));
//  return (p_data[dim1])[dim2];
}

template <class TYPE>
TYPE& Array2dT<TYPE>::operator()( const size_type dim1, const size_type dim2 )
{
#if RANGE_CHECKING==1
  if (/*(dim1) < 0 ||*/ (dim1 ) >= dimension[0] ||
                        /*(dim2) < 0 ||*/ (dim2 ) >= dimension[1] )
    throw 1;

#endif

  return VectorT<TYPE>::operator[]((dim1)*dimension[1] + (dim2));
//  return (p_data[dim1])[dim2];
}

template <class TYPE>
typename Array2dT<TYPE>::size_type Array2dT<TYPE>::Dimension( const size_type dimnum ) const
{
#if RANGE_CHECKING==1
  if ( dimnum>1 )
    throw 1;
#endif
  return dimension[dimnum];
}

#endif
