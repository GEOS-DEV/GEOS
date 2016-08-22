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
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _ARRAY_3D_T_H_
#define _ARRAY_3D_T_H_

//#include "VectorT.h"
#include "VectorT_derived.h"
#include <memory>

// *****************************************************************************
// **** Class Declaration ******************************************************
// *****************************************************************************
template<class TYPE>
class Array3dT: public VectorT<TYPE>
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
Array3dT<TYPE>::Array3dT(const int dim1, const int dim2, const int dim3) :
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
{
}

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
  if ( dimnum>2 ) throw 1;
#endif
  return dimension[dimnum];
}

#endif 
