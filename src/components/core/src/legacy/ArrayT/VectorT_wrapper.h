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
 * @brief This file contains the definition of the VectorT class
 * @file VectorT.h
 * @author Randolph Settgast
 */


#ifndef _VECTOR_T_H_
#define _VECTOR_T_H_

#include <vector>
#include <string>

/**
 * @brief VectorT is a wrapper for std::vector.
 * @author Randolph Settgast
 * @tparam TYPE type of data that is contained.
 *
 * VectorT is a class to add some operators to std::vector. It has a std::vector
 * as a member, and
 * provides some wrapper function as well.
 */
template<class TYPE>
class VectorT
{
private:
  /// data_vector is is the actual data container.
  std::vector<TYPE> data_vector;
public:
//***** Constructors & Destructors ********************************************
  /// default constructor
  VectorT(void);

  /// constructor that initializes the size
  VectorT( const int num_elem );

  /// constructor that initializes the name
  VectorT( const std::string& name );

  /// constructor that initializes size and name
  VectorT( const int num_elem, const std::string& name );

  /// copy constructor
  VectorT( const VectorT& source );

  /// default destructor
  virtual ~VectorT(void);


// ***** Memory Allocation ****************************************************
  /// wrapper for std::vector.resize()
  inline void resize(const int num_elem)    { data_vector.resize(num_elem);  }

  /// wrapper for std::vector.clear()
  inline void clear(void)                   { data_vector.clear();  }

  /// wrapper for std::vector.push_back()
  inline void push_back( const TYPE& val )  { data_vector.push_back(val);  }

  /// wrapper for std::vector.size()
  inline int size(void) const               { return data_vector.size();  }


  //***** Assignment Operators ************************************************
  /// equals operator for other vectors of any type
  template<class rTYPE> VectorT& operator=( const VectorT<rTYPE>& rhs );

  /// plus equals operator for other vectors of any type
  template<class rTYPE> VectorT& operator+=( const VectorT<rTYPE>& rhs );

  /// minus equals operator for other vectors of any type
  template<class rTYPE> VectorT& operator-=( const VectorT<rTYPE>& rhs );

  /// multiply equals operator for other vectors of any type
  template<class rTYPE> VectorT& operator*=( const VectorT<rTYPE>& rhs );

  /// divide equals operator for other vectors of any type
  template<class rTYPE> VectorT& operator/=( const VectorT<rTYPE>& rhs );

  /// equals operator for individual value of any type
  template<class rTYPE> VectorT& operator=( const rTYPE& rhs );

  /// plus equals operator for individual value of any type
  template<class rTYPE> VectorT& operator+=( const rTYPE& rhs );

  /// minus equals operator for individual value of any type
  template<class rTYPE> VectorT& operator-=( const rTYPE& rhs );

  /// multiply equals operator for individual value of any type
  template<class rTYPE> VectorT& operator*=( const rTYPE& rhs );

  /// divide equals operator for individual value of any type
  template<class rTYPE> VectorT& operator/=( const rTYPE& rhs );


//***** Range Checking ********************************************************
#if RANGE_CHECKING==1
  /// wrapper for std::vector.at()
  inline TYPE& operator[](const int index)  { return (data_vector.at(index));}
  /// const wrapper for std::vector.at()
  inline const TYPE& operator[](const int index) const { return (data_vector.at(index));}
#else
  /// wrapper for std::vector.operator[]()
  inline TYPE& operator[](const int index)  { return (data_vector[index]);}
  /// wrapper for std::vector.operator[]()
  inline const TYPE& operator[](const int index) const { return (data_vector[index]);}
#endif


//***** Access Pointers *******************************************************
  /// returns a TYPE* to the first element of the std::vector
  inline TYPE* Pointer(void)
  {
    TYPE* rval=0;
    if( this->data_vector.size() > 0 )
      rval = (&data_vector[0]);
    return rval;
  }

  /// returns a const TYPE* to the first element of the std::vector
  inline const TYPE* Pointer(void) const
  {
    if( this->data_vector.size() > 0 )
      return (&data_vector[0]);
    else
      return 0;
  }

};

// *****************************************************************************
// **** Class Implementation ***************************************************
// *****************************************************************************


//*****************************************************************************
//***** CONSTRUCTOR/DESTRUCTOR ************************************************
//*****************************************************************************


template<class TYPE>
VectorT<TYPE>::VectorT(void):
  data_vector()
{}


template<class TYPE>
VectorT<TYPE>::VectorT(const int num_elem):
  data_vector(num_elem)
{}

template<class TYPE>
VectorT<TYPE>::VectorT( const std::string& name ):
  data_vector()
{}

template<class TYPE>
VectorT<TYPE>::VectorT(const int num_elem, const std::string& name):
  data_vector(num_elem)
{}



template<class TYPE>
VectorT<TYPE>::~VectorT(void)
{}


template<class TYPE>
VectorT<TYPE>::VectorT( const VectorT& source ):
  data_vector()
{
  operator=(source);
}



//*****************************************************************************
//***** ASSIGNMENT OPERATORS **************************************************
//*****************************************************************************

/**
 * @author Randolph Settgast
 * @tparam rTYPE type of rhs
 * @param[in] rhs VectorT to set *this equal to
 * @return reference to *this
 *
 * Templated operator=() to set the values contained in *this to those of rhs on
 * an element
 * by element basis. If there is no valid cast from rTYPE to TYPE, then
 * compilation should fail.
 */
template<class TYPE>
template<class rTYPE>
inline VectorT<TYPE>& VectorT<TYPE>::operator=( const VectorT<rTYPE>& rhs )
{
  if( this->size() == rhs.size() && this->size() )
  {
    TYPE* ptr = &data_vector[0];
    const rTYPE* rhs_ptr = &(rhs[0]);
    for( unsigned int a=0 ; a<data_vector.size() ; ++a )
      *(ptr++) = *(rhs_ptr++);
  }

  return (*this);
}

/**
 * @author Randolph Settgast
 * @tparam rTYPE type of rhs
 * @param[in] rhs VectorT to set *this equal to
 * @return reference to *this
 *
 * Templated operator+=() to add the values contained in rhs to those of *this
 * on an element
 * by element basis. If there is no valid cast from rTYPE to TYPE, then
 * compilation should fail.
 */
template<class TYPE>
template<class rTYPE>
inline VectorT<TYPE>& VectorT<TYPE>::operator+=( const VectorT<rTYPE>& rhs )
{
  if( this->size() == rhs.size() && this->size() )
  {
    TYPE* ptr = &data_vector[0];
    const rTYPE* rhs_ptr = &(rhs[0]);
    for( unsigned int a=0 ; a<data_vector.size() ; ++a )
      *(ptr++) += *(rhs_ptr++);
  }
  else
    throw 1; //1;//eOutOfRange;
  return (*this);
}


/**
 * @author Randolph Settgast
 * @tparam rTYPE type of rhs
 * @param[in] rhs VectorT to set *this equal to
 * @return reference to *this
 *
 * Templated operator-=() to subtract the values contained in rhs from those of
 **this on an element
 * by element basis. If there is no valid cast from rTYPE to TYPE, then
 * compilation should fail.
 */
template<class TYPE>
template<class rTYPE>
inline VectorT<TYPE>& VectorT<TYPE>::operator-=( const VectorT<rTYPE>& rhs )
{
  if( this->size() == rhs.size() && this->size() )
  {
    TYPE* ptr = &data_vector[0];
    const rTYPE* rhs_ptr = &(rhs[0]);
    for( unsigned int a=0 ; a<(*this).size() ; ++a )
      *(ptr++) -= *(rhs_ptr++);
  }
  else
    throw 1; //eOutOfRange;
  return (*this);
}

/**
 * @author Randolph Settgast
 * @tparam rTYPE type of rhs
 * @param[in] rhs VectorT to set *this equal to
 * @return reference to *this
 *
 * Templated operator*=() to multiply the values contained in *this by those of
 * rhs on an element
 * by element basis. If there is no valid cast from rTYPE to TYPE, then
 * compilation should fail.
 */
template<class TYPE>
template<class rTYPE>
inline VectorT<TYPE>& VectorT<TYPE>::operator*=( const VectorT<rTYPE>& rhs )
{
  if( this->size() == rhs.size() && this->size() )
  {
    TYPE* ptr = &data_vector[0];
    const rTYPE* rhs_ptr = &(rhs[0]);
    for( unsigned int a=0 ; a<(*this).size() ; ++a )
      *(ptr++) *= *(rhs_ptr++);
  }
  else
    throw 1; //eOutOfRange;
  return (*this);
}


/**
 * @author Randolph Settgast
 * @tparam rTYPE type of rhs
 * @param[in] rhs VectorT to set *this equal to
 * @return reference to *this
 *
 * Templated operator/=() to divide the values contained in *this by those of
 * rhs on an element
 * by element basis. If there is no valid cast from rTYPE to TYPE, then
 * compilation should fail.
 */
template<class TYPE>
template<class rTYPE>
inline VectorT<TYPE>& VectorT<TYPE>::operator/=( const VectorT<rTYPE>& rhs )
{
  if( this->size() == rhs.size() && this->size() )
  {
    TYPE* ptr = &data_vector[0];
    const rTYPE* rhs_ptr = &(rhs[0]);
    for( unsigned int a=0 ; a<data_vector.size() ; ++a )
      *(ptr++) /= *(rhs_ptr++);
  }
  else
    throw 1; //eOutOfRange;
  return (*this);
}


/**
 * @author Randolph Settgast
 * @tparam rTYPE type of rhs
 * @param[in] rhs VectorT to set *this equal to
 * @return reference to *this
 *
 * Templated operator=() to set the values contained in *this to that of rhs.
 * If there is no valid cast from rTYPE to TYPE, then compilation should fail.
 */
template<class TYPE>
template<class rTYPE>
inline VectorT<TYPE>& VectorT<TYPE>::operator=( const rTYPE& rhs )
{
  if( data_vector.size() )
  {
    TYPE* ptr = &data_vector[0];
    for( unsigned int a=0 ; a<data_vector.size() ; ++a )
      *(ptr++) = rhs;
  }
  return (*this);
}

/**
 * @author Randolph Settgast
 * @tparam rTYPE type of rhs
 * @param[in] rhs VectorT to set *this equal to
 * @return reference to *this
 *
 * Templated operator+=() to add rhs the values contained in *this.
 * If there is no valid cast from rTYPE to TYPE, then compilation should fail.
 */
template<class TYPE>
template<class rTYPE>
inline VectorT<TYPE>& VectorT<TYPE>::operator+=( const rTYPE& rhs )
{
  if( data_vector.size() )
  {
    TYPE* ptr = &data_vector[0];
    for( unsigned int a=0 ; a<data_vector.size() ; ++a )
      *(ptr++) += rhs;
  }
  return (*this);
}

/**
 * @author Randolph Settgast
 * @tparam rTYPE type of rhs
 * @param[in] rhs VectorT to set *this equal to
 * @return reference to *this
 *
 * Templated operator-=() to subtract rhs from the values contained in *this.
 * If there is no valid cast from rTYPE to TYPE, then compilation should fail.
 */
template<class TYPE>
template<class rTYPE>
inline VectorT<TYPE>& VectorT<TYPE>::operator-=( const rTYPE& rhs )
{
  if( data_vector.size() )
  {
    TYPE* ptr = &data_vector[0];
    for( unsigned int a=0 ; a<data_vector.size() ; ++a )
      *(ptr++) -= rhs;
  }
  return (*this);
}

/**
 * @author Randolph Settgast
 * @tparam rTYPE type of rhs
 * @param[in] rhs VectorT to set *this equal to
 * @return reference to *this
 *
 * Templated operator*=() to set multiply the values contained in *this by rhs.
 * If there is no valid cast from rTYPE to TYPE, then compilation should fail.
 */
template<class TYPE>
template<class rTYPE>
inline VectorT<TYPE>& VectorT<TYPE>::operator*=( const rTYPE& rhs )
{
  if( data_vector.size() )
  {
    TYPE* ptr = &data_vector[0];
    for( unsigned int a=0 ; a<data_vector.size() ; ++a )
      *(ptr++) *= rhs;
  }
  return (*this);
}

/**
 * @author Randolph Settgast
 * @tparam rTYPE type of rhs
 * @param[in] rhs VectorT to set *this equal to
 * @return reference to *this
 *
 * Templated operator*=() to set divide the values contained in *this by rhs.
 * If there is no valid cast from rTYPE to TYPE, then compilation should fail.
 */
template<class TYPE>
template<class rTYPE>
inline VectorT<TYPE>& VectorT<TYPE>::operator/=( const rTYPE& rhs )
{
  if( data_vector.size() )
  {
    TYPE* ptr = &data_vector[0];
    for( unsigned int a=0 ; a<data_vector.size() ; ++a )
      *(ptr++) /= rhs;
  }
  return (*this);
}



#endif
