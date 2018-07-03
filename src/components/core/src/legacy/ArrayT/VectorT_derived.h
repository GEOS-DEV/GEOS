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
 * @brief This file contains the definition of the VectorT class
 * @file VectorT.h
 * @author Randolph Settgast
 */


#ifndef VECTOR_T_H_
#define VECTOR_T_H_

#include <vector>


#ifndef RANGE_CHECKING
#define RANGE_CHECKING 0
#endif

/**
 * @brief VectorT is a wrapper for std::vector.
 * @author Randolph Settgast
 * @tparam TYPE type of data that is contained.
 *
 * VectorT is a class to add some operators to std::vector. It derives from
 * std::vector, which is a big "no-no" if
 * you plan on using VectorT as a derived class....so don't EVER use a
 * std::vector* to allocate a VectorT object!
 */
template<typename TYPE>
class VectorT : public std::vector<TYPE>
{
public:
  typedef long size_type;

  //***** Constructors & Destructors
  // ********************************************
  /// default constructor
  VectorT(void):
    std::vector<TYPE>() {}

  /// constructor that initializes the size
  VectorT( const size_type num_elem ):
    std::vector<TYPE>(num_elem) {}

  /// copy constructor
  VectorT( const VectorT& source ):
    std::vector<TYPE>( static_cast< std::vector<TYPE> >(source)) {}

  /// repetitive sequence constructor creates a vector with num_elem copies of
  // value
  VectorT(const size_type num_elem, const TYPE value):
    std::vector<TYPE>(num_elem,value) {}

  /// default destructor
  virtual ~VectorT(void) {}



  //***** Assignment Operators ************************************************
  /// equals operator that sets *this to vectors of any type
  template<class rTYPE> VectorT& operator=( const VectorT<rTYPE>& rhs )
  {
    this->std::vector<TYPE>::assign(rhs.begin(), rhs.end());
    return (*this);
  }

  /// equals operator that sets *this to vectors of any type
  VectorT& operator=( const VectorT& rhs )
  {
    this->std::vector<TYPE>::operator=( rhs );
    return (*this);
  }

  /// equals operator that sets *this to std::vectors of the same type
  /// Fixme check if OK
  VectorT& operator=( const std::vector<TYPE>& rhs )
  {
    this->std::vector<TYPE>::operator=( rhs );
    return (*this);
  }

  /// equals operator that sets *this to a single value of any type
  template<class rTYPE> VectorT& operator=( const rTYPE& rhs )
  {
//    TYPE temp;
//    temp = rhs;
    this->std::vector<TYPE>::assign(this->size(), static_cast<TYPE>(rhs) );
    return (*this);
  }

  /// equals operator that sets *this to a single value of type TYPE
  VectorT& operator=( const TYPE& rhs )
  {
    this->std::vector<TYPE>::assign(this->size(), rhs );
    return (*this);
  }

  /// plus equals operator for other vectors of any type
  template<class rTYPE> VectorT& operator+=( const VectorT<rTYPE>& rhs );

  /// minus equals operator for other vectors of any type
  template<class rTYPE> VectorT& operator-=( const VectorT<rTYPE>& rhs );

  /// multiply equals operator for other vectors of any type
  template<class rTYPE> VectorT& operator*=( const VectorT<rTYPE>& rhs );

  /// divide equals operator for other vectors of any type
  template<class rTYPE> VectorT& operator/=( const VectorT<rTYPE>& rhs );

  /// plus equals operator for individual value of any type
  template<class rTYPE> VectorT& operator+=( const rTYPE& rhs );

  /// minus equals operator for individual value of any type
  template<class rTYPE> VectorT& operator-=( const rTYPE& rhs );

  /// multiply equals operator for individual value of any type
  template<class rTYPE> VectorT& operator*=( const rTYPE& rhs );

  /// divide equals operator for individual value of any type
  template<class rTYPE> VectorT& operator/=( const rTYPE& rhs );

  void SetValue( const TYPE* const rhs, const size_type n );

  /*
     void Trim()
     {
     std::vector<TYPE>( this->begin(), this->end()).swap(*this);

     }*/

  long int size() const { return std::vector<TYPE>::size(); }

//***** Range Checking ********************************************************
#if RANGE_CHECKING==1
  inline       TYPE& operator[](const int index)       { return std::vector<TYPE>::at(index);}
  inline const TYPE& operator[](const int index) const { return std::vector<TYPE>::at(index);}
#else
  inline       TYPE& operator[](const int index)       { return std::vector<TYPE>::operator[]( static_cast<unsigned int>(index) );}
  inline const TYPE& operator[](const int index) const { return std::vector<TYPE>::operator[]( static_cast<unsigned int>(index) );}
#endif


};

// *****************************************************************************
// **** Class Implementation ***************************************************
// *****************************************************************************

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
  if( this->size() == rhs.size() && !this->empty() )
  {
    TYPE* ptr = this->data();
    const rTYPE* rhs_ptr = &(rhs[0]);
    for( unsigned int a=0 ; a<(*this).size() ; ++a )
      *(ptr++) += *(rhs_ptr++);
  }
  else
//    std::__throw_length_error(__N("VectorT<TYPE>&
// VectorT<TYPE>::operator+=()"));
    throw std::exception();
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
  if( this->size() == rhs.size() && !this->empty() )
  {
    TYPE* ptr = this->data();
    const rTYPE* rhs_ptr = &(rhs[0]);
    for( unsigned int a=0 ; a<(*this).size() ; ++a )
      *(ptr++) -= *(rhs_ptr++);
  }
  else
//    std::__throw_length_error(__N("VectorT<TYPE>&
// VectorT<TYPE>::operator+=()"));
    throw std::exception();

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
  if( this->size() == rhs.size() && !this->empty() )
  {
    TYPE* ptr = this->data();
    const rTYPE* rhs_ptr = &(rhs[0]);
    for( unsigned int a=0 ; a<(*this).size() ; ++a )
      *(ptr++) *= *(rhs_ptr++);
  }
  else
    throw std::exception();

//    std::__throw_length_error(__N("VectorT<TYPE>&
// VectorT<TYPE>::operator+=()"));
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
  if( this->size() == rhs.size() && !this->empty() )
  {
    TYPE* ptr = this->data();
    const rTYPE* rhs_ptr = &(rhs[0]);
    for( unsigned int a=0 ; a<this->size() ; ++a )
      *(ptr++) /= *(rhs_ptr++);
  }
  else
    throw std::exception();

//    std::__throw_length_error(__N("VectorT<TYPE>&
// VectorT<TYPE>::operator+=()"));
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
  if( !this->empty() )
  {
    TYPE* ptr = this->data();
    for( unsigned int a=0 ; a<this->size() ; ++a )
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
  if( !this->empty() )
  {
    TYPE* ptr = this->data();
    for( unsigned int a=0 ; a<this->size() ; ++a )
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
  if( !this->empty() )
  {
    TYPE* ptr = this->data();
    for( unsigned int a=0 ; a<this->size() ; ++a )
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
  if( !this->empty() )
  {
    TYPE* ptr = this->data();
    for( unsigned int a=0 ; a<this->size() ; ++a )
      *(ptr++) /= rhs;
  }
  return (*this);
}

template<class TYPE>
inline void VectorT<TYPE>::SetValue( const TYPE* const rhs, const size_type n )
{
  this->resize(n);
  this->std::vector<TYPE>::assign(rhs, rhs+n);
}



#endif
