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
 * @file ReferenceWrapper.hpp
 * This file contains the class definition of ReferenceWrapper.
 */

#ifndef SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_REFERENCEWRAPPER_HPP_
#define SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_REFERENCEWRAPPER_HPP_

#include "SFINAE_Macros.hpp"
#include <type_traits>

namespace geosx
{

/**
 * @class ReferenceWrapper
 * @tparam Type that is wrapped
 *
 * This class manages a pointer to the templated type, but provides a reference-like interface,
 * thus negating the requirement to dereference the object. The primary use for this is for nested
 * arrays that hold pointers at the last level, but allows for reference-like usage. For instance,
 * consider a collection of object that you would like to refer to thought an array of pointers.
 *
 * array< ReferenceWrapper< array< double > > > arr;
 *
 * where the array::operator[] exists. The ReferenceWrapper allows
 *
 * arr[index1][index2]
 *
 * note: this is really only useful for heavy array that hold their own data. For light arrays that
 * hold pointers to their data, then this is unnecessary as a copy of the array does not trigger a
 * deep copy.
 */
template< typename T >
class ReferenceWrapper
{
public:

  /**
   * @brief default constructor sets m_ref to nullptr
   */
  ReferenceWrapper():
    m_ref(nullptr)
  {}


  /**
   * @brief constructor sets m_ref to address of input
   * @param source object to wrap
   */
  ReferenceWrapper( T & source ) noexcept:
    m_ref(&source)
  {}


  /**
   * @brief default destructor
   */
  ~ReferenceWrapper() = default;

  /**
   * @brief copy constructor copies the source m_ref to the new m_ref
   * @param source
   */
  ReferenceWrapper( ReferenceWrapper const & source ):
    m_ref( source.m_ref )
  {}


  /**
   * @brief move constructor copies the source m_ref to the new m_ref
   * @param source
   */
  ReferenceWrapper( ReferenceWrapper && source ):
    m_ref( source.m_ref )
  {
    source.m_ref = nullptr;
  }

  /**
   * @brief assigment operator sets the value that m_ref refers to to the value of the rhs
   * @param rhs the rhs value to be copied
   * @return
   */
  inline ReferenceWrapper & operator=( T const & rhs )
  {
    *m_ref = rhs;
    return *this;
  }

  /**
   * @tparam T_RHS type of the rhs
   * @brief assignment operator calls m_ref->operator=() to allow for any type on the rhs
   *        if m_ref->operator=() has a valid overload for T_RHS
   * @tparam U dummy template parameter to enable SFINAE stuff
   * @param rhs value to be copied
   * @return *this
   */
  template< typename T_RHS, typename U=T >
  inline
  typename std::enable_if< !std::is_const<U>::value,ReferenceWrapper &>::type
  operator=( T_RHS const & rhs )
  {
    *m_ref = rhs;
    return *this;
  }

  /**
   * @brief move assignment operator sets the value that m_ref refers to to the value of the rhs
   * @param rhs the rhs value to be moved
   * @return
   */
  inline ReferenceWrapper & operator=( T && source )
  {
    *m_ref = std::move(source);
    return *this;
  }

  /**
   * @brief user defined conversion to T &
   */
  inline operator T & ()
  {
    return *m_ref;
  }

  /**
   * @brief user defined conversion to T const &
   */
  inline operator T const & () const
  {
    return *m_ref;
  }

  /**
   * @brief accessor function to set the address that m_ref points to
   * @param source reference to object that wrapper will refer to
   */
  inline void set( T & source)
  {
    m_ref = &source;
  }

  /**
   * @brief accessor function to set the address that m_ref points to
   * @param source pointer to object that wrapper will refer to
   */
  inline void set( T * source)
  {
    m_ref = source;
  }

  /**
   * @brief accessor for m_ref
   * @return reference to wrapped value
   */
  inline T & get()
  {
    return *m_ref;
  }

  /**
   * @brief const accessor for m_ref
   * @return const reference to wrapped value
   */
  inline T const & get() const
  {
    return *m_ref;
  }

  /**
   * @brief a pass thru square bracket operator which calls underlying T::operator[]
   * @tparam U dummy type to allow for SFINAE evaluation of availability of T::operator[]
   * @param i index to pass into the T::operator[]
   * @return the return type of T::operator[]
   */
  template< typename INDEX_TYPE, typename U = T >
  inline decltype( std::declval< U >()[1] )
  operator[]( INDEX_TYPE const i )
  {
    return (*m_ref)[i];
  }

  /**
   * @brief a const pass thru square bracket operator which calls underlying T::operator[]
   * @tparam U dummy type to allow for SFINAE evaluation of availability of T::operator[]
   * @param i index to pass into the T::operator[] const
   * @return the return type of T::operator[] const
   */
  template< typename INDEX_TYPE, typename U = T >
  inline decltype( std::declval< U const >()[1] )
  operator[]( INDEX_TYPE const i ) const
  {
    return (*m_ref)[i];
  }


  /**
   * @brief a pass thru parenthesis  operator which calls underlying T::operator()
   * @tparam variadic types to pass through to T::operator()
   * @param args variadic params to pass through to T::operator()
   * @return the return type of T::operator()
   */
  template<typename... ARGS>
  inline typename std::result_of<T&(ARGS&&...)>::type
  operator()(ARGS&&... args)
  {
    return m_ref->operator()( std::forward<ARGS>(args)...) ;
  }

  /**
   * @brief a pass thru parenthesis  operator which calls underlying T::operator()
   * @tparam variadic types to pass through to T::operator()
   * @param args variadic params to pass through to T::operator()
   * @return the return type of T::operator() const
   */
  template<typename... ARGS>
  inline typename std::result_of<T const&(ARGS&&...)>::type
  operator()(ARGS&&... args) const
  {
    return m_ref->operator()( std::forward<ARGS>(args)...) ;
  }


private:
  /// pointer to the address of the object that we would like to wrap
  T * m_ref;
};




} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_REFERENCEWRAPPER_HPP_ */
