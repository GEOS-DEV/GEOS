/*
 * ArrayWrapper.hpp
 *
 *  Created on: Jul 30, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_ARRAY_MULTIDIMENSIONALARRAY_HPP_
#define SRC_COMPONENTS_CORE_SRC_ARRAY_MULTIDIMENSIONALARRAY_HPP_
#include<vector>
#include<iostream>
#if 0
#include "common/DataTypes.hpp"
#else
using  int32     = std::int32_t;
using uint32     = std::uint32_t;
using  int64     = std::int64_t;
using uint64     = std::uint64_t;
using std_size_t = std::size_t;
using string     = std::string;
#endif
namespace multidimensionalArray
{

template< typename T, int NDIM >
class ArrayAccessor
{
public:

  using rtype = ArrayAccessor<T,NDIM-1> &;

  ArrayAccessor() = delete;

  ArrayAccessor( ArrayAccessor const & source ):
    m_data(source.m_data),
    m_lengths(source.m_lengths),
    m_stride(source.m_stride),
    m_childInterface(source.m_childInterface)
  {}

  ArrayAccessor( ArrayAccessor && source ):
  m_data(source.m_data),
  m_lengths(source.m_lengths),
  m_stride(source.m_stride),
  m_childInterface(source.m_childInterface)
  {}
  
  ArrayAccessor & operator=( ArrayAccessor const & ) = delete;
  ArrayAccessor & operator=( ArrayAccessor && ) = delete;

  ArrayAccessor( T * const data, int64 const * const length ):
    m_data(data),
    m_lengths(length),
    m_stride(1),
    m_childInterface(ArrayAccessor<T,NDIM-1>( m_data, &(length[1]) ) )
  {
    for( int a=1 ; a<NDIM ; ++a )
    {
      m_stride *= m_lengths[a];
    }
  }

  inline rtype operator[](int64 const index)
  {
    m_childInterface.m_data = &(m_data[index*m_stride]);
    return m_childInterface;
  }

  T * m_data;
  int64 const * m_lengths;
  int64 m_stride = 1;
  ArrayAccessor<T,NDIM-1> m_childInterface;
};

template< typename T >
class ArrayAccessor<T,1>
{
public:
  using rtype = T &;

  ArrayAccessor() = delete;

  ArrayAccessor( ArrayAccessor const & source ):
  m_data(source.m_data),
  m_lengths(source.m_lengths)
  {}

  ArrayAccessor( ArrayAccessor && source ):
  m_data(source.m_data),
  m_lengths(source.m_lengths)
  {}

  ArrayAccessor & operator=( ArrayAccessor const & ) = delete;
  ArrayAccessor & operator=( ArrayAccessor && ) = delete;

  ArrayAccessor( T * const data, int64 const * const length ):
    m_data(data),
    m_lengths(length)
  {}


  inline T& operator[](int64 const index)
  {
    return m_data[index];
  }

  T * m_data;
  int64 const * m_lengths;
};














template< typename T, int NDIM, typename memBlock = std::vector<T> >
class Array
{
public:

  using rtype = ArrayAccessor<T,NDIM-1> &;

  Array() = delete;


  template< class U=T>
  Array( int64 const lengths[NDIM] ):
  m_memory(),
  m_lengths(),
  m_interface( ArrayAccessor<T,NDIM>(nullptr,lengths) )
  {
    int64 size = 1;
    for( int a=0 ; a<NDIM ; ++a )
    {
      m_lengths[a] = lengths[a];
      size *= lengths[a];
    }
    m_memory.resize(size);

    m_interface.m_data    = m_memory.data();
    m_interface.m_lengths = m_lengths;
  }

  ~Array() = default;

  Array( Array const & ) = delete;
  Array( Array && source ) = delete;

  Array& operator=( Array const & rhs ) = delete;
  Array& operator=( Array && rhs ) = delete;

  //***** Accessors **********************************************************


  inline rtype operator[](int64 const index)
  {
    return m_interface[index];
  }

private:
  memBlock m_memory;
  int64 m_lengths[NDIM] = {0};
  ArrayAccessor<T,NDIM> m_interface;

};

} /* namespace arraywrapper */

#endif /* SRC_COMPONENTS_CORE_SRC_ARRAY_MULTIDIMENSIONALARRAY_HPP_ */
