/*
 * ArrayWrapper.hpp
 *
 *  Created on: Jul 30, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_ARRAY_ARRAYWRAPPER_HPP_
#define SRC_COMPONENTS_CORE_SRC_ARRAY_ARRAYWRAPPER_HPP_

#include "common/DataTypes.hpp"

namespace arraywrapper
{

template< typename T, int NDIM >
class ArrayInterface
{
public:

  using rtype = ArrayInterface<T,NDIM-1> &;

  ArrayInterface() = delete;

  ArrayInterface( ArrayInterface const & source ):
    m_data(source.m_data),
    m_lengths(source.m_lengths),
    m_stride(source.m_stride),
    m_childInterface(source.m_childInterface)
  {}

  ArrayInterface( ArrayInterface && source ):
  m_data(source.m_data),
  m_lengths(source.m_lengths),
  m_stride(source.m_stride),
  m_childInterface(source.m_childInterface)
  {}
  
  ArrayInterface & operator=( ArrayInterface const & ) = delete;
  ArrayInterface & operator=( ArrayInterface && ) = delete;

  ArrayInterface( T * const data, int64 const * const length ):
    m_data(data),
    m_lengths(length),
    m_stride(1),
    m_childInterface(ArrayInterface<T,NDIM-1>( m_data, &(length[1]) ) )
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
  ArrayInterface<T,NDIM-1> m_childInterface;
};

template< typename T >
class ArrayInterface<T,1>
{
public:
  using rtype = T &;

  ArrayInterface() = delete;

  ArrayInterface( ArrayInterface const & source ):
  m_data(source.m_data),
  m_lengths(source.m_lengths)
  {}

  ArrayInterface( ArrayInterface && source ):
  m_data(source.m_data),
  m_lengths(source.m_lengths)
  {}

  ArrayInterface & operator=( ArrayInterface const & ) = delete;
  ArrayInterface & operator=( ArrayInterface && ) = delete;

  ArrayInterface( T * const data, int64 const * const length ):
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
class marray
{
public:

  using rtype = ArrayInterface<T,NDIM-1> &;

  marray() = delete;


  template< class U=T>
  marray( int64 const lengths[NDIM] ):
  m_memory(),
  m_lengths(),
  m_interface( ArrayInterface<T,NDIM>(nullptr,lengths) )
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

  ~marray() {}

  marray( marray const & ) = delete;
  marray( marray && source );

  marray& operator=( marray const & rhs );
  marray& operator=( marray && rhs );

  //***** Accessors **********************************************************


  inline rtype operator[](int64 const index)
  {
    return m_interface[index];
  }

private:
  memBlock m_memory;
  int64 m_lengths[NDIM] = {0};
  ArrayInterface<T,NDIM> m_interface;

};

} /* namespace arraywrapper */

#endif /* SRC_COMPONENTS_CORE_SRC_ARRAY_ARRAYWRAPPER_HPP_ */
