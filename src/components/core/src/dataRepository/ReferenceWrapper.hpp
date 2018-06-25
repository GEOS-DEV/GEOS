/*
 * ReferenceWrapper.hpp
 *
 *  Created on: Jun 24, 2018
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_REFERENCEWRAPPER_HPP_
#define SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_REFERENCEWRAPPER_HPP_

#include "SFINAE_Macros.hpp"
#include <type_traits>

namespace geosx
{

template< typename T >
class ReferenceWrapper
{
public:
  ReferenceWrapper() = delete;

  ReferenceWrapper( T & source ) noexcept:
    m_ref(source)
  {}

  ReferenceWrapper( T && source ) = delete;

  ~ReferenceWrapper() = default;

  operator T & ()
  {
    return m_ref;
  }

  operator T const & () const
  {
    return m_ref;
  }

  T       & get()       { return m_ref; }
  T const & get() const { return m_ref; }

  ReferenceWrapper & operator=( T const & rhs )
  {
    m_ref = rhs;
    return *this;
  }

  template< typename U = T >
  decltype( std::declval< U >()[1] )
  operator[]( int i )
  {
    return m_ref[i];
  }

  template< typename U = T >
  decltype( std::declval< U const >()[1] )
  operator[]( int i ) const
  {
    return m_ref[i];
  }

  template<typename... ARGS>
  typename std::result_of<T&(ARGS&&...)>::type
  operator()(ARGS&&... args) const
  {
    return std::__invoke(get(), std::forward<ARGS>(args)...);
  }


private:
  T & m_ref;
};




} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_REFERENCEWRAPPER_HPP_ */
