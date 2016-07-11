/*
 * DataObject.h
 *
 *  Created on: Jun 8, 2016
 *      Author: settgast
 */

#ifndef CORE_SRC_DATAREPOSITORY_DATAOBJECT_HPP_
#define CORE_SRC_DATAREPOSITORY_DATAOBJECT_HPP_

#include "sidre/sidre.hpp"
#include "../DataTypes.hpp"
#include "codingUtilities/sfinae.hpp"
#include "WrapperBase.hpp"
namespace geosx
{
namespace dataRepository
{

template< typename T >
class Wrapper : public WrapperBase
{

public:
  Wrapper( std::string const & name ):
    WrapperBase(name)
  {}

  virtual ~Wrapper() noexcept override final{}

  Wrapper( Wrapper const & source ):
    m_data(source.m_data)
  {}

  Wrapper( Wrapper&& source ):
    m_data( std::move(source.m_data) )
  {}

  Wrapper& operator=( Wrapper const & source )
  {
    m_data = source.m_data;
  }

  Wrapper& operator=( Wrapper && source )
  {
    m_data = std::move(source.m_data);
  }

  virtual const std::type_info& get_typeid() const noexcept override final
  {
    return typeid(T);
  }



  HAS_MEMBER_FUNCTION(empty,)
  CONDITIONAL_VIRTUAL_FUNCTION0(Wrapper<T>,empty,bool,const)

  HAS_MEMBER_FUNCTION(size,)
  CONDITIONAL_VIRTUAL_FUNCTION0(Wrapper<T>,size,std::size_t,const)

  HAS_MEMBER_FUNCTION(reserve, std::size_t(1) )
  CONDITIONAL_VIRTUAL_FUNCTION( Wrapper<T>,reserve , void,, VA_LIST(std::size_t a), VA_LIST(a) )

  HAS_MEMBER_FUNCTION(capacity,)
  CONDITIONAL_VIRTUAL_FUNCTION0(Wrapper<T>,capacity,std::size_t,const)

  HAS_MEMBER_FUNCTION(max_size,)
  CONDITIONAL_VIRTUAL_FUNCTION0(Wrapper<T>,max_size,std::size_t,const)

  HAS_MEMBER_FUNCTION(clear,)
  CONDITIONAL_VIRTUAL_FUNCTION0(Wrapper<T>,clear,void,)

  HAS_MEMBER_FUNCTION(insert,)
  CONDITIONAL_VIRTUAL_FUNCTION0(Wrapper<T>,insert,void,)

  HAS_MEMBER_FUNCTION(resize, std::size_t(1) )
  CONDITIONAL_VIRTUAL_FUNCTION( Wrapper<T>,resize , void,, VA_LIST(std::size_t a), VA_LIST(a) )


#if 0
  using rtype = T &;
  using rtype_const = T const &;

  rtype_const getObjectData() const
  { return m_data; }

  rtype getObjectData()
  { return m_data; }


#else


  template<class U=T, bool has = has_pointer_type<U>::value >
  struct Get_Type
  {
    typedef U&       type;
    typedef const U& const_type;
  };


  template<class U>
  struct Get_Type<U, true>
  {
    typedef typename U::pointer       type;
    typedef typename U::const_pointer const_type;
  };


  using rtype       = typename Get_Type<T>::type;
  using rtype_const = typename Get_Type<T>::const_type;


  HAS_MEMBER_FUNCTION(data,)
  template<class U = T>
  typename std::enable_if<has_memberfunction_data<U>::value, rtype>::type getObjectData()
  {
    return m_data.data();
  }
  template<class U = T>
  typename std::enable_if<!has_memberfunction_data<U>::value, rtype>::type getObjectData()
  {
      return m_data;
  }


#endif

private:
public:
  T m_data;

  Wrapper() = delete;
};


}
} /* namespace geosx */

#endif /* CORE_SRC_DATAREPOSITORY_DATAOBJECT_HPP_ */
