/*
 * DataObject.h
 *
 *  Created on: Jun 8, 2016
 *      Author: settgast
 */

#ifndef CORE_SRC_DATAREPOSITORY_DATAOBJECT_HPP_
#define CORE_SRC_DATAREPOSITORY_DATAOBJECT_HPP_

#include "sidre/sidre.hpp"

#include "../../common/DataTypes.hpp"
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
  explicit Wrapper( std::string const & name ):
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

  struct empty_wrapper
  {
    HAS_MEMBER_FUNCTION(empty,bool,const,,)
    template<class U = T>
    static typename std::enable_if<has_memberfunction_empty<U>::value, bool>::type empty(Wrapper const * parent)
    {
      return parent->m_data.empty();
    }
    template<class U = T>
    static typename std::enable_if<!has_memberfunction_empty<U>::value, bool>::type empty(Wrapper const * parent)
    {
      return parent;
    }
  };
  virtual bool empty() const override final
  {
    return empty_wrapper::empty(this);
  }

  struct size_wrapper
  {
    HAS_MEMBER_FUNCTION(size,std::size_t,const,,)
    template<class U = T>
    static typename std::enable_if<has_memberfunction_size<U>::value, std::size_t>::type size(Wrapper const * parent)
    {
      return parent->m_data.size();
    }
    template<class U = T>
    static typename std::enable_if<!has_memberfunction_size<U>::value, std::size_t>::type size(Wrapper const * )
    {
      return 0;//parent->m_data;
    }
  };
  virtual std::size_t size() const override final
  {
    return size_wrapper::size(this);
  }


  struct reserve_wrapper
  {
    HAS_MEMBER_FUNCTION(reserve, void, ,VA_LIST(std::size_t),VA_LIST(std::size_t(1)) )
    template<class U = T>
    static typename std::enable_if<has_memberfunction_reserve<U>::value, void>::type reserve(Wrapper * const parent, std::size_t new_cap)
    {
      return parent->m_data.reserve(new_cap);
    }
    template<class U = T>
    static typename std::enable_if<!has_memberfunction_reserve<U>::value, void>::type reserve(Wrapper * const, std::size_t )
    {
      return ;//parent->m_data;
    }
  };
  virtual void reserve( std::size_t new_cap ) override final
  {
    return reserve_wrapper::reserve(this, new_cap);
  }
//  CONDITIONAL_VIRTUAL_FUNCTION( Wrapper<T>,reserve , void,, VA_LIST(std::size_t a), VA_LIST(a) )


  HAS_MEMBER_FUNCTION(capacity,std::size_t,const,,)
  CONDITIONAL_VIRTUAL_FUNCTION0(Wrapper<T>,capacity,std::size_t,const)

  HAS_MEMBER_FUNCTION(max_size,std::size_t,const,,)
  CONDITIONAL_VIRTUAL_FUNCTION0(Wrapper<T>,max_size,std::size_t,const)

  HAS_MEMBER_FUNCTION(clear,void,,,)
  CONDITIONAL_VIRTUAL_FUNCTION0(Wrapper<T>,clear,void,)

  HAS_MEMBER_FUNCTION(insert,void,,,)
  CONDITIONAL_VIRTUAL_FUNCTION0(Wrapper<T>,insert,void,)

//  HAS_MEMBER_FUNCTION(resize, void, std::size_t(1) )
  HAS_MEMBER_FUNCTION(resize, void, , VA_LIST(std::size_t), VA_LIST(std::size_t(1)) )
  CONDITIONAL_VIRTUAL_FUNCTION( Wrapper<T>,resize , void,, VA_LIST(std::size_t a), VA_LIST(a) )

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


  HAS_MEMBER_FUNCTION(data,rtype,,,)
  template<class U = T>
  typename std::enable_if<has_memberfunction_data<U>::value, rtype>::type data()
  {
    return m_data.data();
  }
  template<class U = T>
  typename std::enable_if<!has_memberfunction_data<U>::value, rtype>::type data()
  {
      return m_data;
  }

  T& dataRef()
  { return m_data; }

  T const & dataRef() const
  { return m_data; }

private:
public:
  T m_data;

  Wrapper() = delete;
};


}
} /* namespace geosx */

#endif /* CORE_SRC_DATAREPOSITORY_DATAOBJECT_HPP_ */
