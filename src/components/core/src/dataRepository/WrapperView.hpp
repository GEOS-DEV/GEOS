/*
 * DataObject.h
 *
 *  Created on: Jun 8, 2016
 *      Author: settgast
 */

#ifndef GEOSX_DATAREPOSITORY_WRAPPERVIEW_HPP_
#define GEOSX_DATAREPOSITORY_WRAPPERVIEW_HPP_

#include "sidre/sidre.hpp"
#include "KeyNames.hpp"
#include "common/DataTypes.hpp"
#include "SFINAE_Macros.hpp"
#include <type_traits>

#include "WrapperViewBase.hpp"

namespace geosx
{
namespace dataRepository
{

template< typename T >
class WrapperView : public WrapperViewBase
{

public:
  explicit WrapperView( std::string const & name,
                    SynchronizedGroup * const parent ) :
    WrapperViewBase(name,parent)
  {
    // set up properties of sidre::DataView
    if( std::is_array<T>::value )
    {}
    else
    {
      getSidreView()->setExternalDataPtr( nullptr );
    }
  }

  virtual ~WrapperView() noexcept override final {}

  WrapperView( WrapperView const & source ) :
    m_data(source.m_data)
  {}

  WrapperView( WrapperView&& source ) :
    m_data( std::move(source.m_data) )
  {}

  WrapperView& operator=( WrapperView const & source )
  {
    m_data = source.m_data;
  }

  WrapperView& operator=( WrapperView && source )
  {
    m_data = std::move(source.m_data);
  }


  static std::unique_ptr<WrapperViewBase> Factory( std::string const & name,
                                               SynchronizedGroup * const parent )
  {
    return std::move(std::make_unique<WrapperView<T> >( name, parent ) );
  }


  virtual const std::type_info& get_typeid() const noexcept override final
  {
    return typeid(T);
  }

  struct empty_wrapper
  {
    HAS_MEMBER_FUNCTION(empty,bool,const,,)
    template<class U = T>
    static typename std::enable_if<has_memberfunction_empty<U>::value, bool>::type empty(WrapperView const * parent)
    {
      return parent->m_data.empty();
    }
    template<class U = T>
    static typename std::enable_if<!has_memberfunction_empty<U>::value, bool>::type empty(WrapperView const * parent)
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
    static typename std::enable_if<has_memberfunction_size<U>::value, std::size_t>::type size(WrapperView const * parent)
    {
      return parent->m_data.size();
    }
    template<class U = T>
    static typename std::enable_if<!has_memberfunction_size<U>::value, std::size_t>::type size(WrapperView const * )
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
    static typename std::enable_if<has_memberfunction_reserve<U>::value, void>::type reserve(WrapperView * const parent, std::size_t new_cap)
    {
      return parent->m_data.reserve(new_cap);
    }
    template<class U = T>
    static typename std::enable_if<!has_memberfunction_reserve<U>::value, void>::type reserve(WrapperView * const, std::size_t )
    {
      return; //parent->m_data;
    }
  };
  virtual void reserve( std::size_t new_cap ) override final
  {
    return reserve_wrapper::reserve(this, new_cap);
  }
//  CONDITIONAL_VIRTUAL_FUNCTION( Wrapper<T>,reserve , void,, VA_LIST(std::size_t a), VA_LIST(a) )


  HAS_MEMBER_FUNCTION(capacity,std::size_t,const,,)
  CONDITIONAL_VIRTUAL_FUNCTION0(WrapperView<T>,capacity,std::size_t,const)

  HAS_MEMBER_FUNCTION(max_size,std::size_t,const,,)
  CONDITIONAL_VIRTUAL_FUNCTION0(WrapperView<T>,max_size,std::size_t,const)

  HAS_MEMBER_FUNCTION(clear,void,,,)
  CONDITIONAL_VIRTUAL_FUNCTION0(WrapperView<T>,clear,void,)

  HAS_MEMBER_FUNCTION(insert,void,,,)
  CONDITIONAL_VIRTUAL_FUNCTION0(WrapperView<T>,insert,void,)

//  HAS_MEMBER_FUNCTION(resize, void, std::size_t(1) )
  HAS_MEMBER_FUNCTION(resize, void, , VA_LIST(std::size_t), VA_LIST(std::size_t(1)) )
  CONDITIONAL_VIRTUAL_FUNCTION( WrapperView<T>,resize, void,, VA_LIST(std::size_t a), VA_LIST(a) )

  template<class U=T, bool has = cxx_utilities::has_pointer_type<U>::value >
  struct Get_Type
  {
    typedef U*       type;
    typedef const U* const_type;
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
  typename std::enable_if<has_memberfunction_data<U>::value && !std::is_same<U,string>::value, rtype>::type data()
  {
    return m_data.data();
  }
  template<class U = T>
  typename std::enable_if<has_memberfunction_data<U>::value && std::is_same<U,string>::value, U&>::type data()
  {
    return m_data;
  }
  template<class U = T>
  typename std::enable_if<!has_memberfunction_data<U>::value, U*>::type data()
  {
    return &m_data;
  }

  template<class U = T>
  typename std::enable_if<has_memberfunction_data<U>::value && !std::is_same<U,string>::value, rtype_const>::type data() const
  {
    return m_data.data();
  }
  template<class U = T>
  typename std::enable_if<has_memberfunction_data<U>::value && std::is_same<U,string>::value, U const &>::type data() const
  {
    return m_data;
  }
  template<class U = T>
  typename std::enable_if<!has_memberfunction_data<U>::value, U const *>::type data() const
  {
    return &m_data;
  }

  T& reference()
  { return m_data; }

  T const & reference() const
  { return m_data; }

private:
public:
  T m_data;

  WrapperView() = delete;
};


}
} /* namespace geosx */

#endif /* CORE_SRC_DATAREPOSITORY_DATAOBJECT_HPP_ */
