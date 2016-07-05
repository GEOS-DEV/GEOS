/*
 * DataObject.h
 *
 *  Created on: Jun 8, 2016
 *      Author: settgast
 */

#ifndef CORE_SRC_DATAREPOSITORY_DATAOBJECT_HPP_
#define CORE_SRC_DATAREPOSITORY_DATAOBJECT_HPP_

#include "sidre/sidre.hpp"
#include "DataTypes.hpp"
#include "DataObjectBase.hpp"
#include "codingUtilities/sfinae.hpp"
namespace geosx
{

template< typename T >
class DataObject : public DataObjectBase
{

public:
  DataObject( std::string const & name ):
    DataObjectBase(name)
  {}

  virtual ~DataObject() noexcept override final{}

  DataObject( DataObject const & source ):
    m_data(source.m_data)
  {}

  DataObject( DataObject&& source ):
    m_data( std::move(source.m_data) )
  {}

  DataObject& operator=( DataObject const & source )
  {
    m_data = source.m_data;
  }

  DataObject& operator=( DataObject && source )
  {
    m_data = std::move(source.m_data);
  }

  virtual const std::type_info& get_typeid() const noexcept override final
  {
    return typeid(T);
  }



  HAS_MEMBER_FUNCTION(size,)
  CONDITIONAL_VIRTUAL_FUNCTION0(size,std::size_t)


  HAS_MEMBER_FUNCTION(resize, std::size_t(1) )
  CONDITIONAL_VIRTUAL_FUNCTION( resize , void, VA_LIST(std::size_t a), VA_LIST(a) )




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

  DataObject() = delete;
};



} /* namespace geosx */

#endif /* CORE_SRC_DATAREPOSITORY_DATAOBJECT_HPP_ */
