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



#define DEFINE_WRAPPING_FUNCTION0(NAME,RTYPE)\
template<typename U=T, bool has = has_memberfunction_##NAME<U>::value >\
struct wrapper##NAME\
{\
  RTYPE f( void * ) {return 0;}\
};\
template<typename U>\
struct wrapper##NAME<U,true>\
{\
  RTYPE f( DataObject<U> * const obj )\
  {\
    return (*obj).m_data.NAME();\
  }\
};\
virtual RTYPE NAME() override final\
{\
  wrapper##NAME<T> temp;\
  return temp.f(this);\
}

#define DEFINE_WRAPPING_FUNCTION1(NAME,ARG1)\
template<typename U=T, bool has = has_memberfunction_##NAME<U>::value >\
struct wrapper##NAME\
{\
  void f( void *, ARG1 ) {}\
};\
template<typename U>\
struct wrapper##NAME<U,true>\
{\
  void f( DataObject<U> * const obj, ARG1 arg1 )\
  {\
    (*obj).m_data.NAME( arg1 );\
  }\
};\
virtual void NAME( ARG1 arg1 ) override final\
{\
  wrapper##NAME<T> temp;\
  temp.f(this,arg1);\
}



template< typename T >
class DataObject : public DataObjectBase
{

public:
  DataObject( const std::string& name ):
    DataObjectBase(name)
  {}

  virtual ~DataObject() noexcept override final{}


  virtual const std::type_info& get_typeid() const noexcept override final
  {
    return typeid(T);
  }



  HAS_MEMBER_FUNCTION(size,)
  HAS_MEMBER_FUNCTION(resize, std::size_t(1) )
#if 1
  DEFINE_WRAPPING_FUNCTION0(size,std::size_t)
  DEFINE_WRAPPING_FUNCTION1( resize , std::size_t )
#else
  template<typename U=T, bool has = has_memberfunction_size<U>::value >
  struct size_wrapper
  {
    int f( void * ) {return 0;}
  };
  template<typename U>
  struct size_wrapper<U,true>
  {
    int f( DataObject<U> * const obj )
    {
      return (*obj).m_data.size();
    }
  };
  virtual std::size_t size() override final
  {
    size_wrapper<T> temp;
    return temp.f(this);
  }

  template<typename U=T, bool has = has_memberfunction_resize<U>::value >
  struct resize_wrapper
  {
    void f( void *, std::size_t ) {}
  };
  template<typename U>
  struct resize_wrapper<U,true>
  {
    void f( DataObject<U> * const obj, std::size_t arg1 )
    {
      (*obj).m_data.resize( arg1 );
    }
  };
  virtual void resize( std::size_t arg1 ) override final
  {
    resize_wrapper<T> temp;
    temp.f(this,arg1);
  }
#endif




#ifndef OBJECTDATA_PTR_RETURN
  using rtype = T &;
  using const_rtype = T const &;

  const_rtype getObjectData() const
  { return m_data; }

  rtype getObjectData()
  { return m_data; }


#else

  /*
  template<class TYPE>
  struct has_pointer_type
  {
    template<class U> static char (&test(typename U::pointer const*))[1];
    template<class U> static char (&test(...))[2];
    static const bool value = (sizeof(test<TYPE>(0)) == 1);
  };
*/

  template<class U=T, bool has = sfinae::has_pointer_type<U>::value >
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
  using const_rtype = typename Get_Type<T>::const_type;

  rtype getObjectData()
  { return m_data.data(); }

  const_rtype getObjectData() const
  { return m_data.data(); }
#endif

private:
public:
  T m_data;

  DataObject();
  DataObject( const DataObject&);

};



} /* namespace geosx */

#endif /* CORE_SRC_DATAREPOSITORY_DATAOBJECT_HPP_ */
