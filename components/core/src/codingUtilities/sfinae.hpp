/*
 * sfinae.hpp
 *
 *  Created on: Dec 18, 2014
 *      Author: rrsettgast
 */

#ifndef SRC_CODINGUTILITIES_SFINAE_HPP_
#define SRC_CODINGUTILITIES_SFINAE_HPP_

namespace geosx
{


#define HAS_MEMBER_DATA(NAME)\
template<typename TT>\
struct has_datamember_##NAME\
{\
private:\
  template<typename U> static constexpr auto test(int) -> decltype( std::declval<U>().NAME, bool() )\
  {\
    return std::is_member_object_pointer<decltype(&U::NAME)>::value;\
  }\
  template<typename U> static constexpr auto test(...) -> bool\
  {\
    return false;\
  }\
public:\
  static constexpr bool value = test<TT>(0);\
};

#define HAS_STATIC_MEMBER_DATA(NAME)\
template<typename TT>\
struct has_staticdatamember_##NAME\
{\
private:\
  template<typename U> static constexpr auto test(int) -> decltype( U::NAME, bool() )\
  {\
    return !std::is_function< decltype(U::NAME) >::value && !std::is_member_object_pointer<decltype(&U::NAME)>::value;\
  }\
  template<typename U> static constexpr auto test(...) -> bool\
  {\
    return false;\
  }\
public:\
  static constexpr bool value = test<TT>(0);\
};



#define HAS_MEMBER_FUNCTION(NAME, ...)\
template<typename TT >\
struct has_memberfunction_##NAME\
{\
private:\
  template<typename U> static constexpr auto test(int) -> decltype( std::declval<U>().NAME(__VA_ARGS__), bool() )\
  {\
    return std::is_member_function_pointer< decltype( &U::NAME ) >::value;\
  }\
  template<typename U> static constexpr auto test(...) -> bool \
  {\
    return false;\
  }\
public:\
  static constexpr bool value = test<TT>(0);\
};

#define HAS_STATIC_MEMBER_FUNCTION(NAME, ...)\
template<typename TT>\
struct has_staticmemberfunction_##NAME\
{\
private:\
  template<typename U> static constexpr auto test(int) -> decltype( U::NAME(__VA_ARGS__) , bool() )\
  {\
    return true;\
  }\
  template<typename U> static constexpr auto test(...) -> bool\
  {\
    return false;\
  }\
public:\
  static constexpr bool value = test<TT>(0);\
};


#define HAS_ENUM(NAME)\
template<typename TT >\
struct has_enum_##NAME\
{\
private:\
  template<typename U> static auto test(int) -> typename std::enable_if<std::is_enum< typename U::NAME >::value, std::true_type>::type;\
  template<typename U> static auto test(...) -> std::false_type;\
public:\
  static constexpr bool value = decltype(test<TT>(0))::value;\
};

#define HAS_ALIAS(NAME)\
template<typename TT>\
struct has_alias_##NAME\
{\
private:\
  template<typename U> static auto test(int) -> typename std::enable_if<!std::is_enum< typename U::NAME >::value, std::true_type>::type;\
  template<typename U> static auto test(...) -> std::false_type ;\
public:\
  static constexpr bool value = !(std::is_same<decltype(test<TT>(0)),std::false_type>::value);\
};



#define VA_LIST(...) __VA_ARGS__

#define CONDITIONAL_VIRTUAL_FUNCTION(NAME,RTYPE,PARAMS,ARGS)\
template<typename U=T, bool has = has_memberfunction_##NAME<U>::value >\
struct wrapper##NAME\
{\
  RTYPE f(...) {return RTYPE();}\
};\
template<typename U>\
struct wrapper##NAME<U,true>\
{\
  RTYPE f( DataObject<U> * const obj ,PARAMS )\
  {\
    return (*obj).m_data.NAME(ARGS);\
  }\
};\
virtual RTYPE NAME(PARAMS) override final\
{\
  wrapper##NAME<T> temp;\
  return temp.f(this ,ARGS);\
}

#define CONDITIONAL_VIRTUAL_FUNCTION0(NAME,RTYPE)\
template<typename U=T, bool has = has_memberfunction_##NAME<U>::value >\
struct wrapper##NAME\
{\
  RTYPE f(...) {return 0;}\
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


template<class TT>
struct has_pointer_type
{
  template<class U> static char (&test(typename U::pointer const*))[1];
  template<class U> static char (&test(...))[2];
  static const bool value = (sizeof(test<TT>(0)) == 1);
};



}

#endif /* SRC_CODINGUTILITIES_SFINAE_HPP_ */
