/*
 * sfinae.hpp
 *
 *  Created on: Dec 18, 2014
 *      Author: rrsettgast
 */

#ifndef SRC_CODINGUTILITIES_SFINAE_HPP_
#define SRC_CODINGUTILITIES_SFINAE_HPP_


namespace SFINAE
{
template<class T>
struct has_pointer_type
{
  template<class U> static char (&test(typename U::pointer const*))[1];
  template<class U> static char (&test(...))[2];
  static const bool value = (sizeof(test<T>(0)) == 1);
};


template<typename T>
struct has_resize_method
{
private:

  template<typename U> static auto test(int) -> decltype( std::declval<U>().resize(1), std::true_type() );
  template<typename U> static auto test(...) -> std::false_type ;
public:
  static constexpr bool value = !(std::is_same<decltype(test<T>(0)),std::false_type>::value);
//  static constexpr bool value = decltype(test<T>(0))::value;
};

#define RESIZE()\
template<typename U=T, bool has = SFINAE::has_resize_method<U>::value >\
struct resized\
{\
  void f( const std::size_t, void *  ) {}\
};\
template<typename U>\
struct resized<U,true>\
{\
  void f( const std::size_t newsize, DataObject<U> * const obj )\
  {\
    (*obj).m_data.resize(newsize);\
  }\
};\
virtual void resize(const std::size_t newsize ) override final\
{\
  resized<T> temp;\
  temp.f(newsize,this);\
}


}

#endif /* SRC_CODINGUTILITIES_SFINAE_HPP_ */
