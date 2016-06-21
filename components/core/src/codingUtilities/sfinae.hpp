/*
 * sfinae.hpp
 *
 *  Created on: Dec 18, 2014
 *      Author: rrsettgast
 */

#ifndef SRC_CODINGUTILITIES_SFINAE_HPP_
#define SRC_CODINGUTILITIES_SFINAE_HPP_






#define HAS_MEMBER_FUNCTION(NAME,args...)\
template<typename TT >\
struct has_##NAME\
{\
private:\
template<typename U> static auto test(int) -> decltype( std::declval<U>().NAME(args), std::true_type() );\
template<typename U> static auto test(...) -> std::false_type ;\
public:\
static constexpr bool value = decltype(test<TT>(0))::value;\
};





namespace SFINAE
{
template<class T>
struct has_pointer_type
{
  template<class U> static char (&test(typename U::pointer const*))[1];
  template<class U> static char (&test(...))[2];
  static const bool value = (sizeof(test<T>(0)) == 1);
};



}

#endif /* SRC_CODINGUTILITIES_SFINAE_HPP_ */
