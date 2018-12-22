/*
 * WrapperDefaultValueHelper.hpp
 *
 *  Created on: Dec 17, 2018
 *      Author: settgast
 */

#ifndef CORECOMPONENTS_DATAREPOSITORY_DEFAULTVALUE_HPP_
#define CORECOMPONENTS_DATAREPOSITORY_DEFAULTVALUE_HPP_

#include "SFINAE_Macros.hpp"

namespace geosx
{
namespace dataRepository
{

namespace wrapperDefaultValue
{


template< typename T >
struct is_arithmetic
{
  static constexpr bool value = std::is_same<T,int>::value ||
      std::is_same<T,long int>::value ||
      std::is_same<T,long long int>::value ||
      std::is_same<T,unsigned int>::value ||
      std::is_same<T,unsigned long int>::value ||
      std::is_same<T,unsigned long long int>::value ||
      std::is_floating_point<T>::value;
};

template< typename T, typename ENABLE=void >
struct Helper
{
  static constexpr bool has_default_value = false;
};

template< typename T >
struct Helper< T, typename std::enable_if<is_arithmetic<T>::value >::type >
{
  static constexpr bool has_default_value = true;
  using value_type = T;
  value_type value = value_type();
};

HAS_ALIAS( value_type )
template< typename T >
struct Helper< T, typename std::enable_if< has_alias_value_type<T>::value &&
                                          ( is_arithmetic<typename T::value_type>::value) >::type >
{
  static constexpr bool has_default_value = true;
  using value_type = typename T::value_type;
  value_type value = value_type();
};

}

template< typename T >
using DefaultValue = wrapperDefaultValue::Helper<T>;


}
}


#endif /* CORECOMPONENTS_DATAREPOSITORY_DEFAULTVALUE_HPP_ */
