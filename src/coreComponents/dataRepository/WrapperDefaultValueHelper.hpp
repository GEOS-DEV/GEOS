/*
 * WrapperDefaultValueHelper.hpp
 *
 *  Created on: Dec 17, 2018
 *      Author: settgast
 */

#ifndef CORECOMPONENTS_DATAREPOSITORY_WRAPPERDEFAULTVALUEHELPER_HPP_
#define CORECOMPONENTS_DATAREPOSITORY_WRAPPERDEFAULTVALUEHELPER_HPP_

#include "SFINAE_Macros.hpp"

namespace geosx
{
namespace dataRepository
{

namespace wrapperDefaultValue
{
template< typename T, typename ENABLE=void >
struct Helper
{
  static constexpr bool has_value_type = false;
};

template< typename T >
struct Helper< T, typename std::enable_if<std::is_integral<T>::value || std::is_floating_point<T>::value>::type >
{
  static constexpr bool has_value_type = true;
  using value_type = T;
  value_type m_default;
};

HAS_ALIAS( value_type )
template< typename T >
struct Helper< T, typename std::enable_if<has_alias_value_type<T>::value>::type >
{
  static constexpr bool has_value_type = true;
  using value_type = typename T::value_type;
  value_type m_default;
};

}


}
}


#endif /* CORECOMPONENTS_DATAREPOSITORY_WRAPPERDEFAULTVALUEHELPER_HPP_ */
