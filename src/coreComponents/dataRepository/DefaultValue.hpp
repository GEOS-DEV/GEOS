/*
 * WrapperDefaultValueHelper.hpp
 *
 *  Created on: Dec 17, 2018
 *      Author: settgast
 */

#ifndef CORECOMPONENTS_DATAREPOSITORY_DEFAULTVALUE_HPP_
#define CORECOMPONENTS_DATAREPOSITORY_DEFAULTVALUE_HPP_

#include "common/DataTypes.hpp"
#include "SFINAE_Macros.hpp"

namespace geosx
{
namespace dataRepository
{

/**
 * @namespace wrapperDefaultValue
 *
 * namespace to scope traits that are used for default values
 */
namespace wrapperDefaultValue
{

/**
 * @struct is_defaultable
 * @tparam T type to check
 * @brief trait to determine if type \p T should have a default value
 */
template< typename T >
struct is_defaultable
{
  /// attribute to set what type is able to contain a default value
  static constexpr bool value = std::is_same<T,int>::value ||
                                std::is_same<T,long int>::value ||
                                std::is_same<T,long long int>::value ||
                                std::is_same<T,unsigned int>::value ||
                                std::is_same<T,unsigned long int>::value ||
                                std::is_same<T,unsigned long long int>::value ||
                                std::is_floating_point<T>::value ||
                                std::is_same<T,string>::value ||
                                std::is_same<T,R1Tensor>::value ||
                                std::is_same<T,R2Tensor>::value ||
                                std::is_same<T,R2SymTensor>::value ;
};

/**
 * @struct Helper
 * @tparam T type to check
 * @tparam ENABLE template parameter for use with SFINAE
 *
 * default implementation of struct to return if a type \p T has a default value.
 */
template< typename T, typename ENABLE=void >
struct Helper
{
  /// attribute to indicate whether type \p T has a default value
  static constexpr bool has_default_value = false;
};

/**
 * @struct Helper
 * @tparam T type to check
 *
 * Specialization of Helper struct to return if a type \p T has a default
 * value. This specialization specifically tests the type itself. Contains
 * a member to hold a default value.
 */
template< typename T >
struct Helper< T, typename std::enable_if<is_defaultable<T>::value >::type >
{
  /// attribute to indicate whether type \p T has a default value
  static constexpr bool has_default_value = true;

  /// alias for the type T
  using value_type = T;

  /// a member to hold a default value for the type \p T
  value_type value = value_type();
};

HAS_ALIAS( value_type )
/**
 * @struct Helper
 * @tparam T type to check
 *
 * Specialization of Helper struct to return if a type \p T has a default
 * value. This specialization specifically tests the type has an alias
 * named "value_type" as is the case for stl containers and GEOSX
 * containers.
 */
template< typename T >
struct Helper< T, typename std::enable_if< has_alias_value_type<T>::value &&
                                          ( is_defaultable<typename T::value_type>::value) >::type >
{
  /// attribute to indicate whether type \p T has a default value
  static constexpr bool has_default_value = true;

  /// alias for the type T
  using value_type = typename T::value_type;

  /// a member to hold a default value for the type \p T
  value_type value = value_type();
};

}

/**
 * @tparam T the type to check
 * A templated alias to hold default values.
 */
template< typename T >
using DefaultValue = wrapperDefaultValue::Helper<T>;


}
}


#endif /* CORECOMPONENTS_DATAREPOSITORY_DEFAULTVALUE_HPP_ */
