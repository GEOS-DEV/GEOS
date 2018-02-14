/*
 * static_if.hpp
 *
 *  Created on: Feb 6, 2018
 *      Author: settgast1
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CODINGUTILITIES_STATIC_IF_HPP_
#define SRC_COMPONENTS_CORE_SRC_CODINGUTILITIES_STATIC_IF_HPP_


namespace geosx
{

template<bool CONDITION>
struct static_if_wrapper
{
  template<typename LAMBDA_BODY>
  constexpr inline static void if_function(LAMBDA_BODY&&) { }
};

template<>
struct static_if_wrapper<true>
{
  template<typename LAMBDA_BODY>
  constexpr inline static void if_function(LAMBDA_BODY&& lambda)
  {
    //lambda();
  }
};

}

#define static_if( condition ) \
  geosx::static_if_wrapper<condition>::if_function( [&] () -> void

#endif /* SRC_COMPONENTS_CORE_SRC_CODINGUTILITIES_STATIC_IF_HPP_ */
