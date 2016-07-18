/*
 * StringUtilities.hpp
 *
 *  Created on: Nov 12, 2014
 *      Author: rrsettgast
 */

#ifndef STRINGUTILITIES_HPP_
#define STRINGUTILITIES_HPP_

#include <cxxabi.h>
#include <string>

namespace geosx
{
namespace stringutilities
{
inline std::string demangle( const std::string& name )
{

  int status = -4; // some arbitrary value to eliminate the compiler warning

  // enable c++11 by passing the flag -std=c++11 to g++
  std::unique_ptr<char, void (*)( void* )> res
  {
    abi::__cxa_demangle( name.c_str(), nullptr, nullptr, &status ),
    std::free
  };

  return ( status == 0 ) ? res.get() : name;
}

}
}

#endif /* STRINGUTILITIES_HPP_ */
