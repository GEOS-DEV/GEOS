/*
 * DataTypes.hpp
 *
 *  Created on: Jun 16, 2016
 *      Author: rrsettgast
 */

#ifndef COMPONENTS_CORE_SRC_DATAREPOSITORY_DATATYPES_HPP_
#define COMPONENTS_CORE_SRC_DATAREPOSITORY_DATATYPES_HPP_

#include <cassert>
#include <cstdint>
#include <typeinfo>
#include "codingUtilities/macros.hpp"
#include <unordered_map>
#include <vector>
#include <string>
#include <typeindex>
#include<iostream>

using  int32     = std::int32_t;
using uint32     = std::uint32_t;
using  int64     = std::int64_t;
using uint64     = std::uint64_t;
using std_size_t = std::size_t;
using string     = std::string;

using real32 = float;
using real64 = double;

template< typename T >
using ptr = T*;

template< typename T >
using c_ptr = T const *;

template< typename T >
using array = std::vector<T>;

using int32_ptr        = ptr<int32>;
using int32_const_ptr  = c_ptr<int32>;

using uint32_ptr        = ptr<uint32>;
using uint32_const_ptr  = c_ptr<uint32>;

using int64_ptr        = ptr<int64>;
using int64_const_ptr  = c_ptr<int64>;

using uint64_ptr        = ptr<uint64>;
using uint64_const_ptr  = c_ptr<uint64>;

using real32_ptr        = ptr<real32>;
using real32_const_ptr  = c_ptr<real32>;

using real64_ptr        = ptr<real64>;
using real64_const_ptr  = c_ptr<real64>;



using int32_array        = array<int32>;
using int32_const_array  = array<int32 const>;

using uint32_array        = array<uint32>;
using uint32_const_array  = array<uint32 const>;

using int64_array        = array<int64>;
using int64_const_array  = array<int64 const>;

using uint64_array        = array<uint64>;
using uint64_const_array  = array<uint64 const>;

using real32_array        = array<real32>;
using real32_const_array  = array<real32 const>;

using real64_array        = array<real64>;
using real64_const_array  = array<real64 const>;

class rtTypes
{
public:

  static std::string typeNames( std::type_index const key )
  {
    const std::unordered_map<std::type_index, std::string> type_names =
    {
     {std::type_index(typeid(int32)), "int32"},
     {std::type_index(typeid(uint32)), "uint32"},
     {std::type_index(typeid(int64)), "int64"},
     {std::type_index(typeid(uint64)), "uint64"},
     {std::type_index(typeid(real32)), "real32"},
     {std::type_index(typeid(real64)), "real64"},
     {std::type_index(typeid(int32_array)), "int32_array"},
     {std::type_index(typeid(uint32_array)), "uint32_array"},
     {std::type_index(typeid(int64_array)), "int64_array"},
     {std::type_index(typeid(uint64_array)), "uint64_array"},
     {std::type_index(typeid(real32_array)), "real32_array"},
     {std::type_index(typeid(real64_array)), "real64_array"},
     {std::type_index(typeid(std_size_t)), "std_size_t"},
     {std::type_index(typeid(string)), "string"}
    };
    return type_names.at(key);
  }

  enum class TypeIDs
  {
    int32_id,
    uint32_id,
    int64_id,
    uint64_id,
    real32_id,
    real64_id,
    int32_array_id,
    uint32_array_id,
    int64_array_id,
    uint64_array_id,
    real32_array_id,
    real64_array_id,
    std_size_t_id,
    string_id,
    number_of_ids
  };


  template< typename LAMBDA >
  static auto ApplyTypeLambda( const TypeIDs type,
                               LAMBDA lambda )
  {
    assert( type>=TypeIDs::int32_id && type<TypeIDs::number_of_ids );
    switch( type )
    {
      case( TypeIDs::int32_id ):
      {
        return lambda( int32(1) );
        break;
      }
      case( TypeIDs::uint32_id ):
      {
        return lambda( uint32(1) );
        break;
      }
      case( TypeIDs::int64_id ):
      {
        return lambda( int64(1) );
        break;
      }
      case( TypeIDs::uint64_id ):
      {
        return lambda( uint64(1) );
        break;
      }
      case( TypeIDs::real32_id ):
      {
        return lambda( real32(1) );
        break;
      }
      case( TypeIDs::real64_id ):
      {
        return lambda( real64(1) );
        break;
      }
      case( TypeIDs::int32_array_id ):
      {
        return lambda( int32_array(1) );
        break;
      }
      case( TypeIDs::uint32_array_id ):
      {
        return lambda( uint32_array(1) );
        break;
      }
      case( TypeIDs::int64_array_id ):
      {
        return lambda( int64_array(1) );
        break;
      }
      case( TypeIDs::uint64_array_id ):
      {
        return lambda( uint64_array(1) );
        break;
      }
      case( TypeIDs::real32_array_id ):
      {
        return lambda( real32_array(1) );
        break;
      }
      case( TypeIDs::real64_array_id ):
      {
        return lambda( real64_array(1) );
        break;
      }
      case( TypeIDs::std_size_t_id ):
      {
        return lambda( std_size_t(1) );
        break;
      }
      case( TypeIDs::string_id ):
      {
        return lambda( string("") );
        break;
      }
      default:
      {
        std::cout<<LOCATION<<std::endl;
        throw std::exception();
      }
    }
  }

  /*
  template< typename LAMBDA, typename...ArgsF >
  static auto ApplyTypeLambda( const TypeIDs type,
                               LAMBDA lambda,
                               ArgsF&&... args ) -> decltype( lambda(args...) )
  {
    decltype(lambda(args...)) rval;
    switch( type )
    {
      case( TypeIDs::int32_id ):
      {
        rval = lambda(args...);
        break;
      }
      default:
      {
        std::cout<<LOCATION<<std::endl;
        throw std::exception();
      }
    }

    return rval;
  }
*/

};

#endif /* COMPONENTS_CORE_SRC_DATAREPOSITORY_DATATYPES_HPP_ */
