/*
 * DataTypes.hpp
 *
 *  Created on: Jun 16, 2016
 *      Author: rrsettgast
 */

#ifndef COMPONENTS_CORE_SRC_DATAREPOSITORY_DATATYPES_HPP_
#define COMPONENTS_CORE_SRC_DATAREPOSITORY_DATATYPES_HPP_

#include <cstdint>
#include <typeinfo>
#include "codingUtilities/macros.hpp"
#include <unordered_map>
#include <string>
#include <typeindex>
#include<iostream>

using  int32 = std::int32_t;
using uint32 = std::uint32_t;
using  int64 = std::int64_t;
using uint64 = std::uint64_t;

using real32 = float;
using real64 = double;

template< typename T >
using ptr = T*;

template< typename T >
using c_ptr = T const *;

using real32_ptr       = ptr<real32>;
using real32_const_ptr = ptr<real32 const>;

using real64_ptr       = ptr<real64>;
using real64_const_ptr = ptr<real64 const>;

using real64_array        = ptr<real64>;
using real64_const_array  = ptr<real64 const>;

class rtTypes
{
public:

  static std::string typeNames( std::type_index const key )
  {
    static std::unordered_map<std::type_index, std::string> type_names =
    {
     {std::type_index(typeid(int32)), "int32"},
     {std::type_index(typeid(uint32)), "uint32"},
     {std::type_index(typeid(int64)), "int64"},
     {std::type_index(typeid(uint64)), "uint64"},
     {std::type_index(typeid(real32)), "real32"},
     {std::type_index(typeid(real64)), "real64"}
    };
    return type_names[key];
  }

  enum class TypeIDs
  {
    int32_id,
    uint32_id,
    int64_id,
    uint64_id,
    real32_id,
    real64_id
  };


  template< typename LAMBDA >
  static auto ApplyTypeLambda( const TypeIDs type,
                               LAMBDA lambda )
  {
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
