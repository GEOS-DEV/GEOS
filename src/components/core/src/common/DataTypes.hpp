/*
 * DataTypes.hpp
 *
 *  Created on: Jun 16, 2016
 *      Author: rrsettgast
 */

#ifndef DATATYPES_HPP
#define DATATYPES_HPP

#include <cassert>
#include <cstdint>
#include <iostream>
#include <set>
#include <string>
#include <typeindex>
#include <typeinfo>
#include <map>
#include <unordered_map>
#include <vector>

#include "common/GeosxConfig.hpp"

#include "Macros.hpp"

#include "ManagedArray.hpp"

//#include "legacy/ArrayT/ArrayT.h"
#include "math/TensorT/TensorT.h"

#ifndef CONTAINERARRAY_RETURN_PTR
#define CONTAINERARRAY_RETURN_PTR 1
#endif


namespace geosx
{

// underlying types not for general use!!
using int32 = std::int32_t;
using int64 = std::int64_t;
using uint32 = std::uint32_t;
using uint64 = std::uint64_t;

using integer     = std::int32_t;
using localIndex  = std::int_fast32_t;
using globalIndex = std::int64_t;

using string      = std::string;



using real32 = float;
using real64 = double;
using real   = double;

template< typename T >
using ptr = T*;

template< typename T >
using c_ptr = T const *;


using integer_ptr        = ptr<integer>;
using integer_const_ptr  = c_ptr<integer>;

using localIndex_ptr         = ptr<localIndex>;
using localIndex_const_ptr   = c_ptr<localIndex>;

using real32_ptr        = ptr<real32>;
using real32_const_ptr  = c_ptr<real32>;

using real64_ptr        = ptr<real64>;
using real64_const_ptr  = c_ptr<real64>;

//***** BEGIN LEGACY TYPEDEFS *****

using realT    = double;


template< typename T >
//using array = std::vector<T>;
//using array = Array1dT<T>;
using array = multidimensionalArray::ManagedArray<T,1,localIndex>;


template< typename T >
using set = std::set<T>;

template< typename TKEY, typename TVAL >
using map = std::map<TKEY,TVAL>;

template< typename TKEY, typename TVAL >
using unordered_map = std::unordered_map<TKEY,TVAL>;


//***** END LEGACY TYPEDEFS *****



using integer_array        = array<integer>;
//using integer_array = multidimensionalArray::ManagedArray<integer,1,localIndex>;
using integer_const_array  = array<integer const>;

using real32_array        = array<real32>;
using real32_const_array  = array<real32 const>;

using real64_array        = array<real64>;
using real64_const_array  = array<real64 const>;

using string_array        = array<string>;
using string_const_array  = array<string const>;

using localIndex_array        = array<localIndex>;
using localIndex_const_array  = array<localIndex const>;

using globalIndex_array        = array<globalIndex>;
using globalIndex_const_array  = array<globalIndex const>;



using integer_set        = set<integer>;
using integer_const_set  = set<integer const>;

//using uint32_set        = set<uint32>;
//using uint32_const_set  = set<uint32 const>;
//
//using int64_set        = set<int64>;
//using int64_const_set  = set<int64 const>;
//
//using uint64_set        = set<uint64>;
//using uint64_const_set  = set<uint64 const>;

using real32_set        = set<real32>;
using real32_const_set  = set<real32 const>;

using real64_set        = set<real64>;
using real64_const_set  = set<real64 const>;

using string_set        = set<string>;
using string_const_set  = set<string const>;

using localIndex_set        = set<localIndex>;
using localIndex_const_set  = set<localIndex const>;

using globalIndex_set        = set<globalIndex>;
using globalIndex_const_set  = set<globalIndex const>;



//***** BEGIN LEGACY TYPEDEFS *****
//using rArray1d = Array1dT<real64>;
//using iArray1d = Array1dT<integer>;
//using lArray1d = Array1dT<localIndex>;
//using gArray1d = Array1dT<globalIndex>;

//typedef Array1dT<std::string> sArray1d;

typedef std::set<localIndex> lSet;
typedef std::set<globalIndex> gSet;

typedef int FieldKey;

template< typename T >
using Array1dT = multidimensionalArray::ManagedArray<T,1,localIndex>;

template< typename T >
using Array2dT = multidimensionalArray::ManagedArray<T,2,localIndex>;

typedef multidimensionalArray::ManagedArray<localIndex,2,localIndex> lArray2d;
typedef Array1dT<std::pair<int,localIndex> > pArray1d;
typedef std::set<std::pair<int,localIndex> > pSet;

using rArray1d = array<real64>;
using iArray1d = array<integer>;
using sArray1d = array<string>;
using r1_array = array<R1Tensor>;
using r2_array = array<R2Tensor>;
using r2Sym_array = array<R2SymTensor>;

//using mapPair = std::pair<integer, localIndex>;
using mapPair_array = std::pair<integer_array, integer_array>;


constexpr static auto GLOBALINDEX_MAX = std::numeric_limits<globalIndex>::max();

//***** END LEGACY TYPEDEFS *****

class rtTypes
{
public:

  static std::string typeNames( std::type_index const key )
  {
    const std::unordered_map<std::type_index, std::string> type_names =
    {
      {std::type_index(typeid(integer)), "integer"},
      {std::type_index(typeid(real32)), "real32"},
      {std::type_index(typeid(real64)), "real64"},
      {std::type_index(typeid(localIndex)), "localIndex"},
      {std::type_index(typeid(globalIndex)), "globalIndex"},
      {std::type_index(typeid(R1Tensor)), "r1Tensor"},
      {std::type_index(typeid(R2Tensor)), "r2Tensor"},
      {std::type_index(typeid(R2SymTensor)), "r2SymTensor"},
      {std::type_index(typeid(integer_array)), "integer_array"},
      {std::type_index(typeid(real32_array)), "real32_array"},
      {std::type_index(typeid(real64_array)), "real64_array"},
      {std::type_index(typeid(localIndex_array)), "localIndex_array"},
      {std::type_index(typeid(globalIndex_array)), "globalIndex_array"},
      {std::type_index(typeid(r1_array)), "r1_array"},
      {std::type_index(typeid(r2_array)), "r2_array"},
      {std::type_index(typeid(r2Sym_array)), "r2Sym_array"},
      {std::type_index(typeid(string)), "string"},
      {std::type_index(typeid(mapPair_array)), "mapPair_array"}
    };
    return type_names.at(key);
  }



  enum class TypeIDs
  {
    integer_id,
    localIndex_id,
    globalIndex_id,
    real32_id,
    real64_id,
    r1Tensor_id,
    r2Tensor_id,
    r2SymTensor_id,
    integer_array_id,
    localIndex_array_id,
    globalIndex_array_id,
    real32_array_id,
    real64_array_id,
    r1_array_id,
    r2_array_id,
    r2Sym_array_id,
    std_size_t_id,
    string_id,
    string_array_id,
    mapPair_array_id,
    none_id
  };

  static TypeIDs typeID( string const & name )
  {
    const std::unordered_map<string,TypeIDs> type_names =
    {
      { "integer",        TypeIDs::integer_id },
      { "localIndex",   TypeIDs::localIndex_id },
      { "globalIndex",  TypeIDs::globalIndex_id },
      { "real32",       TypeIDs::real32_id },
      { "real64",       TypeIDs::real64_id },
      { "R1Tensor",     TypeIDs::r1Tensor_id },
      { "R2Tensor",     TypeIDs::r2Tensor_id },
      { "R2SymTensor",  TypeIDs::r2SymTensor_id },
      { "integer_array",  TypeIDs::integer_array_id },
      { "localIndex_array",   TypeIDs::localIndex_array_id },
      { "globalIndex_array",  TypeIDs::globalIndex_array_id },
      { "real32_array", TypeIDs::real32_array_id },
      { "real64_array", TypeIDs::real64_array_id },
      { "r1_array",     TypeIDs::r1_array_id },
      { "r2_array",     TypeIDs::r2_array_id },
      { "r2Sym_array",  TypeIDs::r2Sym_array_id },
      { "std_size_t",   TypeIDs::std_size_t_id },
      { "string",       TypeIDs::string_id },
      { "string_array", TypeIDs::string_array_id },
      { "mapPair_array",      TypeIDs::mapPair_array_id },
      { "",             TypeIDs::none_id }
    };
    return type_names.at(name);
  }

  static TypeIDs typeID( std::type_index typeIndex )
  {
    const std::unordered_map<std::type_index,TypeIDs> type_names =
    {
      { std::type_index(typeid(integer)),        TypeIDs::integer_id },
      { std::type_index(typeid(localIndex)),   TypeIDs::real32_id },
      { std::type_index(typeid(globalIndex)),  TypeIDs::real64_id },
      { std::type_index(typeid(real32)),       TypeIDs::real32_id },
      { std::type_index(typeid(real64)),       TypeIDs::real64_id },
      { std::type_index(typeid(R1Tensor)),     TypeIDs::r1Tensor_id },
      { std::type_index(typeid(R2Tensor)),     TypeIDs::r2Tensor_id },
      { std::type_index(typeid(R2SymTensor)),  TypeIDs::r2SymTensor_id },
      { std::type_index(typeid(integer_array)),  TypeIDs::integer_array_id },
      { std::type_index(typeid(localIndex_array)),  TypeIDs::localIndex_array_id },
      { std::type_index(typeid(globalIndex_array)),  TypeIDs::globalIndex_array_id },
      { std::type_index(typeid(real32_array)), TypeIDs::real32_array_id },
      { std::type_index(typeid(real64_array)), TypeIDs::real64_array_id },
      { std::type_index(typeid(r1_array)),     TypeIDs::r1_array_id },
      { std::type_index(typeid(r2_array)),     TypeIDs::r2_array_id },
      { std::type_index(typeid(r2Sym_array)),  TypeIDs::r2Sym_array_id },
      { std::type_index(typeid(string)),       TypeIDs::string_id },
      { std::type_index(typeid(string_array)), TypeIDs::string_array_id },
      { std::type_index(typeid(mapPair_array)),TypeIDs::mapPair_array_id }
    };
    return type_names.at(typeIndex);
  }



  // Matching regex for data types in xml
  class typeRegex
  {
  private:
    std::string ru = "[0-9]*";
    std::string ri = "[+-]?[0-9]*";    
    std::string rr = "[0-9]*\\.?([0-9]*)?[eE]?[-+]?([0-9]*)?";
    std::string rs = "[a-zA-Z0-9_,\\(\\)+-/\\*]*";
    std::string r1 = rr + ",? " + rr + ",? " + rr;
    std::string r2 = rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr;
    std::string r2s = rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr;
    
    std::unordered_map<std::string, std::string> regexMap = 
    {
      {"integer", ri},
      {"real32", rr},
      {"real64", rr},
      {"R1Tensor", r1},
      {"R2Tensor", r2},
      {"R2SymTensor", r2s},
      {"integer_array", "((" + ri + ",? )*)?" + ri},
      {"real32_array", "((" + rr + ",? )*)?" + rr},
      {"real64_array", "((" + rr + ",? )*)?" + rr},
      {"r1_array", "((" + r1 + "; )*)?" + r1},
      {"r2_array", "((" + r2 + "; )*)?" + r2},
      {"r2Sym_array", "((" + r2s + "; )*)?" + r2s},
      {"std_size_t", ru},
      {"string", rs},
      {"string_array", "((" + rs + ",? )*)?" + rs},
      {"mapPair", rs},
      {"mapPair_array", "((" + rs + ",? )*)?" + rs}
    };

  public:
    std::unordered_map<std::string, std::string>::iterator begin(){return regexMap.begin();}
    std::unordered_map<std::string, std::string>::iterator end(){return regexMap.end();}
    std::unordered_map<std::string, std::string>::const_iterator begin() const {return regexMap.begin();}
    std::unordered_map<std::string, std::string>::const_iterator end() const {return regexMap.end();}
  };



  template< typename LAMBDA >
  static auto ApplyIntrinsicTypeLambda1( const TypeIDs type,
                                         LAMBDA lambda )
  {
    switch( type )
    {
    case ( TypeIDs::integer_id ):
    {
      return lambda( integer(1) );
      break;
    }
    case ( TypeIDs::real32_id ):
    {
      return lambda( real32(1) );
      break;
    }
    case ( TypeIDs::real64_id ):
    {
      return lambda( real64(1) );
      break;
    }
    case ( TypeIDs::r1Tensor_id ):
    {
      return lambda( R1Tensor() );
      break;
    }
    case ( TypeIDs::r2Tensor_id ):
    {
      return lambda( R2Tensor() );
      break;
    }
    case ( TypeIDs::r2SymTensor_id ):
    {
      return lambda( R2SymTensor() );
      break;
    }
    case ( TypeIDs::integer_array_id ):
    {
      return lambda( integer_array(1) );
      break;
    }
    case ( TypeIDs::real32_array_id ):
    {
      return lambda( real32_array(1) );
      break;
    }
    case ( TypeIDs::real64_array_id ):
    {
      return lambda( real64_array(1) );
      break;
    }
    case ( TypeIDs::string_id ):
    {
      return lambda( string("") );
      break;
    }
    default:
    {
      std::cout<<LOCATION<<std::endl;
      assert( false );
    }
    }
  }



  template< typename LAMBDA >
  static auto ApplyArrayTypeLambda1( const TypeIDs type,
                                         LAMBDA lambda )
  {
    switch( type )
    {
    case ( TypeIDs::integer_array_id ):
    {
      return lambda( integer_array(1) );
      break;
    }
    case ( TypeIDs::localIndex_array_id ):
    {
      return lambda( localIndex_array(1) );
      break;
    }
    case ( TypeIDs::globalIndex_array_id ):
    {
      return lambda( globalIndex_array(1) );
      break;
    }
    case ( TypeIDs::real32_array_id ):
    {
      return lambda( real32_array(1) );
      break;
    }
    case ( TypeIDs::real64_array_id ):
    {
      return lambda( real64_array(1) );
      break;
    }
    case ( TypeIDs::r1_array_id ):
    {
      return lambda( r1_array(1) );
      break;
    }
    case ( TypeIDs::r2_array_id ):
    {
      return lambda( r2_array(1) );
      break;
    }
    case ( TypeIDs::r2Sym_array_id ):
    {
      return lambda( r2Sym_array(1) );
      break;
    }

    default:
    {
      std::cout<<LOCATION<<std::endl;
      assert( false );
    }
    }
  }



  template< typename LAMBDA >
  static auto ApplyArrayTypeLambda2( const TypeIDs type,
                                         LAMBDA lambda )
  {
    switch( type )
    {
    case ( TypeIDs::integer_array_id ):
    {
      return lambda( integer_array(1), integer(1) );
      break;
    }
    case ( TypeIDs::localIndex_array_id ):
    {
      return lambda( localIndex_array(1), localIndex(1) );
      break;
    }
    case ( TypeIDs::globalIndex_array_id ):
    {
      return lambda( globalIndex_array(1), globalIndex() );
      break;
    }
    case ( TypeIDs::real32_array_id ):
    {
      return lambda( real32_array(1), real32(1) );
      break;
    }
    case ( TypeIDs::real64_array_id ):
    {
      return lambda( real64_array(1), real64(1) );
      break;
    }
    case ( TypeIDs::r1_array_id ):
    {
      return lambda( r1_array(1), R1Tensor() );
      break;
    }
    case ( TypeIDs::r2_array_id ):
    {
      return lambda( r2_array(1), R2Tensor() );
      break;
    }
    case ( TypeIDs::r2Sym_array_id ):
    {
      return lambda( r2Sym_array(1), R2SymTensor()  );
      break;
    }

    default:
    {
      std::cout<<LOCATION<<std::endl;
      assert( false );
    }
    }
  }

  template< typename LAMBDA >
  static auto ApplyTypeLambda1( const TypeIDs type,
                               LAMBDA lambda )
  {
    switch( type )
    {
    case ( TypeIDs::integer_id ):
    {
      return lambda( integer(1) );
      break;
    }
    case ( TypeIDs::real32_id ):
    {
      return lambda( real32(1) );
      break;
    }
    case ( TypeIDs::real64_id ):
    {
      return lambda( real64(1) );
      break;
    }
    case ( TypeIDs::r1Tensor_id ):
    {
      return lambda( R1Tensor() );
      break;
    }
    case ( TypeIDs::r2Tensor_id ):
    {
      return lambda( R2Tensor() );
      break;
    }
    case ( TypeIDs::r2SymTensor_id ):
    {
      return lambda( R2SymTensor() );
      break;
    }
    case ( TypeIDs::integer_array_id ):
    {
      return lambda( integer_array(1) );
      break;
    }
    case ( TypeIDs::localIndex_array_id ):
    {
      return lambda( localIndex_array(1) );
      break;
    }
    case ( TypeIDs::globalIndex_array_id ):
    {
      return lambda( globalIndex_array(1) );
      break;
    }
    case ( TypeIDs::real32_array_id ):
    {
      return lambda( real32_array(1) );
      break;
    }
    case ( TypeIDs::real64_array_id ):
    {
      return lambda( real64_array(1) );
      break;
    }
    case ( TypeIDs::r1_array_id ):
    {
      return lambda( r1_array(1) );
      break;
    }
    case ( TypeIDs::r2_array_id ):
    {
      return lambda( r2_array(1) );
      break;
    }
    case ( TypeIDs::r2Sym_array_id ):
    {
      return lambda( r2Sym_array(1) );
      break;
    }
    case ( TypeIDs::string_id ):
    {
      return lambda( string("") );
      break;
    }
    case ( TypeIDs::string_array_id ):
    {
      return lambda( string_array(1) );
      break;
    }
    case ( TypeIDs::mapPair_array_id ):
    {
      return lambda( mapPair_array() );
      break;
    }
    default:
    {
      std::cout<<LOCATION<<std::endl;
      assert( false );
    }
    }
  }


  template< typename LAMBDA >
  static auto ApplyTypeLambda2( const TypeIDs type,
                               LAMBDA lambda )
  {
    switch( type )
    {
    case ( TypeIDs::integer_id ):
    {
      return lambda( integer(1), integer(1) );
      break;
    }
    case ( TypeIDs::real32_id ):
    {
      return lambda( real32(1), real32(1) );
      break;
    }
    case ( TypeIDs::real64_id ):
    {
      return lambda( real64(1), real64(1) );
      break;
    }
    case ( TypeIDs::r1Tensor_id ):
    {
      return lambda( R1Tensor(), R1Tensor() );
      break;
    }
    case ( TypeIDs::r2Tensor_id ):
    {
      return lambda( R2Tensor(), R2Tensor() );
      break;
    }
    case ( TypeIDs::r2SymTensor_id ):
    {
      return lambda( R2SymTensor(), R2SymTensor() );
      break;
    }
    case ( TypeIDs::integer_array_id ):
    {
      return lambda( integer_array(1), integer(1) );
      break;
    }
    case ( TypeIDs::localIndex_array_id ):
    {
      return lambda( localIndex_array(1), localIndex(1) );
      break;
    }
    case ( TypeIDs::globalIndex_array_id ):
    {
      return lambda( globalIndex_array(1), globalIndex(1) );
      break;
    }
    case ( TypeIDs::real32_array_id ):
    {
      return lambda( real32_array(1), real32(1) );
      break;
    }
    case ( TypeIDs::real64_array_id ):
    {
      return lambda( real64_array(1), real64(1) );
      break;
    }
    case ( TypeIDs::r1_array_id ):
    {
      return lambda( r1_array(1), R1Tensor() );
      break;
    }
    case ( TypeIDs::r2_array_id ):
    {
      return lambda( r2_array(1), R2Tensor() );
      break;
    }
    case ( TypeIDs::r2Sym_array_id ):
    {
      return lambda( r2Sym_array(1), R2SymTensor() );
      break;
    }
    case ( TypeIDs::string_id ):
    {
      return lambda( string(""), string("") );
      break;
    }
    case ( TypeIDs::string_array_id ):
    {
      return lambda( string_array(1), string("") );
      break;
    }
//    case ( TypeIDs::mapPair_array_id ):
//    {
//      return lambda( mapPair_array(1), mapPair({}) );
//      break;
//    }

    default:
    {
      std::cout<<LOCATION<<std::endl;
      assert( false );
    }
    }
  }



  template< typename LAMBDA >
  static auto ApplyIntrinsicTypeLambda2( const TypeIDs type,
                                         LAMBDA lambda )
  {
    switch( type )
    {
    case ( TypeIDs::integer_id ):
    {
      return lambda( integer(1), integer(1) );
      break;
    }
    case ( TypeIDs::real32_id ):
    {
      return lambda( real32(1), real32(1) );
      break;
    }
    case ( TypeIDs::real64_id ):
    {
      return lambda( real64(1), real64(1) );
      break;
    }
    case ( TypeIDs::r1Tensor_id ):
    {
      return lambda( R1Tensor(), R1Tensor() );
      break;
    }
    case ( TypeIDs::r2Tensor_id ):
    {
      return lambda( R2Tensor(), R2Tensor() );
      break;
    }
    case ( TypeIDs::r2SymTensor_id ):
    {
      return lambda( R2SymTensor(), R2SymTensor() );
      break;
    }
    case ( TypeIDs::integer_array_id ):
    {
      return lambda( integer_array(1), integer(1) );
      break;
    }
    case ( TypeIDs::real32_array_id ):
    {
      return lambda( real32_array(1), real32(1) );
      break;
    }
    case ( TypeIDs::real64_array_id ):
    {
      return lambda( real64_array(1), real64(1) );
      break;
    }
    case ( TypeIDs::string_id ):
    {
      return lambda( string(""), string("") );
      break;
    }
    case ( TypeIDs::string_array_id ):
    {
      return lambda( string_array(1), string("") );
      break;
    }
//    case ( TypeIDs::mapPair_array_id ):
//    {
//      return lambda( mapPair_array(1), mapPair({}) );
//      break;
//    }
    default:
    {
      std::cout<<LOCATION<<std::endl;
      assert( false );
    }
    }
  }







  inline static void equate( R1Tensor & lhs, integer const component, real64 const & rhs )
  {
    lhs[component] = rhs;
  }

  template< typename TLHS, typename TRHS >
  inline static void equate( TLHS & lhs, integer const component, TRHS const & rhs )
  {
    lhs = rhs;
  }


  inline static void add( R1Tensor & lhs, integer const component, real64 const & rhs )
  {
    lhs[component] += rhs;
  }

  template< typename TLHS, typename TRHS >
  inline static void add( TLHS & lhs, integer const component, TRHS const & rhs )
  {
    lhs += rhs;
  }



  inline static real64 value( R1Tensor & lhs, integer const component )
  {
    return lhs[component];
  }

  inline static real64 value( R2Tensor & lhs, integer const component )
  {
    return lhs.Data()[component];
  }

  inline static real64 value( R2SymTensor & lhs, integer const component )
  {
    return lhs.Data()[component];
  }

  template< typename TLHS >
  inline static TLHS value( TLHS & lhs, integer const component )
  {
    return lhs;
  }


  struct equateValue
  {
    inline static void f( R1Tensor & lhs, integer const component, real64 const & rhs )
    {
      lhs[component] = rhs;
    }

    template< typename T >
    inline static void f( R1Tensor & lhs, integer const component, T const & rhs )
    {
      lhs[component] = rhs;
    }

    template< typename TLHS, typename TRHS >
    inline static void f( TLHS & lhs, integer const component, TRHS const & rhs )
    {
      lhs = rhs;
    }

  };

  struct addValue
  {
    inline static void f( R1Tensor & lhs, integer const component, real64 const & rhs )
    {
      lhs[component] += rhs;
    }

    inline static void f( R2Tensor & lhs, integer const component, real64 const & rhs )
    {
      lhs.Data()[component] += rhs;
    }

    inline static void f( R2SymTensor & lhs, integer const component, real64 const & rhs )
    {
      lhs.Data()[component] += rhs;
    }

    template< typename TLHS, typename TRHS >
    inline static void f( TLHS & lhs, integer const component, TRHS const & rhs )
    {
      lhs += rhs;
    }

  };


  enum class operationType
  {
    add,
    multiply,
    equate
  };
};

}



#endif /* COMPONENTS_CORE_SRC_DATAREPOSITORY_DATATYPES_HPP_ */
