#ifndef GEOSX_HISTORY_DATA_SPEC_HPP_
#define GEOSX_HISTORY_DATA_SPEC_HPP_

#include "dataRepository/Group.hpp"
#include "codingUtilities/traits.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{
  using namespace traits;

  template < typename T >
  constexpr bool can_history_io = std::is_same<std::remove_const_t<T>, char>::value ||
                            std::is_same<std::remove_const_t<T>, signed char>::value ||
                            std::is_same<std::remove_const_t<T>, real32>::value ||
                            std::is_same<std::remove_const_t<T>, real64>::value ||
                            std::is_same<std::remove_const_t<T>, integer>::value ||
                            std::is_same<std::remove_const_t<T>, localIndex>::value ||
                            std::is_same<std::remove_const_t<T>, globalIndex>::value;

  class DataSpec
  {
  public:
    virtual void Append(localIndex const num_items,
                            localIndex const units_per_item,
                            size_t const unit_byte_size,
                            std::type_info const & type,
                            string const & name) = 0;
  };

  template < typename ARRAY_T >
  inline
  typename std::enable_if < is_array<ARRAY_T> && (ARRAY_T::ndim == 1) && can_history_io<typename ARRAY_T::value_type>, void >::type
  AppendArraySpec( DataSpec & spec, ARRAY_T const & arr, string const & name, localIndex vals_per_col = 1 )
  {
    spec.Append(arr.size( ), vals_per_col, sizeof(typename ARRAY_T::value_type), typeid(typename ARRAY_T::value_type), name);
  }

  template < typename ARRAY_T >
  inline
  typename std::enable_if < is_array<ARRAY_T> && (ARRAY_T::ndim > 1) && can_history_io<typename ARRAY_T::value_type>, void >::type
  AppendArraySpec( DataSpec & spec, ARRAY_T const & arr, string const & name)
  {
    spec.Append(arr.size( ) / arr.size( 0 ), arr.size( 0 ), sizeof(typename ARRAY_T::value_type), typeid(typename ARRAY_T::value_type), name);
  }

  template < typename ARRAY_T >
  inline
  typename std::enable_if < is_array<ARRAY_T> && (ARRAY_T::ndim == 1) && can_history_io<typename ARRAY_T::value_type>, void >::type
  AppendArrayIndicesSpec( DataSpec & spec, ARRAY_T const & GEOSX_UNUSED_PARAM( arr ), string const & name, localIndex const num_indices, localIndex vals_per_col = 1 )
  {
    spec.Append(num_indices, vals_per_col, sizeof(typename ARRAY_T::value_type), typeid(typename ARRAY_T::value_type), name);
  }

  template < typename ARRAY_T >
  inline
  typename std::enable_if < is_array<ARRAY_T> && (ARRAY_T::ndim > 1) && can_history_io<typename ARRAY_T::value_type>, void >::type
  AppendArrayIndicesSpec( DataSpec & spec, ARRAY_T const & arr, string const & name, localIndex const num_indices )
  {
    spec.Append(num_indices, arr.size( 0 ), sizeof(typename ARRAY_T::value_type), typeid(typename ARRAY_T::value_type), name);
  }

}

#endif