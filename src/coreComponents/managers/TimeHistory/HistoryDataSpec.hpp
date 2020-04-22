#ifndef GEOSX_HISTORY_DATA_SPEC_HPP_
#define GEOSX_HISTORY_DATA_SPEC_HPP_

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

  class HistoryMetadata
  {
  public:
    HistoryMetadata( const string & name, localIndex rank, size_t * dims, std::type_index type ) :
     m_name(name),
     m_rank(rank),
     m_dims(dims,dims+rank),
     m_type(type)
    {}
    HistoryMetadata( const string & name, size_t count, std::type_index type ) :
      m_name(name),
      m_rank(1),
      m_dims(&count,&count+1),
      m_type(type)
    {}
    const string & getName( ) const
    {
      return m_name;
    }
    std::type_index getType( ) const
    {
      return m_type;
    }
    size_t getTypeCount( ) const
    {
      size_t local_size = 1;
      for( size_t dim : m_dims )
      {
        local_size *= dim;
      }
      return local_size;
    }
    localIndex getRank( ) const
    {
      return m_rank;
    }
    const localIndex * getDims( ) const
    {
      return &m_dims[0];
    }
    localIndex getDimExtent(localIndex dim) const
    {
      return m_dims[dim];
    }
  private:
    std::string m_name;
    localIndex m_rank;
    std::vector<localIndex> m_dims;
    std::type_index m_type;
  };

  template < typename ARRAY_T >
  inline
  typename std::enable_if < is_array<ARRAY_T> && (ARRAY_T::ndim == 1) && can_history_io<typename ARRAY_T::value_type>, HistoryMetadata >::type
  getHistoryMetadata( string const & name, const ARRAY_T & arr, size_t size_overwrite = 0 )
  {
    size_t size = size_overwrite == 0 ? arr.size( ) : size_overwrite;
    return HistoryMetadata(name, 1, &size, std::type_index(typeid(typename ARRAY_T::value_type)));
  }

  template < typename ARRAY_T >
  inline
  typename std::enable_if < is_array<ARRAY_T> && (ARRAY_T::ndim > 1) && can_history_io<typename ARRAY_T::value_type>, HistoryMetadata >::type
  getHistoryMetadata( string const & name, ARRAY_T const & arr, size_t size_overwrite = 0 )
  {
    size_t sizes[2] = {size_overwrite == 0 ? integer_conversion<size_t>(arr.size( ) / arr.size(0)) : size_overwrite, integer_conversion<size_t>(arr.size( )) };
    return HistoryMetadata(name, 2, &sizes[0], std::type_index(typeid(typename ARRAY_T::value_type)));
  }

  template < typename T >
  inline typename std::enable_if < can_history_io< T >, HistoryMetadata >::type
  getHistoryMetadata( string const & name, const T & GEOSX_UNUSED_PARAM( type ), size_t size_overwrite = 0 )
  {
    size_t size = size_overwrite == 0 ? 1 : size_overwrite;
    return HistoryMetadata(name, size, std::type_index(typeid(T)));
  }

  template < typename T >
  inline typename std::enable_if < ! is_array<T> && ! can_history_io< T >, HistoryMetadata >::type
  getHistoryMetadata( string const & GEOSX_UNUSED_PARAM(name), const T & GEOSX_UNUSED_PARAM(type), size_t size_overwrite = 0 )
  {
    GEOSX_ERROR("Trying to use time history output on an unsupported type.");
    return HistoryMetadata("NULL", size_overwrite, std::type_index(typeid(NULL)));
  }

}

#endif
