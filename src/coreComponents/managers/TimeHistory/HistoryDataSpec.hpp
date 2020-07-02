#ifndef GEOSX_HISTORY_DATA_SPEC_HPP_
#define GEOSX_HISTORY_DATA_SPEC_HPP_

#include "codingUtilities/traits.hpp"
#include "common/DataTypes.hpp"
#include "LvArray/src/Array.hpp"

namespace geosx
{
  template < typename T >
  constexpr bool can_history_io = std::is_same<std::remove_reference_t<std::remove_const_t<T>>, char>::value ||
    std::is_same<std::remove_reference_t<std::remove_const_t<T>>, signed char>::value ||
    std::is_same<std::remove_reference_t<std::remove_const_t<T>>, real32>::value ||
    std::is_same<std::remove_reference_t<std::remove_const_t<T>>, real64>::value ||
    std::is_same<std::remove_reference_t<std::remove_const_t<T>>, integer>::value ||
    std::is_same<std::remove_reference_t<std::remove_const_t<T>>, int>::value ||
    std::is_same<std::remove_reference_t<std::remove_const_t<T>>, double>::value ||
    std::is_same<std::remove_reference_t<std::remove_const_t<T>>, localIndex>::value ||
    std::is_same<std::remove_reference_t<std::remove_const_t<T>>, globalIndex>::value;

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
    void setName( const string & name )
    {
      m_name = name;
    }
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

  template < typename T >
  constexpr bool is_array_type = traits::is_array_view< T > || traits::is_array< T >;

  template < typename T >
  constexpr bool is_sorted_array_type = traits::is_sorted_array_view< T > || traits::is_sorted_array< T >;

  template < typename T >
  constexpr bool can_history_io_container = is_array_type< T > || is_sorted_array_type< T >;

  template < typename ARRAY_T >
  inline
  HistoryMetadata getFlatArrayMetadata( string const & name, const ARRAY_T & arr, localIndex size_override = -1)
  {
    localIndex size = size_override < 0 ? arr.size( ) : size_override;
    size_t sz = LvArray::integerConversion<size_t>(size);
    return HistoryMetadata(name, 1, &sz, std::type_index(typeid(typename ARRAY_T::value_type)));
  }

  template < typename ARRAY_T >
  inline
  typename std::enable_if < is_array_type< ARRAY_T > && (ARRAY_T::ndim == 1) && can_history_io<typename ARRAY_T::value_type>, HistoryMetadata >::type
  getHistoryMetadata( string const & name, const ARRAY_T & arr, localIndex size_override = -1 )
  {
    return getFlatArrayMetadata< ARRAY_T >( name, arr, size_override );
  }

  template < typename ARRAY_T >
  inline
  typename std::enable_if < is_sorted_array_type< ARRAY_T > && can_history_io< typename ARRAY_T::value_type >, HistoryMetadata >::type
  getHistoryMetadata( string const & name, const ARRAY_T & arr, localIndex size_override = -1 )
  {
    return getFlatArrayMetadata< ARRAY_T >( name, arr, size_override );
  }

  template < typename ARRAY_T >
  inline
  typename std::enable_if < ( is_array_type< ARRAY_T >) && (ARRAY_T::ndim > 1) && can_history_io<typename ARRAY_T::value_type>, HistoryMetadata >::type
  getHistoryMetadata( string const & name, ARRAY_T const & arr, localIndex size_override = -1 )
  {
    size_t sizes[2] = {size_override < 0 ? LvArray::integerConversion<size_t>(arr.size( ) / arr.size(0)) : size_override, LvArray::integerConversion<size_t>(arr.size( )) };
    return HistoryMetadata(name, 2, &sizes[0], std::type_index(typeid(typename ARRAY_T::value_type)));
  }

  template < typename T >
  inline typename std::enable_if < can_history_io< T >, HistoryMetadata >::type
  getHistoryMetadata( string const & name, const T & GEOSX_UNUSED_PARAM( type ), localIndex size_override = -1 )
  {
    size_t size = size_override < 0 ? 0 : LvArray::integerConversion<size_t>(size_override);
    return HistoryMetadata(name, size, std::type_index(typeid(T)));
  }

  template < typename ARRAY_T >
  inline typename std::enable_if < ( can_history_io_container< ARRAY_T > ) && !can_history_io< typename ARRAY_T::value_type >, HistoryMetadata >::type
  getHistoryMetadata( string const & GEOSX_UNUSED_PARAM(name), const ARRAY_T & GEOSX_UNUSED_PARAM(type), localIndex size_override )
  {
    GEOSX_ERROR("Trying to use time history output on an array containing an unsupported type.");
    GEOSX_UNUSED_VAR( size_override );
    return HistoryMetadata("NULL", 0, std::type_index(typeid(NULL)));
  }

  template < typename T >
  inline typename std::enable_if < ! ( can_history_io_container< T > ) && ! can_history_io< T >, HistoryMetadata >::type
  getHistoryMetadata( string const & GEOSX_UNUSED_PARAM(name), const T & GEOSX_UNUSED_PARAM(type), localIndex size_override )
  {
    GEOSX_ERROR("Trying to use time history output on an unsupported type.");
    GEOSX_UNUSED_VAR( size_override );
    return HistoryMetadata("NULL", 0, std::type_index(typeid(NULL)));
  }


}

#endif
