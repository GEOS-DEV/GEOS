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
  ArrayMetadata( string const & name, const ARRAY_T & arr )
  {
    size_t size = arr.size( );
    return HistoryMetadata(name, 1, &size, std::type_index(typeid(typename ARRAY_T::value_type)));
  }

  template < typename ARRAY_T >
  inline
  typename std::enable_if < is_array<ARRAY_T> && (ARRAY_T::ndim > 1) && can_history_io<typename ARRAY_T::value_type>, HistoryMetadata >::type
  ArrayMetadata( string const & name, ARRAY_T const & arr )
  {
    size_t sizes[2] = {integer_conversion<size_t>(arr.size( ) / arr.size(0)), integer_conversion<size_t>(arr.size( )) };
    return HistoryMetadata(name, 2, &sizes[0], std::type_index(typeid(typename ARRAY_T::value_type)));
  }

  template < typename ARRAY_T >
  inline
  typename std::enable_if < is_array<ARRAY_T> && (ARRAY_T::ndim == 1) && can_history_io<typename ARRAY_T::value_type>, HistoryMetadata >::type
  ArrayIndicesMetadata( ARRAY_T const & GEOSX_UNUSED_PARAM( arr ), string const & name, localIndex const num_indices )
  {
    size_t size = integer_conversion<size_t>(num_indices);
    return HistoryMetadata(name, 1, &size, std::type_index(typeid(typename ARRAY_T::value_type)));
  }

  template < typename ARRAY_T >
  inline
  typename std::enable_if < is_array<ARRAY_T> && (ARRAY_T::ndim > 1) && can_history_io<typename ARRAY_T::value_type>, HistoryMetadata >::type
  ArrayIndicesMetadata( string const & name, ARRAY_T const & arr, localIndex const num_indices )
  {
    size_t sizes[2] = { integer_conversion<size_t>(num_indices), integer_conversion<size_t>(arr.size(0)) };
    return HistoryMetadata(name, 2, &sizes[0], std::type_index(typeid(typename ARRAY_T::value_type)));
  }

  // class DataSpec
  // {
  // public:
  //   DataSpec( ):
  //     m_title( ),
  //     m_id( ),
  //     m_is_final(false),
  //     m_data_size(0),
  //     m_data_count(0),
  //     m_data_names(),
  //     m_data_offsets(),
  //     m_data_sizes(),
  //     m_type_byte_size(),
  //     m_data_types()
  //   {}

  //   template <typename T>
  //   void Append(localIndex const num_items,
  //               localIndex const units_per_item,
  //               string const & name)
  //   {
  //     localIndex data_size = num_items * units_per_item * sizeof(T);
  //     m_data_names.push_back( name );
  //     m_data_counts.push_back( num_items );
  //     m_data_offsets.push_back( m_data_size );
  //     m_type_byte_size.push_back( sizeof(T) );
  //     m_data_sizes.push_back( data_size );
  //     m_data_types.push_back( std::type_index(typeid(T)) );
  //     m_data_size += data_size;
  //     m_data_count++;
  //   }

  //   inline void Finalize() { m_is_final = true; }

  //   inline void SetTitleID( string const & title, string const & id )
  //   {
  //     m_title = title;
  //     m_id = id;
  //   }

  //   inline const string & getTitle() const { return m_title; }
  //   inline const string & getID() const { return m_id; }
  //   inline bool isFinal() const { return m_is_final; }

  //   inline size_t getTotalDataSize() const { warn_not_final(); return m_data_size;}
  //   inline localIndex getDiscreteDataCount() const { warn_not_final(); return m_data_count; }

  //   inline size_t getDataSubcount(localIndex idx) const { warn_not_final(); return m_data_counts[idx]; }
  //   inline size_t getDataSize(localIndex idx) const { warn_not_final(); return m_data_sizes[idx]; }
  //   inline size_t getTypeSize(localIndex idx) const { warn_not_final(); return m_type_byte_size[idx]; }
  //   inline std::type_index  getDataType(localIndex idx) const { warn_not_final(); return m_data_types[idx]; }
  //   inline size_t getDataOffset(localIndex idx) const { warn_not_final(); return m_data_offsets[idx]; }
  //   inline string const & getDataName(localIndex idx) const { warn_not_final(); return m_data_names[idx]; }

  //   inline size_t const * getDataSubcounts() const { warn_not_final(); return &m_data_counts[0]; }
  //   inline size_t const * getDataSizes() const { warn_not_final(); return &m_data_sizes[0]; }
  //   inline size_t const * getTypeSizes() const { warn_not_final(); return &m_type_byte_size[0]; }
  //   inline std::type_index  const * getDataTypes() const { warn_not_final(); return &m_data_types[0]; }
  //   inline size_t const * getDataOffsets() const { warn_not_final(); return &m_data_offsets[0]; }
  //   inline string const * getDataNames() const { warn_not_final(); return &m_data_names[0]; }

  // protected:
  //   void warn_not_final() const
  //  {
  //     if(!m_is_final)
  //     {
  //       GEOSX_WARNING("DataSpec: The data specification being operated on has not been finalized.");
  //     }
  //   }

  //   string m_title;
  //   string m_id;

  //   bool m_is_final;

  //   size_t m_data_size;
  //   localIndex m_data_count;

  //   std::vector<string> m_data_names;
  //   std::vector<size_t> m_data_counts;
  //   std::vector<size_t> m_data_offsets;
  //   std::vector<size_t> m_data_sizes;
  //   std::vector<size_t> m_type_byte_size;
  //   std::vector<std::type_index> m_data_types;
  // };

  // template < typename ARRAY_T >
  // inline
  // typename std::enable_if < is_array<ARRAY_T> && (ARRAY_T::ndim == 1) && can_history_io<typename ARRAY_T::value_type>, void >::type
  // AppendArraySpec( DataSpec & spec, ARRAY_T const & arr, string const & name, localIndex vals_per_col = 1 )
  // {
  //   spec.Append<typename ARRAY_T::value_type>(arr.size( ), vals_per_col, name);
  // }

  // template < typename ARRAY_T >
  // inline
  // typename std::enable_if < is_array<ARRAY_T> && (ARRAY_T::ndim > 1) && can_history_io<typename ARRAY_T::value_type>, void >::type
  // AppendArraySpec( DataSpec & spec, ARRAY_T const & arr, string const & name)
  // {
  //   spec.Append<typename ARRAY_T::value_type>(arr.size( ) / arr.size( 0 ), arr.size( 0 ), name);
  // }

  // template < typename ARRAY_T >
  // inline
  // typename std::enable_if < is_array<ARRAY_T> && (ARRAY_T::ndim == 1) && can_history_io<typename ARRAY_T::value_type>, void >::type
  // AppendArrayIndicesSpec( DataSpec & spec, ARRAY_T const & GEOSX_UNUSED_PARAM( arr ), string const & name, localIndex const num_indices, localIndex vals_per_col = 1 )
  // {
  //   spec.Append<typename ARRAY_T::value_type>(num_indices, vals_per_col, name);
  // }

  // template < typename ARRAY_T >
  // inline
  // typename std::enable_if < is_array<ARRAY_T> && (ARRAY_T::ndim > 1) && can_history_io<typename ARRAY_T::value_type>, void >::type
  // AppendArrayIndicesSpec( DataSpec & spec, ARRAY_T const & arr, string const & name, localIndex const num_indices )
  // {
  //   spec.Append<typename ARRAY_T::value_type>(num_indices, arr.size( 0 ), name);
  // }

}

#endif