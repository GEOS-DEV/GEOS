/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HistoryDataSpec.hpp
 */

#ifndef GEOSX_HISTORY_DATA_SPEC_HPP_
#define GEOSX_HISTORY_DATA_SPEC_HPP_

#include "codingUtilities/traits.hpp"
#include "common/DataTypes.hpp"
#include "LvArray/src/Array.hpp"

namespace geosx
{
  /// A constexpr bool to determine wether a type is compatible with the history collected and IO operations.
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

  /**
   * @class HistoryMetadata
   * @brief A minimal class to specify information about time history information being collected and output.
   */
  class HistoryMetadata
  {
  public:
    /**
     * @brief Constructor for multi-dimensional array types.
     * @param name A name for the underlying information -- used to specify the data name in time history files.
     * @param rank The rank of the array of data being collected.
     * @param dims The extent of each dimension of the array being collected.
     * @param type The std::type_index of the array being collected (std::type_index(typeid(T)))
     */
    HistoryMetadata( const string & name, localIndex rank, size_t * dims, std::type_index type ) :
     m_name(name),
     m_rank(rank),
     m_dims(dims,dims+rank),
     m_type(type)
    {}
    /**
     * @brief Constructor for one-dimensional array types.
     * @param name A name for the underlying information -- used to specify the data name in time history files.
     * @param dims The extent of the one-dimensional array.
     * @param type The std::type_index of the array being collected (std::type_index(typeid(T)))
     */
    HistoryMetadata( const string & name, size_t count, std::type_index type ) :
      m_name(name),
      m_rank(1),
      m_dims(&count,&count+1),
      m_type(type)
    {}
    /**
     * @brief Set the name. Typically used for metadata collectors to avoid writing data with the same name
     *         to the history output files.
     * @param name The name to set.
     */
    void setName( const string & name )
    {
      m_name = name;
    }
    /**
     * @brief Get the name.
     * @return The name of the data being collected.
     */
    const string & getName( ) const
    {
      return m_name;
    }
    /**
     * @brief Set the type. Typically used for metadata collectors where local metadata is used to produce
     *         and output global metadata.
     * @param type The std::type_index of the type being collected. (std::type_index(typeid(T)))
     */
    void setType( std::type_index type )
    {
      m_type = type;
    }
    /**
     * @brief Get the type of the collected data.
     * @return The std::type_index of the collected data.
     */
    std::type_index getType( ) const
    {
      return m_type;
    }
    /**
     * @brief Get the total data count for the data being collected.
     * @return The number of data units of HistoryMetadata::getType( ) to be collected.
     */
    size_t getTypeCount( ) const
    {
      size_t local_size = 1;
      for( size_t dim : m_dims )
      {
        local_size *= dim;
      }
      return local_size;
    }
    /**
     * @brief Get the rank of the array data to be collected.
     * @return The rank.
     */
    localIndex getRank( ) const
    {
      return m_rank;
    }
    /**
     * @brief Get a pointer to the extent of each dimension.
     * @return The head of the array containing the dimensional extent of each dimension of the array data being collected.
     */
    const localIndex * getDims( ) const
    {
      return &m_dims[0];
    }
    /**
     * @brief Get the dimensional extent of the specified dimension.
     * @param dim The dimsion to get the extent off.
     * @return The extend of the dimension in the array data being collected.
     */
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

  /**
   * @brief Whether the type is a supported array type.
   * @tparam T The type to check.
   */
  template < typename T >
  constexpr bool is_array_type = traits::is_array_view< T > || traits::is_array< T >;

  /**
   * @brief Whether the type is a supported sorted array (set) type.
   * @tparam T The type to check.
   */
  template < typename T >
  constexpr bool is_sorted_array_type = traits::is_sorted_array_view< T > || traits::is_sorted_array< T >;

  /**
   * @brief Whether the type is a supported container for history collection and io operations.
   * @tparam T The type to check.
   */
  template < typename T >
  constexpr bool can_history_io_container = is_array_type< T > || is_sorted_array_type< T >;

  /**
   * @brief Produce a HistoryMetadata object for a supported one-dimensional array type.
   * @tparam ARRAY_T A supported array type.
   * @param name The name to give the metadata, usually dataRepository::Wrapper::getName by default.
   * @param arr The array to produce metadata about
   * @param size_override Specified in order to overwrite the actual size of the array with the size specified (used when collecting only a portion of the array data).
   * @return HistoryMetadata for the provided one-dimensional array.
   */
  template < typename ARRAY_T >
  inline
  HistoryMetadata getFlatArrayMetadata( string const & name, const ARRAY_T & arr, localIndex size_override = -1)
  {
    localIndex size = size_override < 0 ? arr.size( ) : size_override;
    size_t sz = LvArray::integerConversion<size_t>(size);
    return HistoryMetadata(name, sz, std::type_index(typeid(typename ARRAY_T::value_type)));
  }

  /**
   * @brief Produce a HistoryMetada object for one-dimensional LvArray::Array/ArrayView types.
   * @copydetails getFlatArrayMetadata( )
   */
  template < typename ARRAY_T >
  inline
  typename std::enable_if < is_array_type< ARRAY_T > && (ARRAY_T::ndim == 1) && can_history_io<typename ARRAY_T::value_type>, HistoryMetadata >::type
  getHistoryMetadata( string const & name, const ARRAY_T & arr, localIndex size_override = -1 )
  {
    return getFlatArrayMetadata< ARRAY_T >( name, arr, size_override );
  }

  /**
   * @brief Produce a HistoryMetadata object for LvArray::SortedArray/SortedArrayView types.
   * @copydetails getFlatArrayMetadata( )
   */
  template < typename ARRAY_T >
  inline
  typename std::enable_if < is_sorted_array_type< ARRAY_T > && can_history_io< typename ARRAY_T::value_type >, HistoryMetadata >::type
  getHistoryMetadata( string const & name, const ARRAY_T & arr, localIndex size_override = -1 )
  {
    return getFlatArrayMetadata< ARRAY_T >( name, arr, size_override );
  }

  /**
   * @brief Produce a HistoryMetadata object for multi-dimensional LvArray::Array/ArrayView types.
   * @copydetails getFlatArrayMetadata( )
   */
  template < typename ARRAY_T >
  inline
  typename std::enable_if < ( is_array_type< ARRAY_T >) && (ARRAY_T::ndim > 1) && can_history_io<typename ARRAY_T::value_type>, HistoryMetadata >::type
  getHistoryMetadata( string const & name, ARRAY_T const & arr, localIndex size_override = -1 )
  {
    size_t per_index_size = arr[ 0 ].size( );
    size_t num_indices = ( size_override >= 0 ? LvArray::integerConversion<size_t>( size_override ) : LvArray::integerConversion<size_t>( arr.size( ) / per_index_size ) );
    size_t sizes[2] = { num_indices, per_index_size };
    return HistoryMetadata(name, 2, &sizes[0], std::type_index(typeid(typename ARRAY_T::value_type)));
  }

  /**
   * @brief Produce a HistoryMetadata object for a fundamental type that can_history_io.
   * @tparam T The type to produce HistoryMetadata for.
   * @param name The name to give the metadata, usually dataRepository::Wrapper::getName by default.
   * @param type The data of type T to being used for history collection/output.
   * @param size_override Specified in order to overwrite the actual size of the data. Really only here to make the getHistoryMetadata overloaded function consistent,
   *                       but is still functional.
   */
  template < typename T >
  inline typename std::enable_if < can_history_io< T >, HistoryMetadata >::type
  getHistoryMetadata( string const & name, const T & GEOSX_UNUSED_PARAM( type ), localIndex size_override = -1 )
  {
    size_t size = size_override < 0 ? 0 : LvArray::integerConversion<size_t>(size_override);
    return HistoryMetadata(name, size, std::type_index(typeid(T)));
  }

  /**
   * @brief Fall-through implementation to catch attempts to collect history information about unsupported data type inside supported container types.
   * @tparam ARRAY_T A history collection/output supported collection type.
   * @param name Unused
   * @param type Unused
   * @param size_overrid Unused
   * @return A HistoryMetadata describing a size-zero array with name "NULL" and type_index(typeid(NULL)), will never actually return.
   */
  template < typename ARRAY_T >
  inline typename std::enable_if < ( can_history_io_container< ARRAY_T > ) && !can_history_io< typename ARRAY_T::value_type >, HistoryMetadata >::type
  getHistoryMetadata( string const & name, const ARRAY_T & type, localIndex size_override )
  {
    GEOSX_ERROR("Trying to use time history output on an array containing an unsupported type.");
    GEOSX_UNUSED_VAR( name );
    GEOSX_UNUSED_VAR( type );
    GEOSX_UNUSED_VAR( size_override );
    return HistoryMetadata("NULL", 0, std::type_index(typeid(NULL)));
  }

  /**
   * @brief Fall-through implementation to catch attempts to collect history that cannot be collected/output.
   * @tparam T A history collection/output unsupported type.
   * @param name Unused
   * @param type Unused
   * @param size_overrid Unused
   * @return A HistoryMetadata describing a size-zero array with name "NULL" and type_index(typeid(NULL)), will never actually return.
   */
  template < typename T >
  inline typename std::enable_if < ! ( can_history_io_container< T > ) && ! can_history_io< T >, HistoryMetadata >::type
  getHistoryMetadata( string const & name, const T & type, localIndex size_override )
  {
    GEOSX_ERROR("Trying to use time history output on an unsupported type.");
    GEOSX_UNUSED_VAR( name );
    GEOSX_UNUSED_VAR( type );
    GEOSX_UNUSED_VAR( size_override );
    return HistoryMetadata("NULL", 0, std::type_index(typeid(NULL)));
  }


}

#endif
