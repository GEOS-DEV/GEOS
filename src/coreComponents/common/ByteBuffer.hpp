/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ByteBuffer.hpp
 * @brief Append values to a byte buffer and read them back in order, for MPI messages.
 */

#ifndef GEOS_COMMON_BYTEBUFFER_HPP
#define GEOS_COMMON_BYTEBUFFER_HPP

#include <cstddef>
#include <cstring>
#include <stdexcept>
#include <type_traits>

namespace geos::bytes
{

/**
 * @brief Append raw bytes to a byte buffer.
 * @tparam BUFFER A contiguous container of a one-byte type (e.g. stdVector< char >).
 * @param buffer The buffer to extend.
 * @param data The bytes to append.
 * @param size The number of bytes to append.
 */
template< typename BUFFER >
void appendRaw( BUFFER & buffer, void const * const data, std::size_t const size )
{
  using Byte = typename BUFFER::value_type;
  static_assert( sizeof( Byte ) == 1, "A byte buffer stores one-byte values" );
  if( size == 0 )
  {
    return;
  }
  std::size_t const position = buffer.size();
  buffer.resize( position + size );
  std::memcpy( buffer.data() + position, data, size );
}

/**
 * @brief Append the bytes of a trivially copyable value to a byte buffer.
 * @tparam T The value type.
 * @tparam BUFFER A contiguous container of a one-byte type.
 * @param buffer The buffer to extend.
 * @param value The value to append.
 */
template< typename T, typename BUFFER >
void append( BUFFER & buffer, T const & value )
{
  static_assert( std::is_trivially_copyable_v< T >, "Only trivially copyable values are appended as bytes" );
  appendRaw( buffer, &value, sizeof( T ) );
}

/**
 * @brief Read values from a byte buffer in the order they were appended.
 * @details Every read checks the remaining size and throws std::invalid_argument
 * on a truncated buffer. The buffer must outlive the reader.
 */
class Reader
{
public:
  /**
   * @brief Read from a contiguous container of a one-byte type.
   * @tparam BUFFER The container type.
   * @param buffer The buffer to read.
   */
  template< typename BUFFER >
  explicit Reader( BUFFER const & buffer ):
    Reader( buffer.data(), buffer.size() )
  {
    static_assert( sizeof( typename BUFFER::value_type ) == 1, "A byte buffer stores one-byte values" );
  }

  /**
   * @brief Read from raw memory.
   * @param data The first byte.
   * @param size The number of bytes.
   */
  Reader( void const * const data, std::size_t const size ):
    m_data( static_cast< unsigned char const * >( data ) ),
    m_size( size )
  {}

  /// @return Whether every byte has been read.
  bool done() const { return m_cursor >= m_size; }

  /// @return The number of bytes not read yet.
  std::size_t remaining() const { return m_size - m_cursor; }

  /**
   * @brief Read the next value.
   * @tparam T A trivially copyable type.
   * @return The value.
   */
  template< typename T >
  T read()
  {
    static_assert( std::is_trivially_copyable_v< T >, "Only trivially copyable values are read as bytes" );
    T value;
    readRaw( &value, sizeof( T ) );
    return value;
  }

  /**
   * @brief Copy the next bytes.
   * @param destination Where to copy the bytes.
   * @param size The number of bytes.
   */
  void readRaw( void * const destination, std::size_t const size )
  {
    char const * const source = view( size );
    if( size > 0 )
    {
      std::memcpy( destination, source, size );
    }
  }

  /**
   * @brief Skip the next bytes and return where they start, without copying.
   * @param size The number of bytes.
   * @return A pointer to the first of these bytes.
   */
  char const * view( std::size_t const size )
  {
    if( size > remaining() )
    {
      throw std::invalid_argument( "Truncated byte buffer" );
    }
    char const * const result = reinterpret_cast< char const * >( m_data + m_cursor );
    m_cursor += size;
    return result;
  }

private:
  unsigned char const * m_data;
  std::size_t m_size;
  std::size_t m_cursor = 0;
};

} // namespace geos::bytes

#endif // GEOS_COMMON_BYTEBUFFER_HPP
