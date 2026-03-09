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
 * @file LogHistory.hpp
 */

#ifndef GEOS_COMMON_LOGGER_MSG_REPORT_DATA_HPP
#define GEOS_COMMON_LOGGER_MSG_REPORT_DATA_HPP

#include "common/DataTypes.hpp"
#include "common/StdContainerWrappers.hpp"
#include "common/format/LogPart.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "DiagnosticMessage.hpp"
#include <string>


namespace geos
{

/**
 * @brief Keep track of all diagnostic message occured during the simulation
 */
class LogHistory
{
public:

  struct LogRecord
  {
    /** @brief Identifier for a diagnostic message (source localization). */
    struct Key
    {
      string m_filename; ///< Source file name
      integer m_lineId;  ///< Line number in the file.
    } m_key;

    /** @brief Content and metadata of the diagnostic message. */
    struct Values
    {
      string m_logPart;  ///< The string logPart.
      MsgType m_msgType; ///< Message type.
      integer m_count;   ///< Number of occurrences detected.
    } m_value;

    /** @brief Calculates the total size required for the serialization.
     * @return Size in bytes.
     */
    size_t getSerializedSize() const;

    /**
     * @brief Construct an empty Log Record object
     */
    LogRecord();

    /**
     * @brief Construct a new Log Record object
     * @param key The log record key
     * @param values The log record values
     */
    LogRecord( Key const & key, Values const & values );

    /** @brief Serializes the record fields into a binary buffer.
     * @param out Destination vector for the serialized data.
     */
    void serialize( stdVector< buffer_unit_type > & out ) const;

    /**
     * @brief Deserializes a complete record and advances the read pointer.
     * @param logRecordBytes Reference to the read pointer.
     * @param end Upper limit of readable memory.
     */
    void deserialize( buffer_unit_type const * & logRecordBytes, buffer_unit_type const * end );

    /**
     * @tparam T The trivial type
     * @return Returns the size occupied by a trivial type in memory.
     */
    template< typename T >
    unsigned long sizeOfField( T ) const
    { return sizeof(T); }

    /**
     * @brief Returns the size of a string (header size + content).
     * @param str The target string
     * @return Size in bytes.
     */
    unsigned long sizeOfField( string_view str ) const
    { return sizeof(string::size_type) + str.size(); }

    /** @brief Reads a trivial value from the buffer and advances the pointer.
     * @param data Destination variable.
     * @param ptr Current read pointer (advanced by sizeof(T)).
     * @param end Safety: maximum buffer limit.
     */
    template< typename T >
    void deserializeField( T & data, buffer_unit_type const * & ptr, buffer_unit_type const * end )
    {
      static_assert( std::is_trivially_copyable_v< T > );
      if( ptr + sizeof(T)> end ) throw std::runtime_error( "Buffer truncated" );
      memcpy( &data, ptr, sizeof(T) );
      ptr += sizeof(T);
    }

    /** @brief Reads a string value from the buffer and advances the pointer.
     * @param data Destination variable.
     * @param ptr Current read pointer (advanced by sizeof(string)).
     * @param end Safety: maximum buffer limit.
     */
    void deserializeField( string & str, buffer_unit_type const * & ptr, buffer_unit_type const * end )
    {
      string::size_type strSize = 0;
      deserializeField( strSize, ptr, end );
      if( std::distance( ptr, end ) < (long) strSize )
      {
        throw std::runtime_error( "Buffer truncated reading string" );
      }
      str.assign( ptr, ptr + strSize );
      ptr += str.size();
    }

  };

  /** * @brief Records a diagnostic message occurrence in the history.
   * @param logPartName The string log part name.
   * @param diagMsg The diagnostic message associated
   */
  void recordDiagnostic( string_view logPartName, DiagnosticMsg const & diagMsg );

  /**
   * @brief Display the diagnostic statistics to the log
   */
  void diagnosticStatsReport();

  /**
   * @return The const historical diagnostic
   */
  auto const & getDiagnosticHistory() const
  { return m_diagnosticHistory; }

private:

  /// @cond DO_NOT_DOCUMENT
  struct LocationKeyHash
  {

    size_t operator()( LogRecord::Key const & key ) const noexcept
    {
      size_t h1 = std::hash< string >{} (key.m_filename);
      size_t h2 = std::hash< integer >{} (key.m_lineId);

      return h1 ^ (h2 << 1);
    }

  };
  /// @endcond

  /// @cond DO_NOT_DOCUMENT
  struct LocationKeyEqual
  {

    bool operator()( LogRecord::Key const & lhs,
                     LogRecord::Key const & rhs ) const
    {
      return lhs.m_filename == rhs.m_filename  &&
             lhs.m_lineId == rhs.m_lineId;
    }
  };
  /// @endcond

  /**
   * @brief Insert a LogRepord in the m_diagnosticHistory
   * @param log The logRecord with all the information
   */
  void insertDiagnosticReport( LogRecord log );

  /**
   * @brief Diagnostic history happened during the simulation
   */
  stdUnorderedMap< LogRecord::Key, LogRecord::Values,
                   LocationKeyHash, LocationKeyEqual > m_diagnosticHistory;
};

/**
 * @brief Template specialisation to convert a LogHistory to a table string.
 * @param logHistory The LogHistory object to convert.
 * @return The CSV string representation of the logHistory.
 */
template<>
string TableTextFormatter::toString< LogHistory >( LogHistory const & logHistory ) const;

}

#endif
