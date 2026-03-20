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
 * @file TableData.hpp
 */

#ifndef GEOS_COMMON_FORMAT_TABLE_TABLEDATA_HPP
#define GEOS_COMMON_FORMAT_TABLE_TABLEDATA_HPP

#include "common/StdContainerWrappers.hpp"
#include "common/Units.hpp"
#include "common/DataTypes.hpp"
#include "common/format/Format.hpp"
#include "TableTypes.hpp"
#include <cstddef>

namespace geos
{

/**
 * @brief Class for managing table data
 */
class TableData
{
public:

  /// @cond DO_NOT_DOCUMENT
  TableData();

  TableData( TableData const & other );

  TableData( TableData && other );

  TableData & operator=( TableData const & other );

  TableData & operator=( TableData && other );
  ///@endcond

  /**
   * @brief Lexicographic sorting
   * @param other The table data to compate
   * @return true
   */
  bool operator<( TableData const & other ) const;

  /**
   * @brief Comparison operator for data rows
   * @param comparingTable The tableData values to compare
   * @return The comparison result
   */
  bool operator==( TableData const & comparingTable ) const;

  /**
   * @brief Representing a data in TableData
   */
  struct CellData
  {
    /// The cell type
    CellType type;
    /// The cell value
    string value;

    size_t getSerializedSize() const;

    void serialize( stdVector< buffer_unit_type > & out ) const;
  };

  /**
   * @tparam T The trivial type
   * @return Returns the size occupied by a trivial type in memory.
   */
  template< typename T >
  static unsigned long sizeOfField( T )
  { return sizeof(T); }

  /**
   * @brief Returns the size of a string (header size + content).
   * @param str The target string
   * @return Size in bytes.
   */
  static unsigned long sizeOfField( string_view str )
  { return sizeof(string::size_type) + str.size(); }

  template< typename T >
  static void serializePrimitive ( T const data, stdVector< buffer_unit_type > & out )
  {
    buffer_unit_type const * begin = reinterpret_cast< buffer_unit_type const * >( &data );
    buffer_unit_type const * end = begin + sizeof(data);
    out.insert( out.end(), begin, end );
  };

  static void serializeString ( string const & data, stdVector< buffer_unit_type > & out )
  {
    serializePrimitive( data.size(), out );
    auto * begin = data.data();
    auto * end = begin + data.size();
    out.insert( out.end(), begin, end );
  };

  template< typename T >
  static void deserializeField( T & data, buffer_unit_type const * & ptr, buffer_unit_type const * end )
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
  static void deserializeField( string & str, buffer_unit_type const * & ptr, buffer_unit_type const * end )
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


  /// Alias for table data rows with cells values
  using DataRows = stdVector< stdVector< CellData > >;

  void serialize( stdVector< buffer_unit_type > localTableData ) const;

  /**
   * @brief Add a row to the table.
   * The values passed to addRow (can be any type).
   * @param args CellData values to be added to the row.
   */
  template< typename ... Args >
  void addRow( Args const & ... args );

  /**
   * @brief Add a row to the table
   * @param row A vector of string representing a row
   */
  void addRow( stdVector< CellData > const & row );

  /**
   * @brief Add a line separator to the table
   * You must have filled values in TableData before using it
   */
  void addSeparator();

  /**
   * @brief Reset data in the table
   */
  void clear();

  /**
   * @brief Remove all errors
   */
  void clearErrors()
  { m_errors->clear(); }

  /**
   * @brief Get all error messages
   * @return The vector of error messages
   */
  stdVector< string > const & getErrorMsgs() const;

  /**
   * @return The const table data rows
   */
  DataRows const & getCellsData() const
  { return m_rows; }

  /**
   * @return The const table data rows
   */
  DataRows & getCellsData()
  { return m_rows; }

  /**
   * @brief Get all error messages
   * @return The list of error messages
   */
  TableErrorListing const & getErrorsList() const
  { return *m_errors; }

  /**
   * @brief Get all error messages
   * @return The list of error messages
   */

  TableErrorListing & getErrorsList()
  { return *m_errors; }

  /**
   * @brief Gather all the TableData rows to the rank 0
   * @param func The callable comparison function object to sort TableData rows, by default none
   */
  template< typename SortingFunc = std::nullptr_t >
  void gatherRowsRank0( SortingFunc && func );

private:
  /// @brief vector containing all rows with cell values
  DataRows m_rows;

  /// @brief Store all errors that can be found during the generation of the TableData
  std::unique_ptr< geos::TableErrorListing > m_errors;

};


/**
 * @brief Class for managing 2D table m_data
 */
class TableData2D
{
public:

  /// Type real64 for a row
  using RowType = real64;
  /// Type real64 for a column
  using ColumnType = real64;

  /// Struct containing conversion informations
  struct TableDataHolder
  {
    /// Vector containing all columns names
    /// A header value is presented as "pressure [K] = {}"
    stdVector< string > headerNames;
    /// TableData to be built
    TableData tableData;
  };

  /**
   * @brief Add a cell to the table. If necessary, create automatically the containing column & row.
   * @tparam T The value passed to addCell (can be any type).
   * @param rowValue The value of the row containing the cell.
   * @param columnValue The value of the column containing the cell.
   * @param value CellData value to be added.
   */
  template< typename T >
  void addCell( RowType rowValue, ColumnType columnValue, T const & value );

  /**
   * @brief Collects all the values needed to build the table
   * @param dim0AxisCoordinates Vector containing all row axis values
   * @param dim1AxisCoordinates Vector containing all column axis values
   * @param values Array containing all table values contiguously
   * @param columnMajorValues Set the row/column major convention
   */
  void collectTableValues( arrayView1d< real64 const > dim0AxisCoordinates,
                           arrayView1d< real64 const > dim1AxisCoordinates,
                           arrayView1d< real64 const > values,
                           bool columnMajorValues );

  /**
   * @brief Convert from 2D axis/values a structure the information needed to build a TableFormatter
   * @param coordX Array containing row axis values
   * @param coordY Array containing column axis values
   * @param rowAxisDescription The description for a row unit value
   * @param columnAxisDescription The description for a column unit value
   * @param values Vector containing all table values
   * @param columnMajorValues Set the row/column major convention
   * @param valueDescription The description of the value (typically, the value unit description)
   * @return A struct containing the tableData converted and all header values ;
   */
  TableData2D::TableDataHolder convertTable2D( arrayView1d< real64 const > coordX, arrayView1d< real64 const > coordY,
                                               string_view rowAxisDescription,
                                               string_view columnAxisDescription,
                                               arrayView1d< real64 const > const values,
                                               bool columnMajorValues,
                                               string_view valueDescription );

  /**
   * @return Convert and return a struct containing a 1D Table, the column names list from a TableData2D and any errors related to the table
   * @param dataDescription The table dataDescription shown at the top left side
   * @param rowFmt The y axis units of the table.
   * @param columnFmt  The x axis units of the table.
   * @note The rows and columns FMT can be customized. The bracket "{}" will be replaced by the axis value.
   * By default it displays the axis value.
   * I.E to display a customized axis to show the pressures in y axis, a rowFmt value can be : "pressure [K] = {}"
   */
  TableDataHolder buildTableData( string_view dataDescription,
                                  string_view rowFmt = "{}", string_view columnFmt = "{}" ) const;

  /**
   * @brief Clear all data stored in TableData
   */
  inline void clear()
  {
    m_data.clear();
    m_columnValues.clear();
    m_errors->clear();
  }

private:
  /// @brief all cell values by their [ row ][ column ]
  stdMap< RowType, stdMap< ColumnType, string > > m_data;

  /// @brief Store all column values when adding cell
  std::set< real64 > m_columnValues;
  /// @brief Store all errors that can be found during the generation of the TableData
  std::unique_ptr< geos::TableErrorListing > m_errors = std::make_unique< geos::TableErrorListing >();
};

/**
 * @brief Trait to check is the args is a special type of cell
 * @tparam T The type of a cell
 */
template< typename T >
constexpr bool isCellType = std::is_same_v< T, CellType >;

template< typename ... Args >
void TableData::addRow( Args const &... args )
{
  stdVector< CellData > cells;
  ( [&] {
    static_assert( has_formatter_v< decltype(args) > || isCellType< std::decay_t< decltype(args) > >, "Argument passed in addRow cannot be converted to string nor a CellType" );
    if constexpr (std::is_same_v< Args, CellType >) {
      cells.push_back( { args, string() } );
    }
    else if constexpr (std::is_floating_point_v< std::decay_t< decltype(args) > >) {
      if( !getErrorsList().hasErrors() && (std::isnan( args ) ||  std::isinf( args )))
      {
        m_errors->addError( "Warning : Invalid values detected (nan/inf)." );
      }
      cells.push_back( {CellType::Value, GEOS_FMT( "{}", args )} );
    }
    else
    {
      cells.push_back( {CellType::Value, GEOS_FMT( "{}", args )} );
    }
  } (), ...);
  addRow( cells );
}

/**
 * @brief Add a cell to the table.
 *
 * @tparam T The value type to insert.
 * @param rowValue The row key.
 * @param columnValue The column key.
 * @param value The value to store in the cell.
 */
template< typename T >
void TableData2D::addCell( real64 const rowValue, real64 const columnValue, T const & value )
{
  static_assert( has_formatter_v< decltype(value) >, "Argument passed in addCell cannot be converted to string" );
  m_columnValues.insert( columnValue );
  m_data.get_inserted( rowValue ).get_inserted( columnValue ) =  GEOS_FMT( "{}", value );
}

// Custom Comp function;
namespace tabledatasorting
{
/**
 * @brief Compare two string number string by  in ascending numerical order.
 * @param a The string to compare
 * @param b The string to compare
 * @return True if a is greater than b
 */
bool positiveNumberStringComp( string_view a, string_view b );
}

}
#endif /* GEOS_COMMON_FORMAT_TABLE_TABLEDATA_HPP */
