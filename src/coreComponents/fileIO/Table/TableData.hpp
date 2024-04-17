/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TableData.hpp
 */

#ifndef GEOS_COMMON_TableData_HPP
#define GEOS_COMMON_TableData_HPP

#include "common/DataTypes.hpp"
#include "common/Format.hpp"

namespace geos
{

/**
 * @brief Class for managing table data
 */
class TableData
{
public:

  /**
   * @brief Add a row to the table.
   * The values passed to addRow (can be any type).
   * @param args Cell values to be added to the row.
   */
  template< typename ... Args >
  void addRow( Args const & ... args );

  /**
   * @brief Add a row to the table
   * @param row A vector of string representing a row
   */
  void addRow( std::vector< string > const & row );

  /**
   * @brief Reset data in the table
   */
  void clear();

  /**
   * @return The rows of the table
   */
  std::vector< std::vector< string > > const & getTableDataRows() const;

  /**
   * @brief Get all error messages
   * @return The vector of error messages
   */
  std::vector< string > const & getErrorMsgs() const;

  /**
   * @brief Set an error message
   * @param msg The error msg to vector
   */
  void addErrorMsgs( string const & msg );

protected:
  /// vector containing all rows with cell values
  std::vector< std::vector< string > > m_rows;

  /// store error if there are any inconsistencies related to the table
  std::vector< string > errorsMsg;

};

/**
 * @brief Class for managing 2D table m_data
 */
class TableData2D
{
public:

  /// Struct containing conversion informations
  struct Conversion1D
  {
    /// Vector containing all columns names
    std::vector< string > headerNames;
    /// TableData to be built
    TableData tableData;
  };

  /**
   * @brief Add a cell to the table. If necessary, create automatically the containing column & row.
   * @tparam T The value passed to addCell (can be any type).
   * @param value Cell value to be added.
   * @param rowValue The value of the row containing the cell.
   * @param columnValue The value of the column containing the cell.
   */
  template< typename T >
  void addCell( real64 rowValue, real64 columnValue, T const & value );

  /**
   * @return Convert and return a struct containing a 1D Table, the column names list from a TableData2D and any errors related to the table
   * @param dataDescription The table dataDescription shown at the top left side
   * @param rowFmt The y axis units of the table.
   * @param columnFmt  The x axis units of the table.
   * @note The rows and columns FMT can be customized. The bracket "{}" will be replaced by the axis value.
   * By default it displays the axis value.
   * I.E to display a customized axis to show the pressures in y axis, a rowFmt value can be : "pressure [K] = {}"
   */
  Conversion1D buildTableData( string_view dataDescription, string_view rowFmt = "{}", string_view columnFmt = "{}" ) const;

private:
  using RowType = real64;
  using ColumnType = real64;

  /// @brief all cell values by their [ row ][ column ]
  std::map< RowType, std::map< ColumnType, string > > m_data;
};

template< typename ... Args >
void TableData::addRow( Args const &... args )
{
  std::vector< string > m_cellsValue;
  ( [&] {
    static_assert( has_formatter< decltype(args) >, "Argument passed in addRow cannot be converted to string" );
    string const cellValue = GEOS_FMT( "{}", args );
    m_cellsValue.push_back( cellValue );
  } (), ...);

  addRow( m_cellsValue );
}

template< typename T >
void TableData2D::addCell( real64 const rowValue, real64 const columnValue, T const & value )
{
  static_assert( has_formatter< decltype(value) >, "Argument passed in addCell cannot be converted to string" );
  m_data[rowValue][columnValue] = GEOS_FMT( "{}", value );
}


}

#endif
