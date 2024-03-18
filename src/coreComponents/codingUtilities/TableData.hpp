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

namespace geos
{

#if __cplusplus < 202002L
template< class T >
static constexpr bool has_formatter = fmt::has_formatter< fmt::remove_cvref_t< T >, fmt::format_context >();
#else
template< typename T >
concept has_formatter = requires ( T& v, std::format_context ctx )
{
  std::formatter< std::remove_cvref_t< T > >().format( v, ctx );
};
#endif

// Class for managing table data
class TableData
{
public:

  /**
   * @brief Add a row to the table.
   * @param Args The values passed to addRow (can be any type).
   * @param args Cell values to be added to the row.
   */
  template< typename ... Args >
  void addRow( Args const & ... args );

  /**
   * @brief Add a row to the table
   * @param row A vector of string who contains cell Values
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

private:

  std::vector< std::vector< string > > m_rows;

};

// Class for managing 2D table m_data
class TableData2D
{
public:

  /**
   * @brief Add a cell to the table.
   * Construct a map of pair<x,y> and cell value
   * @param T The value passed to addCell (can be any type).
   * @param value Cell value to be added.
   */
  template< typename T >
  void addCell( real64 rowValue, real64 columnValue, T const & value );

  /**
   * @brief Construct a TableData from a Table2D
   * @return A TableData
   */
  TableData buildTableData() const;

  /**
   * @return return all columns values for 2D table
   */
  std::set< real64 > const & getColumns() const;

  /**
   * @return return all rows values for 2D table
   */
  std::set< real64 > const & getRows() const;

private:
  std::map< std::pair< real64, real64 >, string > m_data;
  std::set< real64 > m_columns;
  std::set< real64 > m_rows;
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

  m_rows.push_back( m_cellsValue );
}

template< typename T >
void TableData2D::addCell( real64 const rowValue, real64 const columnValue, T const & value )
{
  static_assert( has_formatter< decltype(value) >, "Argument passed in addCell cannot be converted to string" );
  std::pair< real64, real64 > const id = std::pair< real64, real64 >( rowValue, columnValue );
  m_data[id] = GEOS_FMT( "{}", value );
  m_columns.insert( columnValue );
  m_rows.insert( rowValue );
}

}

#endif
