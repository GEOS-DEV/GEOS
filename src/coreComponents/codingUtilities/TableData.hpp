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
template< typename T >
constexpr bool is_string = std::is_same_v< T, std::string >;
// Class for managing table m_data
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
  void addRow( std::vector< string > row );

  /**
   * @brief Reset data in the table
   */
  void clear();

  /**
   * @return The rows of the table
   */
  std::vector< std::vector< string > > & getTableDataRows();

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
  void addCell( real64 rowValue, real64 columnValue, T value );

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
    string cellValue = GEOS_FMT( "{}", args );
    static_assert( is_string< decltype(cellValue) >, "An argunment passed in addRow cannot be converted to string" );
    m_cellsValue.push_back( cellValue );
  } (), ...);

  m_rows.push_back( m_cellsValue );
}

template< typename T >
void TableData2D::addCell( real64 rowValue, real64 columnValue, T value )
{
  std::pair< real64, real64 > id = std::pair< real64, real64 >( rowValue, columnValue );
  m_data[id] = GEOS_FMT( "{}", value );
  m_columns.insert( columnValue );
  m_rows.insert( rowValue );
}

}

#endif
