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

class TableData
{
public:

  /**
   * @brief Add a row to the table.
   * @param Args The values passed to addRow (can be any type).
   * @param args Cell values to be added to the line.
   */
  template< typename ... Args >
  void addRow( Args const & ... args );

  std::vector< std::vector< string > > & getTableDataRows();

private:

  std::vector< std::vector< string > > m_rows;

};

class TableData2D : public TableData
{
public:

  /**
   * @brief Add a cell to the table.
   * Construct a map of pair<x,y> and cell value
   * @param T The value passed to addCell (can be any type).
   * @param value Cell value to be added.
   */
  template< typename T >
  void addCell( real64 x, real64 y, T value );

  /**
   * @brief Construct all rows from all cell values stored in map previously
   * @param tableRows Rows to be built
   */
  void buildRows( std::vector< std::vector< string > > & tableRows );

private:
  std::map< std::pair< real64, real64 >, string > data;
  std::set< real64 > columns;
  std::set< real64 > row;
};

template< typename ... Args >
void TableData::addRow( Args const &... args )
{
  //int idx = 0;
  std::vector< string > m_cellsValue;
  ( [&] {
    string cellValue = GEOS_FMT( "{}", args );
    // if( m_columns[idx].parameter.enabled )
    // {
    m_cellsValue.push_back( cellValue );
    // }
  } (), ...);

  m_rows.push_back( m_cellsValue );
}

template< typename T >
void TableData2D::addCell( real64 rowValue, real64 columnValue, T value )
{
  std::pair< real64, real64 > id = std::pair< real64, real64 >( rowValue, columnValue );
  data[id] = GEOS_FMT( "{}", value );
  columns.insert( columnValue );
  row.insert( rowValue );
}

}

#endif
