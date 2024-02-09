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
 * @file Table.hpp
 */

#ifndef GEOS_COMMON_TABLE_HPP
#define GEOS_COMMON_TABLE_HPP

#include "common/DataTypes.hpp"

namespace geos
{

class Table
{
public:

  enum Alignment { right, left, middle };

  enum MarginColumn {small, medium, large};

  /**
   * @brief Struct for a ColumnParam.
   * @param [in] headerName A vector containing all string for a header name
   * @param [in] columnValues Alignment for a column. By default aligned to the right side
   * @param [in] enabled A boolean to display a colummn
   */
  struct ColumnParam
  {
    std::vector< std::string > headerName;
    Alignment alignment = Alignment::right;
    bool enabled = true;
  };

  /**
   * @brief Construct a new Table object, all values in the table are centered by default
   *
   * @param columnNames
   */
  Table( std::vector< std::string > const & columnNames );

  /**
   * @brief Construct a new Table object by specifying alignment of values in the table and optionally their display according to their log
   *level
   *
   * @param columnParameter
   */
  Table( std::vector< ColumnParam > const & columnParameter );

  /**
   * @brief Add a row the the table. The number of arguments need to match with the number of header values
   *
   * @tparam N
   * @tparam Args
   * @param args
   */
  template< size_t N, typename ... Args >
  void addRow( const Args &... args );

  /**
   * @brief Display the table
   *
   * @param os
   */
  void draw( std::ostringstream & os );

  /**
   * @brief Set the name of the table
   * @param [in] title The name of the table
   */
  void setTitle( const std::string & title );

  /**
   * @brief Get the table title
   *
   * @return return the table name
   */
  std::string const & getTitle();

private:

  enum Section { header, values };

  /**
   * @brief Struct for a column.
   * @param [in] parameter Column parameter
   * @param [in] columnValues All column values
   * @param [in] m_maxStringSize The largest string in the column
   */
  struct Column
  {
    ColumnParam parameter;
    std::vector< std::string > columnValues;
    std::string m_maxStringSize;
  };

  /**
   * @brief Get the name of the section given an enum
   *
   * @param section
   * @return The name of the section
   */
  std::string getStringSection( Section section ) const;

  /**
   * @brief the vector column \p m_column with all values previously stored in the constructor
   *
   */
  void fillColumnsValuesFromCellsRows();

  /**
   * @brief Split all header names by detecting \p '\n' character and store each split header name into a vector.
   *
   */
  void splitHeadersStringAndStore();

  /**
   * @brief Iterate throught all header name vector and set them to the same size by adding an empty string.
   *
   */
  void addSpaceToSplitHeaderAndStore();

  /**
   * @brief For each column find the largest string size
   *
   */
  void findMaxStringSize();

  /**
   * @brief If the title is the largest string size in the table, recalculate for all column the \p m_maxStringSize value
   * by adding an extra number of characters
   *
   */
  void computeAndSetMaxStringSize( string::size_type const & sectionlineLength,
                                   string::size_type const & titleLineLength );
  /**
   * @brief Compute the lines separator
   *
   */
  void computeAndBuildLines();

  /**
   * @brief Build the title section
   *
   */
  void buildTitleRow();

  /**
   * @brief build the header or values section
   *
   * @param nbRows
   * @param sectionName
   */
  void buildSectionRows( integer const & nbRows, Section const & sectionName );

  std::vector< std::vector< std::string > > m_splitHeader;
  std::vector< std::vector< string > > m_cellsRows;
  std::vector< Column > m_columns;

  std::string title;
  std::string topSeparator;
  std::string sectionSeparator;
  size_t maxRowHeader;

  int marginTitle = 2;
  int borderMargin = 2;
  int columnMargin = 5;

  std::string titleRow;
  std::string rows;

};

template< size_t N, typename ... Args >
void Table::addRow( const Args &... args )
{
  constexpr std::size_t nbColumn_ = sizeof...(args);
  static_assert( nbColumn_ == N,
                 "The number of cells per line does not correspond to the number of parameters" );

  std::vector< std::string > rowsValues;
  int idx = 0;
  ( [&] {
    std::string cellValue = GEOS_FMT( "{}", args );
    if( m_columns[idx].parameter.enabled )
    {
      rowsValues.push_back( cellValue );
    }
  } (), ...);

  m_cellsRows.push_back( rowsValues );
}

}

#endif
