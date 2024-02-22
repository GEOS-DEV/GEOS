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

  enum MarginType { border, column };

  enum MarginValue { tiny, small, medium, large };

  /**
   * @brief Struct for a column margin and border margin.
   * @param tiny 0 margin for both;
   * @param small 1 margin from |
   * @param medium 2 margin from |
   * @param large 3 margin from |
   *
   * @param marginValue current margin value
   */
  struct Margin
  {
    integer tiny;
    integer small;
    integer medium;
    integer large;

    integer marginValue;

    void setWorkingValue( integer const & value )
    {
      marginValue = value;
    }
  };

  /**
   * @brief Struct for a ColumnParam.
   */
  struct ColumnParam
  {
    // A vector containing all string for a header name
    std::vector< string > headerName;
    // Alignment for a column. By default aligned to the right side
    Alignment alignment = Alignment::right;
    // A boolean to display a colummn
    bool enabled = true;
  };

  /**
   * @brief Construct a new Table object, all values in the table are centered by default
   * @param columnNames
   */
  Table( std::vector< string > const & columnNames );

  /**
   * @brief Construct a new Table object by specifying alignment of values in the table and optionally their display according to their log
   * level
   * @param columnParameter
   */
  Table( std::vector< ColumnParam > const & columnParameter );

  /**
   * @brief Add a row the the table. The number of arguments need to match with the number of header values
   * @tparam N
   * @tparam Args
   * @param args
   */
  template< size_t N, typename ... Args >
  void addRow( const Args &... args );

  /**
   * @brief Write the the table into specified stream
   * @param os
   */
  void draw( std::ostream & os = std::cout );

  /**
   * @brief Set the name of the table
   * @param [in] title The name of the table
   */
  void setTitle( string_view title );

  /**
   * @brief Sets the minimal margin width between the row content an its borders.
   * @param valueType
   */
  void setMargin( MarginValue valueType );

  /**
   * @return return the table name
   */
  string_view getTitle();

private:

  enum Section { header, values };

  /**
   * @brief Struct for a column.
   */
  struct Column
  {
    ColumnParam parameter;
    // A vector containing all column values
    std::vector< string > columnValues;
    // The largest string in the column
    string m_maxStringSize;
  };

  /**
   * @brief Get the Margin object
   * @param type
   * @return Margin
   */
  Margin getMargin( MarginType const & type );

  /**
   * @brief Get the name of the section given an enum
   * @param section
   * @return The name of the section
   */
  string getStringSection( Section section ) const;

  /**
   * @brief fill the vector \p m_column with values from m_cellsRows who store all value in an unsorted order
   */
  void fillColumnsValuesFromMCellsRows();

  /**
   * @brief Split all header names by detecting \p '\n' character and store each split header name into a vector.
   */
  void splitHeadersStringAndStore();

  /**
   * @brief Iterate throught all header name vector and set them to the same size by adding empty string(s).
   */
  void addSpaceToSplitHeaderAndStore();

  /**
   * @brief For each column find the largest string size
   */
  void findMaxStringSize();

  /**
   * @brief If the title is the largest string size in the table, recalculate for all column, the \p m_maxStringSize value
   * by adding an extra number of characters
   */
  void computeAndSetMaxStringSize( string::size_type sectionlineLength,
                                   string::size_type titleLineLength );
  /**
   * @brief Compute the lines separator
   */
  void computeAndBuildLines();

  /**
   * @brief Build the title section
   *
   */
  void buildTitleRow();

  /**
   * @brief build the header or values section
   * @param nbRows
   * @param sectionName
   */
  void buildSectionRows( integer const & nbRows, Section const & sectionName );

  std::vector< std::vector< string > > m_splitHeader;
  std::vector< std::vector< string > > m_cellsRows;

  std::vector< Column > m_columns;

  Margin borderMargin;
  Margin columnMargin;

  string title;
  string topSeparator;
  string sectionSeparator;
  size_t maxRowHeader;

  int marginTitle = 2;

  string titleRow;
  string rows;

};

template< size_t N, typename ... Args >
void Table::addRow( const Args &... args )
{
  constexpr std::size_t nbColumn_ = sizeof...(args);
  static_assert( nbColumn_ == N,
                 "The number of cells per line does not correspond to the number of parameters" );

  std::vector< string > rowsValues;
  int idx = 0;
  ( [&] {
    string cellValue = GEOS_FMT( "{}", args );
    if( m_columns[idx].parameter.enabled )
    {
      rowsValues.push_back( cellValue );
    }
  } (), ...);

  m_cellsRows.push_back( rowsValues );
}

}

#endif
