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
   * @brief Struct for a column margin and a border margin.
   */
  struct Margin
  {
    // 0 margin from |
    integer tiny;
    // 1 margin from |
    integer small;
    // 2 margin from |
    integer medium;
    // 3 margin from |
    integer large;
    // current margin value
    integer marginValue;

    void setWorkingValue( integer const value )
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
   * @brief Construct a new Table object by specifying value alignment and optionally their displays based to log levels
   * level
   * @param columnParameter
   */
  Table( std::vector< ColumnParam > const & columnParameter );

  /**
   * @brief Add a row the the table. The number of arguments must match with the number of header values
   * @tparam N
   * @tparam Args
   * @param args
   */
  template< size_t N, typename ... Args >
  void addRow( Args const &... args );

  /**
   * @brief Write the the table into specified stream
   * @param os
   */
  void draw( std::ostream & os = std::cout );

  /**
   * @brief Set the name of the table
   * @param tableTitle The name of the table
   */
  void setTitle( string_view tableTitle );

  /**
   * @brief Sets the minimal margin width between row content and borders.
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
   * @brief Get the margin type either a borderMagin or a column
   * @param type
   * @return Margin
   */
  Margin getMargin( MarginType type ) const;

  /**
   * @brief Get the name of the section given an enum
   * @param section
   * @return The name of the section
   */
  string getStringSection( Section section ) const;

  /**
   * @brief Fill the vector \p m_column with values from m_cellsRows who store all values in an unsorted order
   */
  void fillColumnsValuesFromCellsRows( );

  /**
   * @brief Split all header names by detecting the newline \\n character.
   * @param splitHeader A vector containing all split header names
   * @param largestHeaderVectorSize The largest split header vector size
   */
  void parseAndStoreHeaderSections( size_t & largestHeaderVectorSize,
                                    std::vector< std::vector< string > > & splitHeader );

  /**
   * @brief Iterate throught the header names vector
   * Adjust the size of each header vector by adding empty strings if needed.
   * Store the modified header names in the corresponding column parameter.
   * @param largestHeaderVectorSize The largest split header vector size
   * @param splitHeader A vector containing all split header names
   */
  void adjustHeaderSizesAndStore( size_t largestHeaderVectorSize,
                                  std::vector< std::vector< string > > & splitHeader );

  /**
   * @brief For each column find and set the column's longest string
   */
  void findAndSetMaxStringSize();

  /**
   * @brief Compute the largest string size in the table. If the table title is the largest string size in the table, recalculate for all
   * columns the \p m_maxStringSize value by adding extra characters
   * @param sectionlineLength The length of a section line
   * @param titleLineLength The length of a title line
   */
  void computeAndSetMaxStringSize( string::size_type sectionlineLength,
                                   string::size_type titleLineLength );

  /**
   * @brief Compute and build the top and the section line separator
   * @param topSeparator
   * @param sectionSeparator
   */
  void computeAndBuildSeparator( string & topSeparator, string & sectionSeparator );

  /**
   * @brief Build the table title section
   * @param titleRows Rows containing the title section
   * @param topSeparator The top line separator
   * @param sectionSeparator The section line separator
   */
  void buildTitleRow( string & titleRows, string topSeparator, string sectionSeparator );

  /**
   * @brief Build a section by specifying it's type ( header or section )
   * @param sectionSeparator Line separator between sections
   * @param rows A section row
   * @param nbRows Indicates the number of lines in a section
   * @param sectionName The section to be built
   */
  void buildSectionRows( string sectionSeparator,
                         string & rows,
                         integer const nbRows,
                         Section const sectionName );


  std::vector< std::vector< string > > m_cellsRows;
  std::vector< Column > m_columns;

  Margin borderMargin;
  Margin columnMargin;

  string tableTitle;

  int marginTitle = 2;

};

template< size_t N, typename ... Args >
void Table::addRow( Args const &... args )
{
  constexpr std::size_t nbColumn_ = sizeof...(args);
  static_assert( nbColumn_ == N,
                 "The number of cells per line does not correspond to the number of parameters." );

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
