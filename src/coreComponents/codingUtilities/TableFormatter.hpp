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
 * @file TableFormatter.hpp
 */

#ifndef GEOS_COMMON_TABLEFORMATTER_HPP
#define GEOS_COMMON_TABLEFORMATTER_HPP

#include "codingUtilities/TableData.hpp"
#include "codingUtilities/TableLayout.hpp"

namespace geos
{

class TableFormatter
{
public:

  /**
   * @brief Constructor by default
   */
  TableFormatter() = default;

  /**
   * @brief Construct a new Table Formatter from a tableLayout
   * @param tableLayout Contain all column names and optionnaly the table title
   */
  TableFormatter( TableLayout tableLayout );

  TableLayout m_tableLayout;

  /**
   * @brief Fill the vector (m_column) in tableData with values from m_rows in tableLayout who store all values in an unsorted order
   */
  void fillTableColumnsFromRows( std::vector< TableLayout::Column > & columns,
                                 std::vector< std::vector< string > > const & tableData );

};

class TableCSVFormatter : public TableFormatter
{
public:
  /**
   * @brief Constructor by default
   */
  TableCSVFormatter() = default;

  /**
   * @brief Construct a new Table Formatter from a tableLayout
   * @param tableLayout Contain all column names and optionnaly the table title
   */
  TableCSVFormatter( TableLayout tableLayout );

  /**
   * @param tableData A 2-dimensions tabke
   * @return A string of CSV data from a 2-dimensions table 
   */
  string dataToString( TableData2D tableData );

  /**
   * @param tableData A 2-dimensions tabke
   * @return A string of CSV data from a 1-dimensions table 
   */
  string dataToString( TableData tableData );

  /**
   * @param columns 
   * @param nbRows 
   * @return string 
   */
  string headerToString( std::vector< TableLayout::Column > & columns, integer nbRows );

};

class TableTextFormatter : public TableFormatter
{

public:

  TableTextFormatter( TableLayout tableLayout );

  /**
   * @brief return a string following the formatter
   */
  string ToString( TableData & tableData );

  /**
   * @brief return a string following the formatter
   */
  string ToString( TableData2D & tableData );

private:

  /**
   * @brief Split all header names by detecting the newline \\n character.
   * @param splitHeader A empty vector who will contain all split header names
   * @param largestHeaderVectorSize The largest split header vector size
   */
  void parseAndStoreHeaderSections( std::vector< TableLayout::Column > & columns,
                                    size_t & largestHeaderVectorSize,
                                    std::vector< std::vector< string > > & splitHeader );

  string constructTable( std::vector< std::vector< string > > & rowsValues );
  /**
   * @brief Iterate throught the header names vector.
   * Adjust the size of each header vector by adding empty strings if needed.
   * Store the modified header names in the corresponding column parameter.
   * @param largestHeaderVectorSize The largest split header vector size
   * @param splitHeader A vector containing all split header names
   */
  void adjustHeaderSizesAndStore( std::vector< TableLayout::Column > & columns,
                                  size_t largestHeaderVectorSize,
                                  std::vector< std::vector< string > > & splitHeader );

  /**
   * @brief For each column find and set the column's longest string
   */
  void findAndSetMaxStringSize( std::vector< TableLayout::Column > & columns, size_t nbRows );

  /**
   * @brief Compute the largest string size in the table. If the table title is the largest string size in the table, recalculate for all
   * columns the \p m_maxStringSize value by adding extra characters
   * @param sectionlineLength The length of a section line
   * @param titleLineLength The length of a title line
   */
  void computeAndSetMaxStringSize( std::vector< TableLayout::Column > & columns,
                                   string::size_type sectionlineLength,
                                   string::size_type titleLineLength );

  /**
   * @brief Compute and build the top and the section line separator
   * @param topSeparator An empty string to be built
   * @param sectionSeparator An empty string to be built
   */
  void computeAndBuildSeparator( std::vector< TableLayout::Column > & columns,
                                 string & topSeparator,
                                 string & sectionSeparator );

  /**
   * @brief Build the table title section
   * @param titleRows Rows containing the title section.
   * @param topSeparator The top line separator
   * @param sectionSeparator The section line separator
   */
  void buildTitleRow( string & titleRows, string_view topSeparator, string_view sectionSeparator );

  /**
   * @brief Build a section by specifying it's type ( header or section )
   * @param sectionSeparator Line separator between sections
   * @param rows A section row
   * @param nbRows Indicates the number of lines in a section
   * @param section The section to be built
   */
  void buildSectionRows( std::vector< TableLayout::Column > & columns,
                         string_view sectionSeparator,
                         string & rows,
                         integer const nbRows,
                         TableLayout::Section const section );

};
}

#endif
