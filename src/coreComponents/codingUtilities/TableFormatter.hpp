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

// Class for formatting table data
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
  TableFormatter( TableLayout const & tableLayout );

  /**
   * @brief Fill the vector (m_column) in tableData with values from m_rows in tableLayout, storing all values in an unsorted order.
   * @param columns Vector of columns to be filled.
   * @param tableData Vector of table data.
   */
  void fillTableColumnsFromRows( std::vector< TableLayout::Column > & columns,
                                 std::vector< std::vector< string > > const & tableData ) const;

protected:

  TableLayout m_tableLayout;
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
  TableCSVFormatter( TableLayout const & tableLayout );

  /**
   * @brief Convert the table data to a CSV string.
   * @param tableData The 1D table data.
   * @return The CSV string representation of the table data.
   */
  string dataToString( TableData const & tableData ) const;

  /**
   * @return The string with all column names.
   */
  string headerToString() const;

  /**
   * @brief Convert the TableData to a table string.
   * @param tableData The TableData to convert.
   * @return The table string representation of the TableData.
   */
  string toString( TableData const & tableData ) const;

};

class TableTextFormatter : public TableFormatter
{

public:

  TableTextFormatter( TableLayout const & tableLayout );

  TableTextFormatter( std::vector< string > const & columnNames );

  /**
   * @brief Convert the TableData to a table string.
   * @param tableData The TableData to convert.
   * @return The table string representation of the TableData.
   */
  string toString( TableData const & tableData ) const;

  /**
   * @brief Converts a TableLayout into a formatted string representation.
   * @return string
   */
  string layoutToString() const;

private:

  /**
   * @brief Converts a TableLayout into a formatted representation.
   * @param tableOutput The output stream
   * @param columns The vectors of table columns
   * @param nbRows Number of rows in the table
   * @param sectionSeparator An empty string for building the section separator
   */
  void layoutToString( std::ostringstream & tableOutput,
                       std::vector< TableLayout::Column > & columns,
                       integer const nbRows,
                       string & sectionSeparator ) const;

  /**
   * @brief Split all header names by detecting the newline \\n character.
   * @param splitHeader A empty vector who will contain all split header names
   * @param largestHeaderVectorSize The largest split header vector size
   */
  void parseAndStoreHeaderSections( std::vector< TableLayout::Column > const & columns,
                                    size_t & largestHeaderVectorSize,
                                    std::vector< std::vector< string > > & splitHeader ) const;

  /**
   * @brief Set the same vector size for each split header and merge it into columns
   * @param columns The table columns to be merged
   * @param largestHeaderVectorSize The reference value for adjusting splitHeader vector
   * @param splitHeader The vector containing all split headers
   */
  void adjustHeaderSizesAndStore( std::vector< TableLayout::Column > & columns,
                                  size_t const & largestHeaderVectorSize,
                                  std::vector< std::vector< string > > & splitHeader ) const;

  /**
   * @brief For each column find and set the column's longest string
   */
  void findAndSetMaxStringSize( std::vector< TableLayout::Column > & columns, size_t const & nbRows ) const;

  /**
   * @brief Compute the largest string size in the table. If the table title is the largest string size in the table, recalculate for all
   * columns the \p m_maxStringSize value by adding extra characters
   * @param sectionlineLength The length of a section line
   * @param titleLineLength The length of a title line
   */
  void computeAndSetMaxStringSize( std::vector< TableLayout::Column > & columns,
                                   string::size_type const sectionlineLength,
                                   string::size_type const titleLineLength ) const;

  /**
   * @brief Compute and build the top and the section line separator
   * @param topSeparator An empty string to be built
   * @param sectionSeparator An empty string to be built
   */
  void computeAndBuildSeparator( std::vector< TableLayout::Column > & columns,
                                 string & topSeparator,
                                 string & sectionSeparator ) const;

  /**
   * @brief Build the table title section
   * @param titleRows Rows containing the title section.
   * @param topSeparator The top line separator
   * @param sectionSeparator The section line separator
   */
  void buildTitleRow( string & titleRows, string_view topSeparator, string_view sectionSeparator ) const;

  /**
   * @brief Build a section by specifying it's type ( header or section )
   * @param sectionSeparator Line separator between sections
   * @param rows A section row
   * @param nbRows Indicates the number of lines in a section
   * @param section The section to be built
   */
  void buildSectionRows( std::vector< TableLayout::Column > & columns,
                         string_view sectionSeparator,
                         std::ostringstream & rows,
                         integer const nbRows,
                         TableLayout::Section const section ) const;
};
}

#endif
