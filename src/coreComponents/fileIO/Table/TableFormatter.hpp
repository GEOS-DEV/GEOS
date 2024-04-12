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

#include "TableData.hpp"
#include "TableLayout.hpp"

namespace geos
{

/**
 * @brief abstract class for formatting table data
 */
class TableFormatter
{

protected:

  /// Layout for a table
  TableLayout m_tableLayout;

  /**
   * @brief Construct a new Table Formatter from a tableLayout
   * @param tableLayout Contain all column names and optionnaly the table title
   */
  TableFormatter( TableLayout const & tableLayout );

  /**
   * @brief Destroy the Table Formatter object
   */
  virtual ~TableFormatter() = default;

  /**
   * @brief Fill the vector (m_column) in tableData with values from rows stored in tableLayout.
   * @param columns Vector of columns to be filled.
   * @param tableData Vector containing all rows filled with values
   */
  void fillTableColumnsFromRows( std::vector< TableLayout::Column > & columns,
                                 std::vector< std::vector< string > > const & tableData,
                                 std::vector< string > & msgTableError ) const;

};

/**
 * @brief class for CSV formatting
 */
class TableCSVFormatter : public TableFormatter
{
public:

  /**
   * @brief Construct a new Table Formatter from a tableLayout
   * @param tableLayout Contain all column names and optionnaly the table title
   */
  TableCSVFormatter( TableLayout const & tableLayout );

  /**
   * @brief Destroy the TableCSVFormatter object
   */
  virtual ~TableCSVFormatter() = default;

  /**
   * @return The string with all column names.
   */
  string headerToString() const;

  /**
   * @brief Convert the table data to a CSV string.
   * @param tableData The 1D table data.
   * @return The CSV string representation of the table data.
   */
  string dataToString( TableData const & tableData ) const;

  /**
   * @brief Convert the TableData to a table string.
   * @param tableData The TableData to convert.
   * @return The table string representation of the TableData.
   */
  string toString( TableData const & tableData ) const;

};

/**
 * @brief class for log formatting
 */
class TableTextFormatter : public TableFormatter
{

public:

  /**
   * @brief Construct a new TableFormatter from a tableLayout
   * @param tableLayout Contain all column names and optionnaly the table title
   */
  TableTextFormatter( TableLayout const & tableLayout );

  /**
   * @brief Destroy the Table Text Formatter object
   */
  virtual ~TableTextFormatter() = default;

  /**
   * @brief Convert the TableData to a table string.
   * @param tableData The TableData to convert.
   * @return The table string representation of the TableData.
   */
  string toString( TableData const & tableData ) const;

private:

  /**
   * @brief Converts a TableLayout into a formatted representation.
   * @param tableOutput The output stream
   * @param columns The vector containing all table columns
   * @param msgTableError A vector containg all error related to the table
   * @param sectionSeparator An empty string for building the section separator
   */
  void layoutToString( std::ostringstream & tableOutput,
                       std::vector< TableLayout::Column > & columns,
                       std::vector< string > & msgTableError,
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
   * @param columns The vector containg all columns
   */
  void findAndSetMaxStringSize( std::vector< TableLayout::Column > & columns ) const;

  /**
   * @brief recalculate the largest string size for each columns
   * @param extraLines Extra characters to be added to \p m_maxStringSize of each columns
   */
  void recalculateMaxStringSize( std::vector< TableLayout::Column > & columns, integer const extraLines ) const;

  /**
   * @brief Compute the max table line length
   * @param columns Vector of column containing containing the largest string for each column
   * @param msgTableError Vector containing all error messages
   */
  void computeTableMaxLineLength( std::vector< TableLayout::Column > & columns,
                                  std::vector< string > & msgTableError ) const;

  /**
   * @brief Build all separator needed from length information contained in columns vector
   * @param topSeparator Top separator to be built
   * @param sectionSeparator section separator to be built
   */
  void buildTableSeparators( std::vector< TableLayout::Column > & columns,
                             string & topSeparator,
                             string & sectionSeparator ) const;

  /**
   * @brief add a row on top of the table
   * @param tableOutput The output stream
   * @param msg Vector of string to display
   * @param topSeparator The top table separator
   * @param sectionSeparator The section table separator
   */
  void addTopRow( std::ostringstream & tableOutput,
                  std::vector< string > const & msg,
                  string_view topSeparator,
                  string_view sectionSeparator ) const;

  /**
   * @brief Add a row on top of the table
   * @param tableOutput The output stream
   * @param msg The string to display
   * @param topSeparator The top table separator
   * @param sectionSeparator The section table separator
   */
  void addTopRow( std::ostringstream & tableOutput,
                  string const & msg,
                  string_view topSeparator,
                  string_view sectionSeparator ) const;

  /**
   * @brief Build a row at the top of the table
   * @param tableOutput The output stream
   * @param msg The string to display.
   * @param topSeparator The top table separator
   * @param sectionSeparator The section table separator
   */
  void buildTopRow( std::ostringstream & tableOutput,
                    string const & msg,
                    string_view topSeparator,
                    string_view sectionSeparator ) const;

  /**
   * @brief Build a section by specifying it's type ( header or section )
   * @param sectionSeparator Line separator between sections
   * @param rows A section row
   * @param nbRows Indicates the number of lines in a section
   * @param section The section to be built
   */
  void buildSectionRows( std::vector< TableLayout::Column > const & columns,
                         string_view sectionSeparator,
                         std::ostringstream & rows,
                         integer const nbRows,
                         TableLayout::Section const section ) const;
};
}

#endif
