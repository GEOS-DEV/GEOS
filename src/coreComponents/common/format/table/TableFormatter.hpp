/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TableFormatter.hpp
 */

#ifndef GEOS_COMMON_FORMAT_TABLE_TABLEFORMATTER_HPP
#define GEOS_COMMON_FORMAT_TABLE_TABLEFORMATTER_HPP

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

  TableFormatter() = default;

  /**
   * @brief Construct a new Table Formatter from a tableLayout
   * @param tableLayout Contain all column names and optionnaly the table title
   */
  TableFormatter( TableLayout const & tableLayout );

  /**
   * @brief Destroy the Table Formatter object
   */
  virtual ~TableFormatter() = default;
};

/**
 * @brief class for CSV formatting
 */
class TableCSVFormatter : public TableFormatter
{
public:

  /**
   * @brief Construct a new Table Formatter
   */
  TableCSVFormatter():
    TableFormatter( TableLayout() )
  {}

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
   * @brief Convert a data source to a CSV string.
   * @tparam DATASOURCE The source to convert
   * @param tableData The data source to convert
   * @return The CSV string representation of a data source.
   */
  template< typename DATASOURCE >
  string toString( DATASOURCE const & tableData ) const;

};

/**
 * @brief Convert the TableData to a CSV string.
 * @param tableData The TableData to convert.
 * @return The CSV string representation of the TableData.
 */
template<>
string TableCSVFormatter::toString< TableData >( TableData const & tableData ) const;


/**
 * @brief class for log formatting
 */
class TableTextFormatter : public TableFormatter
{
public:

  /**
   * @brief Construct a new TableFormatter
   */
  TableTextFormatter():
    TableFormatter( TableLayout() )
  {}

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
   * @return A TableLayout converted into a formatted representation.
   */
  string layoutToString() const;

  /**
   * @brief Convert a data source to a table string.
   * @param tableData The data source to convert.
   * @return The table string representation of the TableData.
   */
  template< typename DATASOURCE >
  string toString( DATASOURCE const & tableData ) const;

private:

  /// symbol for separator construction
  static constexpr char m_verticalLine = '|';
  ///  for the extremity of a row
  static constexpr char m_horizontalLine = '-';

  /**
   * @brief Prepare all the columns with the appropriate values and formatting separator
   * @param columns The vector containg all columns. Each column contains its own
   *        parameters (such as name, alignment, etc.) and column values.
   * @param tableData Vector containing all rows filled with values
   * @param nbHeaderRows A variable to be calculated which will contain the number of header lines to be displayed
   * @param sectionSeparatingLine The row table separator
   * @param topSeparator The top table separator
   */
  void prepareAndBuildTable( std::vector< TableLayout::Column > & columns,
                             TableData const & tableData,
                             size_t & nbHeaderRows,
                             std::string & sectionSeparatingLine,
                             std::string & topSeparator ) const;
/**
 * @brief Displays the complete table
 * @param tableOutput The output stream
 * @param columns The vector containg all columns
 * @param tableData Vector containing all rows filled with values
 * @param nbHeaderRows A variable to be calculated which will contain the number of header lines to be displayed
 * @param sectionSeparatingLine The row table separator
 * @param topSeparator The top table separator
 */
  void outputTable( std::ostringstream & tableOutput,
                    std::vector< TableLayout::Column > & columns,
                    TableData const & tableData,
                    size_t & nbHeaderRows,
                    std::string const & sectionSeparatingLine,
                    std::string const & topSeparator ) const;

  /**
   * @brief Populate all the column values  with values extracted the table data.
   * @param columns The vector containg all columns
   * @param tableData Vector containing all rows filled with values
   * @param isSubColumn A boolean indicating whether the current column is a subcolumn.
   *
   * @note Use transposition for populating column with no subColumn
   */
  void populateColumnsFromTableData( std::vector< TableLayout::Column > & columns,
                                     std::vector< std::vector< string > > const & tableData,
                                     bool isSubColumn ) const;

  /**
   * @brief Split all header names by detecting the newline \\n character. and
   * set the same vector size for each split header and merge it into columns
   * @param columns The vector containg all columns
   * @param nbHeaderRows Variable which will contain the number of header lines to be displayed
   * @param splitHeaders A empty vector who will contain all split header names
   */
  void splitAndMergeColumnHeaders( std::vector< TableLayout::Column > & columns,
                                   size_t & nbHeaderRows,
                                   std::vector< std::vector< string > > & splitHeaders ) const;

  /**
   * @brief For each column find and set the column's longest string
   * @param column The current column to be processed
   * @param maxStringSize Store the longest string(s) for each column.
   * @param idxColumn The current index of the column
   *
   * @note Compares the longest string from the header with the longest string from the column values.
   * If the column contains subcolumns, it recursively applies the same logic to them
   */
  void findAndSetLongestColumnString( TableLayout::Column & column,
                                      std::vector< std::string > & maxStringSize,
                                      integer const idxColumn ) const;

  /**
   * @brief Compute the max table line length, taking into account the length of : title, error, columns header and content
   * Increase the size of the columns if necessary
   * @param columns Vector of column containing containing the largest string for each column
   */
  void computeTableWidth( std::vector< TableLayout::Column > & columns ) const;

  /**
   * @brief Increase each column size if the title is larger than all the columns
   * @param columns Vector containing all table columns
   * @param extraCharacters ExtraCharacters to be distributed between each columns
   */
  void increaseColumnsSize( std::vector< TableLayout::Column > & columns,
                            integer const extraCharacters ) const;

  /**
   * @brief Build all separators needed from content length contained in the columns vector
   * @param columns Vector containing all table columns
   * @param sectionSeparatingLine The row table separator
   * @param topSeparator The top table separator
   */
  void buildTableSeparators( std::vector< TableLayout::Column > const & columns,
                             std::string & sectionSeparatingLine,
                             std::string & topSeparator ) const;

  /**
   * @brief Output the values rows in the table
   * @param columns Vector containing all table columns
   * @param tableOutput The output stream
   * @param nbRows The total number of rows to output.
   * @param sectionSeparatingLine The row table separator
   */
  void outputValuesSectionRows( std::vector< TableLayout::Column > const & columns,
                                std::ostringstream & tableOutput,
                                size_t const nbRows,
                                std::string const & sectionSeparatingLine ) const;

  /**
   * @brief utput the title row in the table
   * @param topSeparator The top separator string
   */
  void outputTitleRow( std::ostringstream & tableOutput,
                       std::string const & topSeparator ) const;

  /**
   * @brief Output the header rows in the table
   * @param columns Vector containing all table columns
   * @param tableOutput The output stream
   * @param nbRows The total number of rows to output.
   * @param sectionSeparatingLine The row table separator
   */
  void outputHeaderSectionRows( std::vector< TableLayout::Column > const & columns,
                                std::ostringstream & tableOutput,
                                size_t const nbRows,
                                std::string const & sectionSeparatingLine ) const;

  /**
   * @brief Output the sub column row
   * @param columns Vector containing the subcolumn values
   * @param tableOutput The output stream
   * @param The current row in the table
   */
  void outputSubSection( std::vector< TableLayout::Column > const & columns,
                         std::ostringstream & tableOutput,
                         integer idxRow ) const;
};

/**
 * @brief Convert a TableData to a table string.
 * @param tableData The TableData to convert.
 * @return The table string representation of the TableData.
 */
template<>
string TableTextFormatter::toString< TableData >( TableData const & tableData ) const;
}

#endif /* GEOS_COMMON_FORMAT_TABLE_TABLEFORMATTER_HPP */
