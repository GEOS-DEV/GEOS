/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
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
   * @brief
   * @param columns
   * @param nbHeaderRows
   */
  void prepareAndBuildTable( std::vector< TableLayout::Column > & columns,
                             TableData const & tableData,
                             size_t & nbHeaderRows,
                             std::string & sectionSeparatingLine,
                             std::string & topSeparator ) const;
/**
 * @brief
 *
 * @param tableOutput
 * @param columns
 * @param nbHeaderRows
 */
  void outputTable( std::ostringstream & tableOutput,
                    std::vector< TableLayout::Column > & columns,
                    TableData const & tableData,
                    size_t & nbHeaderRows,
                    std::string const & sectionSeparatingLine,
                    std::string const & topSeparator ) const;

  /**
   * @brief Fill the vector (m_column) in tableData with values from rows stored in tableData.
   * @param columns Vector of columns to be filled.
   * @param tableData Vector containing all rows filled with values
   */
  void populateColumnsFromTableData( std::vector< TableLayout::Column > & columns,
                                     std::vector< std::vector< string > > const & tableData,
                                     bool isSubColumn ) const;

  /**
   * @brief Split all header names by detecting the newline \\n character. and
   * set the same vector size for each split header and merge it into columns
   * @param columns The vector containg all columns
   * @param largestHeaderVectorSize The largest split header vector size
   * @param splitHeaders A empty vector who will contain all split header names
   */
  void splitAndMergeColumnHeaders( std::vector< TableLayout::Column > & columns,
                                   size_t & largestHeaderVectorSize,
                                   std::vector< std::vector< string > > & splitHeaders ) const;

  /**
   * @brief For each column find and set the column's longest string
   * @param maxStringSize The largest string(s) in the column
   * @param columns The vector containg all columns
   * @param idxColumn The current index of the column
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
   * @brief Build all separators needed from content length contained in the columns vector
   * @param columns Vector containing all table columns
   */
  void buildTableSeparators( std::vector< TableLayout::Column > const & columns,
                             std::string & sectionSeparatingLine,
                             std::string & topSeparator ) const;

  /**
   * @brief Add a row on top of the table
   * @param tableOutput The output stream
   * @param msg Vector of string(s) to display
   * @param topSeparator The top table separator
   * @param alignment The aligment for a row
   */
  void outputTopRows( std::ostringstream & tableOutput,
                      std::vector< string > const & msg,
                      string_view topSeparator,
                      TableLayout::Alignment alignment,
                      string_view sectionSeparatingLine ) const;

  /**
   * @brief Output a section by specifying it's type ( header or section )
   * @param columns Vector containing all table columns
   * @param sectionSeparatingLine Line separator between sections
   * @param rows A section row
   * @param nbRows Indicates the number of lines in a section
   * @param section The section to be built
   * @note Add the ending line if there are one or more rows
   */
  void outputValuesSectionRows( std::vector< TableLayout::Column > const & columns,
                                std::ostringstream & tableOutput,
                                size_t const nbRows,
                                std::string const & sectionSeparatingLine ) const;

  /**
   * @brief Construct a new output Title Row object
   * @param topSeparator
   */
  void outputTitleRow( std::ostringstream & tableOutput,
                       std::string const & topSeparator ) const;

  /**
   * @brief
   *
   * @param columns
   * @param tableOutput
   * @param nbRows
   * @param sectionSeparatingLine
   */
  void outputHeaderSectionRows( std::vector< TableLayout::Column > const & columns,
                                std::ostringstream & tableOutput,
                                size_t const nbRows,
                                std::string const & sectionSeparatingLine ) const;

  /**
   * @brief
   *
   * @param columns
   * @param tableOutput
   */
  void outputSubSection( std::vector< TableLayout::Column > const & columns,
                         std::ostringstream & tableOutput,
                         integer idxCell ) const;
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
