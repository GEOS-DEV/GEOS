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
   * @param tableLayout Contain all tableColumnData names and optionnaly the table title
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
   * @param tableLayout Contain all tableColumnData names and optionnaly the table title
   */
  TableCSVFormatter( TableLayout const & tableLayout );

  /**
   * @brief Destroy the TableCSVFormatter object
   */
  virtual ~TableCSVFormatter() = default;

  /**
   * @return The string with all tableColumnData names.
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
   * @param tableLayout Contain all tableColumnData names and optionnaly the table title
   */
  TableTextFormatter( TableLayout const & tableLayout );


  /**
   * @brief Destroy the Table Text Formatter object
   */
  virtual ~TableTextFormatter() = default;

  /**
   * @return A TableLayout string representation,
   * The TableTextFormatter receives hasn't receive any data, so only the top part is returned.
   */
  string toString() const;

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
   * @brief Prepare all the tableColumnsData  with the appropriate values and formatting separator
   * @param tableColumnsData  The vector containg all tableColumnsData . Each tableColumnData contains its own
   *        parameters (such as name, alignment, etc.) and column values.
   * @param tableData Vector containing all rows filled with values
   * @param nbHeaderRows Number of header rows, which will be calculated based on tableColumnData headers and their formatting.
   * @param sectionSeparatingLine Separator string used between sections of the table
   * @param topSeparator The table top separator
   */
  void prepareAndBuildTable( std::vector< TableLayout::TableColumnData > & tableColumnsData,
                             TableData const & tableData,
                             size_t & nbHeaderRows,
                             string & sectionSeparatingLine,
                             string & topSeparator ) const;
/**
 * @brief Displays the complete table
 * @param tableOutput The output stream
 * @param tableColumnsData  Vector containg all tableColumnsData
 * @param tableData Vector containing all rows filled with values
 * @param nbHeaderRows A variable to be calculated which will contain the number of header lines
 * @param sectionSeparatingLine Separator string used between sections of the table
 * @param topSeparator The table top separator
 */
  void outputTable( std::ostringstream & tableOutput,
                    std::vector< TableLayout::TableColumnData > & tableColumnsData,
                    TableData const & tableData,
                    size_t & nbHeaderRows,
                    string_view sectionSeparatingLine,
                    string_view topSeparator ) const;

  /**
   * @brief Populate all the tableColumnData values with values extracted from TableData.
   * @param tableColumnsData  Vector of tableColumnsData  to populate.
   * @param tableData Vector containing all rows filled with values
   * @param isSubColumn Boolean indicating if the current tableColumnData is a subcolumn
   */
  void populateColumnsFromTableData( std::vector< TableLayout::TableColumnData > & tableColumnsData,
                                     std::vector< std::vector< string > > const & tableData,
                                     bool isSubColumn ) const;

  /**
   * @brief Split all header names by detecting the newline \\n character. and
   * set the same vector size for each split header and merge it into tableColumnsData
   * @param tableColumnsData  The vector containg all tableColumnsData
   * @param nbHeaderRows Variable which will contain the number of header lines to be displayed
   * @param splitHeaders Vector to store the split header names for each tableColumnData
   */
  void splitAndMergeColumnHeaders( std::vector< TableLayout::TableColumnData > & tableColumnsData,
                                   size_t & nbHeaderRows,
                                   std::vector< std::vector< string > > & splitHeaders ) const;

  /**
   * @brief For each tableColumnData find and set the column's longest string
   * @param tableColumnData The tableColumnData to process.
   * @param maxStringSize Store the longest string(s) for each tableColumnData.
   * @param idxColumn The current index of the tableColumnData
   *
   * @note Compares the longest string from the header with the longest string from the column values.
   * If the column contains subcolumns, it recursively applies the same logic to them
   */
  void findAndSetLongestColumnString( TableLayout::TableColumnData & tableColumnData,
                                      std::vector< string > & maxStringSize,
                                      integer const idxColumn ) const;

  /**
   * @brief Compute the max table line length, taking into account the length of : title, tableColumnsData  header/values
   * Increase the size of the tableColumnsData  if necessary
   * @param tableColumnsData Vector of tableColumnData containing containing the largest string for each tableColumnData
   */
  void computeTableWidth( std::vector< TableLayout::TableColumnData > & tableColumnsData ) const;

  /**
   * @brief Increase each tableColumnData size if the title is larger than all the tableColumnsData
   * @param tableColumnsData  Vector containing all table tableColumnsData
   * @param extraCharacters ExtraCharacters to be distributed between each tableColumnsData
   */
  void increaseColumnsSize( std::vector< TableLayout::TableColumnData > & tableColumnsData,
                            real64 const extraCharacters ) const;

  /**
   * @brief Builds the table's separating lines based on the content length of the tableColumnsData .
   * @param tableColumnsData  Vector containing all table tableColumnsData
   * @param sectionSeparatingLine Separator string used between sections of the table
   * @param topSeparator The table top separator
   */
  void buildTableSeparators( std::vector< TableLayout::TableColumnData > const & tableColumnsData,
                             string & sectionSeparatingLine,
                             string & topSeparator ) const;

  /**
   * @brief Output the values rows in the table
   * @param tableColumnsData  Vector containing all table tableColumnsData
   * @param tableOutput The output stream
   * @param nbRows Total number of rows to output.
   * @param sectionSeparatingLine Separator string used between sections of the table
   */
  void outputValuesSectionRows( std::vector< TableLayout::TableColumnData > const & tableColumnsData,
                                std::ostringstream & tableOutput,
                                size_t const nbRows,
                                string_view sectionSeparatingLine ) const;

  /**
   * @brief Output the title row in the table
   * @param topSeparator The top separator string
   */
  void outputTitleRow( std::ostringstream & tableOutput,
                       string_view topSeparator ) const;

  /**
   * @brief Output the header rows in the table
   * @param tableColumnsData  Vector containing all table tableColumnsData
   * @param tableOutput The output stream
   * @param nbRows The total number of rows to output.
   * @param sectionSeparatingLine Separator string used between sections of the table
   */
  void outputHeaderSectionRows( std::vector< TableLayout::TableColumnData > const & tableColumnsData,
                                std::ostringstream & tableOutput,
                                size_t const nbRows,
                                string_view sectionSeparatingLine ) const;

  /**
   * @brief Outputs subcolumns for the given row in the table.
   * @param tableColumnsData  Vector containing the subcolumn values
   * @param tableOutput The output stream
   * @param idxRow Index of the current row in the table
   */
  void outputSubSection( std::vector< TableLayout::TableColumnData > const & tableColumnsData,
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
