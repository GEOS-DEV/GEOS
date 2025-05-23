/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
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
#include "TableTypes.hpp"

namespace geos
{

/**
 * @brief abstract class for formatting table data
 */
class TableFormatter
{

public:
  /// Represent the TableData values
  using RowsCellInput = stdVector< stdVector< TableData::CellData > >;

  /// Represent a row of the Table (header or values) when structured for formatting
  struct CellLayoutRow
  {
    /// The cell list of the row instance.
    stdVector< TableLayout::CellLayout > cells;

    /// The maximum number of lines in the `cells` texts (no text is considered as one line).
    size_t sublinesCount;
  };

  /// Represent a table section (title + header or values) layout: view on the data and its layout settings.
  using CellLayoutRows = stdVector< CellLayoutRow >;


protected:

  /// Layout for a table
  PreparedTableLayout const m_tableLayout;

  /**
   * @brief Construct a default Table Formatter without layout specification (to only insert data in it,
   * without any column / title). Feature is not tested.
   */
  TableFormatter();

  /**
   * @brief Construct a new Table Formatter from a tableLayout
   * @param tableLayout Contain all tableColumnData names and optionnaly the table title
   */
  TableFormatter( TableLayout const & tableLayout );
};

/**
 * @brief class for CSV formatting
 */
class TableCSVFormatter : public TableFormatter
{
public:

  /**
   * @brief Construct a default Table Formatter without layout specification (to only insert data in it,
   * without any column / title). Feature is not tested.
   */
  TableCSVFormatter():
    TableFormatter()
  {}

  /**
   * @brief Construct a new Table Formatter from a tableLayout
   * @param tableLayout Contain all tableColumnData names and optionnaly the table title
   */
  TableCSVFormatter( TableLayout const & tableLayout );

  /**
   * @return The string with all tableColumnData names.
   */
  string headerToString() const;

  /**
   * @brief Convert the table data to a CSV string..
   * @param tableData The table data
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
   * @brief Construct a default Table Formatter without layout specification (to only insert data in it,
   * without any column / title). Feature is not tested.
   */
  TableTextFormatter():
    TableFormatter()
  {}

  /**
   * @brief Construct a new TableFormatter from a tableLayout
   * @param tableLayout Contain all tableColumnData names and optionnaly the table title
   */
  TableTextFormatter( TableLayout const & tableLayout );

  /**
   * @return A TableLayout string representation,
   * The TableTextFormatter receives hasn't receive any data, so only the header part is returned.
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
  /// for the extremity of a row
  static constexpr char m_horizontalLine = '-';


  /**
   * @brief Initializes the table layout with the given table data and prepares necessary layouts for headers and data cells.
   * @param tableLayout A reference to the `TableLayout` object.
   * @param tableData A constant reference to the `TableData` object, which contains the actual data for the table.
   * @param headerCellsLayout A reference to a `CellLayoutRows` where the header cells will be populated.
   * @param dataCellsLayout A reference to a `CellLayoutRows` where the data cells will be populated.
   * @param separatorLine A string that will be used as the table separator line
   */
  void initalizeTableGrids( PreparedTableLayout const & tableLayout,
                            TableData const & tableData,
                            CellLayoutRows & dataCellsLayout,
                            CellLayoutRows & headerCellsLayout,
                            size_t & tableTotalWidth ) const;

  /**
   * @brief Outputs the formatted table to the provided output stream.
   * @param tableLayout The layout of the table
   * @param tableOutput A reference to an `std::ostringstream` where the formatted table will be written.
   * @param headerCellsLayout The layout of the header rows
   * @param dataCellsLayout The layout of the data rows
   * @param separatorLine The string to be used as the table separator line
   */
  void outputTable( PreparedTableLayout const & tableLayout,
                    std::ostringstream & tableOutput,
                    CellLayoutRows const & headerCellsLayout,
                    CellLayoutRows const & dataCellsLayout,
                    size_t tableTotalWidth ) const;

  /**
   * @brief Populate a grid of CellLayout with the title rows.
   * @param tableLayout The layout of the table, containing information about columns, headers, and their layers.
   * @param headerCellsLayout A reference to the collection of header cells that will be updated with the
   *                          gridified layout.
   */
  void populateTitleCellsLayout( PreparedTableLayout const & tableLayout,
                                 CellLayoutRows & headerCellsLayout ) const;

  /**
   * @brief Populate a grid of CellLayout with all visible columns of the given table layout.
   * @note To produce a grid with the given column tree, there are 2 corner cases:
   *       - A column have less subcolumns layers than its neightboors -> empty "Header" cells  will be added bellow.
   *       - A parent column has 2 or more sub-columns -> it will be subdivised with "MergeNext" cells.
   *         This is why stretchColumnsByMergedCellsWidth() must be called on the grid,
   * @param tableLayout The layout of the table, containing information about columns, headers, and their layers.
   * @param headerCellsLayout A reference to the collection of header cells that will be updated with the
   *                          gridified layout.
   * @param inputDataColumnsCount The number of input data columns count, helps verifying the number of column.
   */
  void populateHeaderCellsLayout( PreparedTableLayout const & tableLayout,
                                  CellLayoutRows & headerCellsLayout,
                                  size_t inputDataColumnsCount ) const;
  /**
   * @brief Populates the data cells layout based on input data values, as a free layout (no columns layout).
   * @param tableLayout The layout of the table,
   * @param dataCellsLayout A reference to the layout for the data cells that will be populated.
   * @param inputDataValues A 2D vector containing the actual input data values.
   */
  void populateDataCellsLayout( PreparedTableLayout const & tableLayout,
                                CellLayoutRows & dataCellsLayout,
                                RowsCellInput const & inputDataValues ) const;

  /**
   * @brief Populates the data cells layout based on input data values, taking into account the columns layout.
   * @param tableLayout The layout of the table,
   * @param dataCellsLayout A reference to the layout for the data cells that will be populated.
   * @param inputDataValues A 2D vector containing the actual input data values.
   * @param nbVisibleColumn The number of columns that are not hidden
   */
  void populateDataCellsLayout( PreparedTableLayout const & tableLayout,
                                CellLayoutRows & dataCellsLayout,
                                RowsCellInput const & inputDataValues,
                                size_t nbVisibleColumn ) const;

  /**
   * @brief Expend the columns width to accomodate with the content of all cells that are not merged.
   * @param columnsWidth The array to store the resulting columns width in.
   * @param tableGrid The grid of cells containing content.
   */
  void stretchColumnsByCellsWidth( stdVector< size_t > & columnsWidth,
                                   TableFormatter::CellLayoutRows const & tableGrid ) const;

  /**
   * @brief Adjust cell widths to accommodate merged cells across multiple columns.
   * @param columnsWidth The array to store the resulting columns width in.
   *                     Initialized by stretchColumnsByCellsWidth().
   * @param tableGrid The grid of cells containing content that is potencially merged.
   *                  The merged cells width will be computed.
   * @param tableLayout Layout information, including column margins and other settings.
   * @param compress Enable a final compression pass instead of only expanding widths.
   */
  void stretchColumnsByMergedCellsWidth( stdVector< size_t > & columnsWidth,
                                         TableFormatter::CellLayoutRows & tableGrid,
                                         PreparedTableLayout const & tableLayout,
                                         bool const compress ) const;

  /**
   * @brief Applies column widths to all rows in the table grid.
   * @param columnsWidth The row containing the finalized column width values.
   * @param tableGrid The grid of cells that will have widths propagated to all rows.
   * @param tableLayout Layout information including spacing and other display settings.
   */
  void applyColumnsWidth( stdVector< size_t > const & columnsWidth,
                          TableFormatter::CellLayoutRows & tableGrid,
                          PreparedTableLayout const & tableLayout ) const;


  /**
   * @brief Formats a table cell and appends it to the table output.
   * @param tableLayout The layout of the table
   * @param tableOutput The output stream
   * @param cell The cell to format
   * @param idxLine The current line index used to access the specific content for the cell.
   */
  void formatCell( std::ostringstream & tableOutput,
                   TableLayout::CellLayout const & cell,
                   size_t idxLine ) const;

  /**
   * @brief Outputs the formatted table lines to the output stream.
   * @param tableLayout The layout of the table
   * @param cellsLayout A collection of rows, each containing a layout of cells to be processed and formatted.
   * @param tableOutput The output stream
   * @param nbLinesRow A vector containing the number of sub-lines for each row.
   * @param sectionType The type of the section being processed (Header, Value, etc.).
   * @param separatorLine The table separator line string
   */
  void outputLines( PreparedTableLayout const & tableLayout,
                    CellLayoutRows const & cellsLayout,
                    std::ostringstream & tableOutput ) const;
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
