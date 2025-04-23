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
  using RowsCellInput = std::vector< std::vector< TableData::CellData > >;

  /// Represent a row of the Table (header or values) when structured for formatting
  struct CellLayoutRow
  {
    /// The cell list of the row instance.
    std::vector< TableLayout::CellLayout > cells;

    /// The maximum number of lines in the `cells` texts (no text is considered as one line).
    size_t sublinesCount;
  };

  /// Represent a table section (title + header or values) layout: view on the data and its layout settings.
  using CellLayoutRows = std::vector< CellLayoutRow >;

  /**
   * @return The Errors List object
   */
  TableErrorListing & getErrorsList() const
  { return *m_errors; }


protected:

  /// Layout for a table
  PreparedTableLayout const m_tableLayout;

  std::unique_ptr< geos::TableErrorListing > m_errors = std::make_unique< geos::TableErrorListing >();

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

  /**
   * @brief Implements the actual writing of content to an output stream.
   *        Adds appropriate messages to the error list when the operation fails.
   * @param outputStream The stream to write the content to.
   * @param content The string view containing data to be written.
   */
  void toStreamImpl( std::ostream & outputStream, string_view content ) const;
};

/**
 * @brief class for CSV formatting
 */
class TableCSVFormatter final : public TableFormatter
{
public:

  /**
   * @brief The column separator for the CSV output.
   */
  static constexpr string_view m_separator = ",";

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
   * @brief Destroy the Table CSV Formatter object
   * We launch GEOS_WARNING if we have encountered any errors
   */
  ~TableCSVFormatter();

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
   * @tparam DATASOURCE The type of the source to convert
   * @param tableData The data source to convert
   * @return The CSV string representation of a data source.
   */
  template< typename DATASOURCE >
  string toString( DATASOURCE const & tableData ) const;

  /**
   * @brief Output the formatted data to a stream. Adds appropriate messages to the error list when the operation fails.
   * @see toString( DATASOURCE const & tableData )
   * @param outputStream The stream to write the content to.
   */
  void headerToStream( std::ostream & outputStream ) const
  { toStreamImpl( outputStream, headerToString() ); }

  /**
   * @brief Output the formatted data to a stream. Adds appropriate messages to the error list when the operation fails.
   * @see toString( DATASOURCE const & tableData )
   * @param tableData The table data
   * @param outputStream The stream to write the content to.
   */
  void dataToStream( std::ostream & outputStream, TableData const & tableData ) const
  { toStreamImpl( outputStream, dataToString( tableData ) ); }

  /**
   * @brief Output the formatted data to a stream. Adds appropriate messages to the error list when the operation fails.
   * @see toString( DATASOURCE const & tableData )
   * @tparam DATASOURCE The source to convert
   * @param tableData The data source to convert
   * @param outputStream The stream to write the content to.
   */
  template< typename DATASOURCE >
  void toStream( std::ostream & outputStream, DATASOURCE const & tableData ) const
  { toStreamImpl( outputStream, toString( tableData ) ); }

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
class TableTextFormatter final : public TableFormatter
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
   * The TableTextFormatter receives hasn't received any data, so only the header part is returned.
   */
  string toString() const;

  /**
   * @brief Convert a data source to a table string.
   * @tparam DATASOURCE The type of the source to convert
   * @param tableData The data source to convert.
   * @return The table string representation of the TableData.
   */
  template< typename DATASOURCE >
  string toString( DATASOURCE const & tableData ) const;

  /**
   * @brief Output the formatted data to a stream. Adds appropriate messages to the error list when the operation fails.
   * @see toString()
   * @param outputStream The stream to write the content to.
   */
  void toStream( std::ostream & outputStream ) const
  { toStreamImpl( outputStream, toString() ); }

  /**
   * @brief Output the formatted data to a stream. Adds appropriate messages to the error list when the operation fails.
   * @see toString( DATASOURCE const & tableData )
   * @tparam DATASOURCE The type of the source to convert
   * @param tableData The data source to convert.
   * @param outputStream The stream to write the content to.
   */
  template< typename DATASOURCE >
  void toStream( std::ostream & outputStream, DATASOURCE const & tableData ) const
  { toStreamImpl( outputStream, toString( tableData ) ); }

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
   * @param errorCellsLayout A reference to a `CellLayoutRows` where the error cells will be populated.
   * @param separatorLine A string that will be used as the table separator line
   */
  void initalizeTableGrids( PreparedTableLayout const & tableLayout,
                            TableData const & tableData,
                            CellLayoutRows & dataCellsLayout,
                            CellLayoutRows & headerCellsLayout,
                            CellLayoutRows & errorCellsLayout,
                            size_t & tableTotalWidth ) const;

  /**
   * @brief Outputs the formatted table to the provided output stream.
   * @param tableLayout The layout of the table
   * @param tableOutput A reference to an `std::ostringstream` where the formatted table will be written.
   * @param headerCellsLayout The layout of the header rows
   * @param dataCellsLayout The layout of the data rows
   * @param errorCellsLayout The layout of the error rows
   * @param separatorLine The string to be used as the table separator line
   */
  void outputTable( PreparedTableLayout const & tableLayout,
                    std::ostringstream & tableOutput,
                    CellLayoutRows const & headerCellsLayout,
                    CellLayoutRows const & dataCellsLayout,
                    CellLayoutRows & errorCellsLayout,
                    size_t tableTotalWidth ) const;

  /**
   * @brief Outputs the formatted table lines to the output stream.
   * @param tableLayout The layout of the table
   * @param cellsLayout A collection of rows, each containing a layout of cells to be processed and formatted.
   * @param tableOutput A reference to an `std::ostringstream` where the formatted table will be written.
   */
  void outputLines( PreparedTableLayout const & tableLayout,
                    CellLayoutRows const & cellsLayout,
                    std::ostringstream & tableOutput ) const;

  /**
   * @brief Outputs the formatted table lines to the output stream.
   * @param tableLayout The layout of the table
   * @param errorCellsLayout The layout of the error rows
   * @param tableOutput A reference to an `std::ostringstream` where the formatted table will be written.
   */
  void outputErrors( PreparedTableLayout const & tableLayout,
                     CellLayoutRows & errorCellsLayout,
                     std::ostringstream & tableOutput ) const;

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
                                  CellLayoutRows & headerCellsLayout ) const;
  /**
   * @brief Populates the data cells layout based on input data values, as a free layout (no columns layout).
   * @param tableLayout The layout of the table, containing information about columns, headers, and their layers.
   * @param dataCellsLayout A reference to the layout for the data cells that will be populated.
   * @param inputDataValues A 2D vector containing the actual input data values.
   */
  void populateDataCellsLayout( PreparedTableLayout const & tableLayout,
                                CellLayoutRows & dataCellsLayout,
                                RowsCellInput const & inputDataValues ) const;

  /**
   * @brief Populates the error cells layout based on input error values
   * @param tableLayout The layout of the table, containing information about columns, headers, and their layers.
   * @param errorCellsLayout A reference to the layout for the error cells that will be populated.
   * @param tableData A constant reference to the `TableData` object, which contains the actual data for the table.
   */
  void populateErrorCellsLayout( PreparedTableLayout const & tableLayout,
                                 CellLayoutRows & errorCellsLayout,
                                 TableData const &  tableData ) const;

  /**
   * @brief Populates the data cells layout based on input data values, taking into account the columns layout.
   * @param tableLayout The layout of the table, containing information about columns, headers, and their layers.
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
  void stretchColumnsByCellsWidth( std::vector< size_t > & columnsWidth,
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
  void stretchColumnsByMergedCellsWidth( std::vector< size_t > & columnsWidth,
                                         TableFormatter::CellLayoutRows & tableGrid,
                                         PreparedTableLayout const & tableLayout,
                                         bool const compress ) const;

  /**
   * @brief Applies column widths to all rows in the table grid.
   * @param columnsWidth The row containing the finalized column width values.
   * @param tableGrid The grid of cells that will have widths propagated to all rows.
   * @param tableLayout Layout information including spacing and other display settings.
   */
  void applyColumnsWidth( std::vector< size_t > const & columnsWidth,
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
