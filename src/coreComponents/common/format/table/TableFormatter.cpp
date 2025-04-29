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
 * @file TableFormatter.cpp
 */

#include "TableFormatter.hpp"
#include <numeric>
#include "common/format/StringUtilities.hpp"
#include "common/logger/Logger.hpp"
#include "TableFormatter.hpp"

namespace geos
{

TableFormatter::TableFormatter():
  m_tableLayout()
{}

TableFormatter::TableFormatter( TableLayout const & tableLayout ):
  m_tableLayout( tableLayout )
{}

///////////////////////////////////////////////////////////////////////
////// CSV Formatter implementation
///////////////////////////////////////////////////////////////////////

TableCSVFormatter::TableCSVFormatter( TableLayout const & tableLayout ):
  TableFormatter( tableLayout )
{}

static constexpr string_view csvSeparator = ",";
string TableCSVFormatter::headerToString() const
{
  string result;

  size_t total_size = 0;
  for( auto const & column : m_tableLayout.getColumns())
  {
    for( auto const & str : column.m_header.m_layout.getLines() )
    {
      total_size += str.size();
    }
    total_size += csvSeparator.size();
  }
  result.reserve( total_size );

  for( std::size_t idxColumn = 0; idxColumn < m_tableLayout.getColumns().size(); ++idxColumn )
  {
    std::ostringstream strValue;
    for( auto const & str :  m_tableLayout.getColumns()[idxColumn].m_header.m_layout.getLines() )
    {
      result.append( str );
    }

    if( idxColumn < m_tableLayout.getColumns().size() - 1 )
    {
      result.append( csvSeparator );
    }
  }
  result.append( "\n" );
  return result;
}

string TableCSVFormatter::dataToString( TableData const & tableData ) const
{

  RowsCellInput const rowsValues( tableData.getTableDataRows() );
  string result;
  size_t total_size = 0;
  for( auto const & row : rowsValues )
  {
    for( auto const & item : row )
    {
      total_size += item.value.size();
    }
    total_size += row.size();
  }

  result.reserve( total_size );

  for( auto const & row : rowsValues )
  {
    stdVector< string > rowConverted;
    for( auto const & item : row )
    {
      std::istringstream strStream( item.value );
      string line;
      bool detectNewLine = false;
      while( getline( strStream, line, '\n' ))
      {
        rowConverted.push_back( line );
        detectNewLine = true;
      }

      if( !detectNewLine )
        rowConverted.push_back( item.value );
    }
    result.append( stringutilities::join( rowConverted.cbegin(), rowConverted.cend(), csvSeparator ));
    result.append( "\n" );
  }

  return result;
}

template<>
string TableCSVFormatter::toString< TableData >( TableData const & tableData ) const
{
  return headerToString() + dataToString( tableData );
}

///////////////////////////////////////////////////////////////////////
////// Log Formatter implementation
///////////////////////////////////////////////////////////////////////

TableTextFormatter::TableTextFormatter( TableLayout const & tableLayout ):
  TableFormatter( tableLayout )
{}

string TableTextFormatter::toString() const
{
  TableData tableData;
  return toString( tableData );
}

template<>
string TableTextFormatter::toString< TableData >( TableData const & tableData ) const
{
  std::ostringstream tableOutput;
  CellLayoutRows headerCellsLayout;
  CellLayoutRows dataCellsLayout;
  size_t tableTotalWidth = 0;

  initalizeTableGrids( m_tableLayout, tableData,
                       headerCellsLayout, dataCellsLayout,
                       tableTotalWidth );
  outputTable( m_tableLayout, tableOutput,
               headerCellsLayout, dataCellsLayout,
               tableTotalWidth );
  return tableOutput.str();
}

void TableTextFormatter::initalizeTableGrids( PreparedTableLayout const & tableLayout,
                                              TableData const & tableInputData,
                                              CellLayoutRows & headerCellsLayout,
                                              CellLayoutRows & dataCellsLayout,
                                              size_t & tableTotalWidth ) const
{
  bool const hasColumnLayout = tableLayout.getColumnLayersCount() > 0;
  RowsCellInput const & inputDataValues( tableInputData.getTableDataRows() );
  size_t const inputDataRowsCount = !inputDataValues.empty() ? inputDataValues.front().size() : 0;
  // this array will store the displayed width of all columns (it will be scaled by data & headers width)
  stdVector< size_t > columnsWidth;

  populateTitleCellsLayout( tableLayout, headerCellsLayout );
  if( hasColumnLayout )
  {
    populateHeaderCellsLayout( tableLayout, headerCellsLayout, inputDataRowsCount );
    size_t nbVisibleColumns = headerCellsLayout.back().cells.size();
    populateDataCellsLayout( tableLayout, dataCellsLayout, inputDataValues, nbVisibleColumns );
    columnsWidth = stdVector< size_t >( nbVisibleColumns, 0 );
  }
  else
  {
    populateDataCellsLayout( tableLayout, dataCellsLayout, inputDataValues );
    columnsWidth = stdVector< size_t >( inputDataRowsCount, 0 );
  }

  stretchColumnsByCellsWidth( columnsWidth, headerCellsLayout );
  stretchColumnsByCellsWidth( columnsWidth, dataCellsLayout );

  // only after all cells that are not merge, we can process the merged cells.
  stretchColumnsByMergedCellsWidth( columnsWidth, headerCellsLayout, tableLayout, false );
  stretchColumnsByMergedCellsWidth( columnsWidth, dataCellsLayout, tableLayout, true );

  // the columns width array is now sized after all the table, we can compute the total table width
  tableTotalWidth = tableLayout.getBorderMargin() * 2 + 2;
  for( size_t columnId = 0; columnId < columnsWidth.size(); ++columnId )
  {
    tableTotalWidth += columnsWidth[columnId] +
                       size_t( columnId > 0 ? tableLayout.getColumnMargin() : 0 );
  }

  // we can now propagate the columns width width to all cells
  applyColumnsWidth( columnsWidth, headerCellsLayout, tableLayout );
  applyColumnsWidth( columnsWidth, dataCellsLayout, tableLayout );
}

void TableTextFormatter::populateTitleCellsLayout( PreparedTableLayout const & tableLayout,
                                                   CellLayoutRows & headerCellsLayout ) const
{
  TableLayout::CellLayout const & titleInput = tableLayout.getTitleLayout();
  if( !titleInput.isEmpty() )
  { // if it exists, we add the title, as a first row with all cells merged in one containing the title text
    headerCellsLayout.reserve( headerCellsLayout.size() + 2 );

    // the title row consists in a row of cells merging with the last cell containing the title text
    headerCellsLayout.emplace_back() = {
      stdVector< TableLayout::CellLayout >( tableLayout.getLowermostColumnsCount(),
                                            TableLayout::CellLayout( CellType::MergeNext ) ),   // cells
      titleInput.getHeight(), // sublinesCount
    };
    headerCellsLayout.back().cells.back() = titleInput;

    headerCellsLayout.emplace_back() = {
      stdVector< TableLayout::CellLayout >( tableLayout.getLowermostColumnsCount(),
                                            TableLayout::CellLayout( CellType::Separator ) ),   // cells
      1, // sublinesCount
    };
  }
}

void TableTextFormatter::populateHeaderCellsLayout( PreparedTableLayout const & tableLayout,
                                                    CellLayoutRows & headerCellsLayout,
                                                    size_t const inputDataColumnsCount ) const
{
  using CellLayout = TableLayout::CellLayout;

  // Number of lines per header layer.
  size_t const columnsLayersCount = tableLayout.getColumnLayersCount();
  // number of visible data columns (we fit the number of data columns if no column layout has been specified)
  size_t const lowermostColumnsCount = tableLayout.getLowermostColumnsCount();
  // Number of rows in headerCellsLayout by taking into account title & separators lines
  size_t const headerRowsCount = ( columnsLayersCount ) * 2;
  // lambda to get the header row id in the headerCellsLayout: [title -> separator ->] column layer 0 ->  separator -> ... -> column layer
  // n-1 -> separator
  size_t const previousRowsCount = headerCellsLayout.size();
  auto const getColumnRowId = [=] ( size_t columnLayer ) { return previousRowsCount + columnLayer * 2; };

  // TODO: integrate this error in the table, and use an equality with the visible+non-visible lowermost column count
  // (PreparedTableLayout should have a visible & nonvisible getLowermostColumnsCount() verion)
  if( inputDataColumnsCount > 0 )
    GEOS_ERROR_IF_GT( lowermostColumnsCount, inputDataColumnsCount );

  headerCellsLayout.resize( previousRowsCount + headerRowsCount );
  for( size_t rowId = previousRowsCount; rowId < headerCellsLayout.size(); rowId++ )
  {
    if( ( rowId % 2 ) == 1 )
    { // each even row is a separator
      headerCellsLayout[rowId] = {
        stdVector< CellLayout >( lowermostColumnsCount, CellLayout( CellType::Separator ) ), // cells
        1 // sublinesCount
      };
    }
    else
    { // the others are headers that will be filled
      headerCellsLayout[rowId].cells.reserve( lowermostColumnsCount );
    }
  }

  // number of times we will divide the each headers in the headerCellsLayout (key is the Column ptr)
  std::map< std::ptrdiff_t, size_t > subdivsCount;

  for( auto it = tableLayout.beginDeepFirst(); it != tableLayout.endDeepFirst(); ++it )
  {
    CellLayout const & currentCell = it->m_header.m_layout;
    if( currentCell.m_cellType != CellType::Hidden )
    {
      size_t const layer = it.getCurrentLayer();
      size_t const rowId = getColumnRowId( layer );

      // we subdivide divide parent cells each time we have a cell which has a neightboor to its right
      if( it->hasParent() )
      {
        // we subdivide the parent by the subdiv count of the current (+1 if cell has a neightboor at its right)
        size_t const incrementsCount = subdivsCount[std::ptrdiff_t( it.getPtr() )] + ( it->hasNext() ? 1 : 0 );
        subdivsCount[std::ptrdiff_t( it->getParent() )] += incrementsCount;
      }

      // the current cell layer number of lines must be the max of all cells of this layer
      headerCellsLayout[rowId].sublinesCount = std::max( headerCellsLayout[rowId].sublinesCount,
                                                         currentCell.getHeight() );

      // we add the current cell and all its subdivision
      for( size_t idxColumn = 0; idxColumn < subdivsCount[std::ptrdiff_t( it.getPtr() )]; idxColumn++ )
      {
        CellLayout mergingCell{ CellType::MergeNext, TableLayout::Alignment::center };
        headerCellsLayout[rowId].cells.push_back( mergingCell );
      }
      headerCellsLayout[rowId].cells.push_back( currentCell );

      // we add up empty cells under the current one to fill the space to the data.
      if( !it->hasChild() )
      {
        for( size_t subLayer = layer + 1; subLayer < columnsLayersCount; subLayer++ )
        {
          CellLayout emptyCell{CellType::Header, TableLayout::Alignment::center};
          emptyCell.setWidth( currentCell.getWidth() );
          headerCellsLayout[getColumnRowId( subLayer )].cells.push_back( emptyCell );
        }
      }
    }
  }
}

void TableTextFormatter::populateDataCellsLayout( PreparedTableLayout const & tableLayout,
                                                  CellLayoutRows & dataCellsLayout,
                                                  RowsCellInput const & inputDataValues,
                                                  size_t const nbVisibleColumn ) const
{
  dataCellsLayout = stdVector< CellLayoutRow >{
    inputDataValues.size(),
    {
      stdVector< TableLayout::CellLayout >( nbVisibleColumn, TableLayout::CellLayout() ), // cells
      0 // sublinesCount
    }
  };

  // TODO: error if inputDataValues size is not consistent with visible headers count / columns count

  for( size_t idxRow = 0; idxRow < inputDataValues.size(); ++idxRow )
  {
    CellLayoutRow & outputRow = dataCellsLayout[idxRow];
    size_t maxLinesInRow = 0;
    size_t idxInputColumn = 0;
    size_t idxOutputColumn = 0;
    for( auto columnIt = tableLayout.beginDeepFirst(); columnIt != tableLayout.endDeepFirst(); ++columnIt )
    {
      if( !columnIt->hasChild())
      {   // we take into account only the (enabled) headers that are the nearest to the data
        if( columnIt->isVisible() )
        {
          TableData::CellData const & inputCell = inputDataValues[idxRow][idxInputColumn];
          CellType const type = !columnIt->isVisible() ?
                                CellType::Hidden :
                                inputCell.type;
          string_view value = inputCell.type == CellType::Separator ?
                              string_view( &m_horizontalLine, 1 ) :
                              string_view( inputCell.value );

          TableLayout::CellLayout & outputCell = outputRow.cells[idxOutputColumn];
          outputCell = TableLayout::CellLayout( type, columnIt->m_alignment.valueAlignment );
          outputCell.prepareLayout( value, tableLayout.getMaxColumnWidth() );

          maxLinesInRow  = std::max( maxLinesInRow, outputCell.getHeight() );
          idxOutputColumn++;
        }
        idxInputColumn++;
      }
    }
    outputRow.sublinesCount = maxLinesInRow;
  }
}

void TableTextFormatter::populateDataCellsLayout( PreparedTableLayout const & tableLayout,
                                                  CellLayoutRows & dataCellsLayout,
                                                  RowsCellInput const & inputDataValues ) const
{
  size_t const nbColumns = !inputDataValues.empty() ? inputDataValues[0].size() : 0;
  dataCellsLayout = stdVector< CellLayoutRow >{
    inputDataValues.size(),
    {
      stdVector< TableLayout::CellLayout >( nbColumns, TableLayout::CellLayout() ), // cells
      0 // sublinesCount
    }
  };

  for( size_t idxRow = 0; idxRow < inputDataValues.size(); ++idxRow )
  {
    CellLayoutRow & outputRow = dataCellsLayout[idxRow];
    size_t maxLinesInRow = 0;
    for( size_t idxColumn = 0; idxColumn < nbColumns; ++idxColumn )
    {
      TableData::CellData const & inputCell = inputDataValues[idxRow][idxColumn];
      string_view value = inputCell.type == CellType::Separator ?
                          string_view( &m_horizontalLine, 1 ) :
                          string_view( inputCell.value );

      TableLayout::CellLayout & outputCell = outputRow.cells[idxColumn];
      outputCell = TableLayout::CellLayout( inputCell.type, TableLayout::defaultValueAlignment );
      outputCell.prepareLayout( value, tableLayout.getMaxColumnWidth() );

      maxLinesInRow  = std::max( maxLinesInRow, outputCell.getHeight() );
    }
    outputRow.sublinesCount = maxLinesInRow;
  }
}

void TableTextFormatter::stretchColumnsByCellsWidth( stdVector< size_t > & columnsWidth,
                                                     TableFormatter::CellLayoutRows const & tableGrid ) const
{
  // first, we reduce by column all regular cells in the first row.
  size_t const numColumns = tableGrid.empty() ? 0 : tableGrid[0].cells.size();
  for( TableFormatter::CellLayoutRow const & currentRow : tableGrid )
  {
    auto const & rowCells = currentRow.cells;
    if( rowCells[0].m_cellType != CellType::Separator )
    {
      for( size_t columnId = 0; columnId < numColumns; columnId++ )
      {
        auto const & cell = rowCells[columnId];
        bool const isToMerge = columnId > 0 && rowCells[columnId-1].m_cellType == CellType::MergeNext;
        bool const isContentCell = cell.m_cellType != CellType::MergeNext;
        if( isContentCell && !isToMerge && columnsWidth[columnId] < cell.getWidth() )
          columnsWidth[columnId] = rowCells[columnId].getWidth();
      }
    }
  }
}

void TableTextFormatter::stretchColumnsByMergedCellsWidth( stdVector< size_t > & columnsWidth,
                                                           TableFormatter::CellLayoutRows & tableGrid,
                                                           PreparedTableLayout const & tableLayout,
                                                           bool const compress ) const
{
  // To get consistent results, we must process column by column.
  size_t const numRows = tableGrid.size();
  size_t const numColumns = tableGrid.empty() ? 0 : tableGrid[0].cells.size();
  size_t const spaceBetweenColumns = size_t( tableLayout.getColumnMargin() );
  stdVector< size_t > flexSpaces = stdVector< size_t >( columnsWidth.size(), 0 );

  for( size_t columnId = 0; columnId < numColumns; columnId++ )
  {
    for( size_t rowId = 0; rowId < numRows; rowId++ )
    {
      TableLayout::CellLayout & cell = tableGrid[rowId].cells[columnId];
      bool const isContentCell = cell.m_cellType != CellType::MergeNext;
      bool const isLastToMerge = columnId > 0 && tableGrid[rowId].cells[columnId-1].m_cellType == CellType::MergeNext;
      if( isContentCell && isLastToMerge )
      { // detected cells to merge, accumulating content width & column number & space
        size_t mergedCellsWidth = cell.getWidth();
        size_t mergedColumnsWidth = columnsWidth[columnId] + flexSpaces[columnId];
        size_t mergedCellsCount = 1;
        for( integer mergedId = columnId - 1; mergedId >= 0; --mergedId )
        {
          TableLayout::CellLayout const & leftCell = tableGrid[rowId].cells[mergedId];
          if( leftCell.m_cellType != CellType::MergeNext )
            break;

          mergedColumnsWidth += columnsWidth[mergedId] + flexSpaces[mergedId];
          mergedColumnsWidth += spaceBetweenColumns;
          ++mergedCellsCount;
        }

        integer const overflowingWidth = mergedCellsWidth - mergedColumnsWidth;
        if( overflowingWidth > 0 )
        { // if the merged content width exceeds the available columns width, we balance the resizing over all cells.
          auto [stretchPerColumn, remainingStretch] = std::div( overflowingWidth, mergedCellsCount );
          for( size_t mergedId = size_t( columnId - mergedCellsCount + 1 ); mergedId <= columnId; ++mergedId )
          {
            flexSpaces[mergedId] += size_t( stretchPerColumn );
            flexSpaces[mergedId] += ( remainingStretch-- > 0 ) ? 1 : 0;
          }
        }
        else
        { //if not overflowing, it means the merged cell needs to be stretch to fit the available column space.
          cell.setWidth( cell.getWidth() + size_t( -overflowingWidth ) );
        }
      }
    }
  }

  // compression of flexible space, mandatory for when we merge cells in intersecting column sets.
  if( compress && numRows > 0 )
  {
    for( size_t columnId = 0; columnId < numColumns; columnId++ )
    {
      size_t const currentColumnWidth = columnsWidth[columnId] + flexSpaces[columnId];
      integer oversize = std::numeric_limits< integer >::max();
      for( size_t rowId = 0; rowId < numRows; rowId++ )
      {
        CellLayoutRow & row = tableGrid[rowId];
        TableLayout::CellLayout const & cell = row.cells[columnId];
        bool const isToMergeWithRight = cell.m_cellType == CellType::MergeNext;
        bool const isToMergeWithLeft = ( columnId > 0 ) ? ( row.cells[columnId - 1].m_cellType == CellType::MergeNext ) : false;
        bool const isFirstToMerge = isToMergeWithRight && !isToMergeWithLeft;
        if( isFirstToMerge )
        { // we detected a set of cells to merge, let's compute its total width
          size_t mergedColumnsWidth = currentColumnWidth;
          for( size_t mergedId = columnId + 1; mergedId < numColumns; mergedId++ )
          {
            TableLayout::CellLayout const & mergedCell = row.cells[mergedId];
            mergedColumnsWidth += columnsWidth[mergedId] + flexSpaces[mergedId];
            mergedColumnsWidth += spaceBetweenColumns;
            if( mergedCell.m_cellType != CellType::MergeNext || mergedId == numColumns - 1 )
            {
              // this is the last cell to merge (which contains the actual merged content),
              // we can compute here the space potencially wasted by the flexible space
              integer const potentialOversize = integer( mergedColumnsWidth ) - integer( mergedCell.getWidth() );
              oversize = std::min( oversize, potentialOversize );
              break;
            }
          }
        }
        else if( isToMergeWithLeft )
        { // to prevent any side effects, we only shrink columns at the left of a merge
          oversize = 0;
          break;
        }
        else
        { // no merge here, so any flexible space may be unnecessary
          integer const potentialOversize = flexSpaces[columnId];
          oversize = std::min( oversize, potentialOversize );
        }
      }
      if( oversize > 0 )
        flexSpaces[columnId] = size_t( std::max( 0, integer( flexSpaces[columnId] ) - oversize ) );
    }
  }

  for( size_t columnId = 0; columnId < numColumns; columnId++ )
  {
    // TODO warning if <0 (wrong computation) -> ignore addition
    columnsWidth[columnId] += size_t( flexSpaces[columnId] );
  }
}

void TableTextFormatter::applyColumnsWidth( stdVector< size_t > const & columnsWidth,
                                            TableFormatter::CellLayoutRows & tableGrid,
                                            PreparedTableLayout const & tableLayout ) const
{
  size_t const numRows = tableGrid.size();
  size_t const numColumns = tableGrid.empty() ? 0 : tableGrid[0].cells.size();
  size_t const spaceBetweenColumns = size_t( tableLayout.getColumnMargin() );
  for( size_t rowId = 0; rowId < numRows; rowId++ )
  {
    for( size_t columnId = 0; columnId < numColumns; columnId++ )
    {
      auto & cell = tableGrid[rowId].cells[columnId];
      bool const isToMerge = columnId > 0 && tableGrid[rowId].cells[columnId-1].m_cellType == CellType::MergeNext;
      if( !isToMerge )
      {
        cell.setWidth( std::max( cell.getWidth(), columnsWidth[columnId] ) );
      }
      else
      {
        size_t mergedColumnsWidth = columnsWidth[columnId];
        for( integer previousColumnId = columnId-1; previousColumnId >= 0; --previousColumnId )
        {
          TableLayout::CellLayout const & leftCell = tableGrid[rowId].cells[previousColumnId];
          if( leftCell.m_cellType != CellType::MergeNext )
            break;

          mergedColumnsWidth += columnsWidth[previousColumnId] + spaceBetweenColumns;
        }
        cell.setWidth( mergedColumnsWidth );
      }
    }
  }
}

void TableTextFormatter::outputTable( PreparedTableLayout const & tableLayout,
                                      std::ostringstream & tableOutput,
                                      CellLayoutRows const & headerCellsLayout,
                                      CellLayoutRows const & dataCellsLayout,
                                      size_t const tableTotalWidth ) const
{
  string const sepLine = string( tableTotalWidth, m_horizontalLine );
  if( tableLayout.isLineBreakEnabled())
  {
    tableOutput << '\n';
  }
  tableOutput << sepLine << '\n';
  outputLines( tableLayout, headerCellsLayout, tableOutput );
  if( !dataCellsLayout.empty())
  {
    outputLines( tableLayout, dataCellsLayout, tableOutput );
    tableOutput << sepLine;
  }
  if( tableLayout.isLineBreakEnabled())
  {
    tableOutput << '\n';
  }
}

/**
 * @brief Build cell given an m_alignment, a value and spaces
 * @param m_alignment The aligment of cell value
 * @param value The cell value
 * @param spaces The number of spaces in the cell
 * @return A formated cell
 */
string buildCell( TableLayout::Alignment const m_alignment, string_view value, size_t const spaces )
{
  switch( m_alignment )
  {
    case TableLayout::right:   return GEOS_FMT( "{:>{}}", value, spaces );
    case TableLayout::left:    return GEOS_FMT( "{:<{}}", value, spaces );
    case TableLayout::center:  return GEOS_FMT( "{:^{}}", value, spaces );
    default:                   return GEOS_FMT( "{:>{}}", value, spaces );
  }
}

void TableTextFormatter::formatCell( std::ostringstream & tableOutput,
                                     TableLayout::CellLayout const & cell,
                                     size_t const idxLine ) const
{
  if( cell.m_cellType == CellType::Separator )
  {
    tableOutput << string( cell.getWidth(), m_horizontalLine );
  }
  else
  {
    string_view value = idxLine < cell.getHeight() && cell.m_cellType != CellType::Hidden ?
                        string_view( cell.getLines()[idxLine] ):
                        string_view( "" );
    tableOutput << buildCell( cell.m_alignment, value, cell.getWidth() );
  }
}

void TableTextFormatter::outputLines( PreparedTableLayout const & tableLayout,
                                      CellLayoutRows const & rows,
                                      std::ostringstream & tableOutput ) const
{
  size_t const nbRows = rows.size();
  size_t const nbColumns = !rows.empty() ? rows[0].cells.size() : 0;
  size_t const nbBorderSpaces = tableLayout.getBorderMargin();
  size_t const nbColumnSpaces = ( tableLayout.getColumnMargin() - 1 ) / 2;

  size_t idxRow = 0;
  for( CellLayoutRow const & row : rows )
  {
    for( size_t idxSubLine = 0; idxSubLine < row.sublinesCount; idxSubLine++ )
    {
      bool isLeftBorderCell = true;

      for( size_t idxColumn = 0; idxColumn < nbColumns; ++idxColumn )
      {
        auto & cell = row.cells[idxColumn];
        bool const isRightBorderCell = idxColumn == nbColumns - 1;
        if( cell.m_cellType != CellType::MergeNext || isRightBorderCell )
        {
          bool const isSeparator = cell.m_cellType == CellType::Separator;
          char const cellSpaceChar = isSeparator ? m_horizontalLine : ' ';

          if( isLeftBorderCell )
          { // left table border
            isLeftBorderCell=false;
            tableOutput << m_verticalLine << string( nbBorderSpaces, cellSpaceChar );
          }
          else
          { // left side of a cell that have a neightboor
            tableOutput << string( nbColumnSpaces, cellSpaceChar );
          }

          // cell content / fill
          formatCell( tableOutput, cell, idxSubLine );

          if( !isRightBorderCell )
          { // right side of a cell that have a neightboor
            bool const isNextSeparator = row.cells[idxColumn + 1].m_cellType == CellType::Separator;
            bool const upMerged = idxRow > 0 && rows[idxRow - 1].cells[idxColumn].m_cellType == CellType::MergeNext;
            bool const downMerged = idxRow < nbRows - 1 && rows[idxRow + 1].cells[idxColumn].m_cellType == CellType::MergeNext;
            tableOutput << string( nbColumnSpaces, cellSpaceChar );
            tableOutput << ( isSeparator && isNextSeparator && (upMerged || downMerged) ?
                             m_horizontalLine : m_verticalLine );
          }
          else
          { // right table border
            tableOutput << string( nbBorderSpaces, cellSpaceChar ) << m_verticalLine << "\n";
          }
        }
      }
    }
    idxRow++;
  }
}
}
