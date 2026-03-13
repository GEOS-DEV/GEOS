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
 * @file TableMpiComponents.cpp
 */

#include "TableMpiComponents.hpp"
#include "common/MpiWrapper.hpp"

namespace geos
{

TableTextMpiOutput::TableTextMpiOutput( TableMpiLayout mpiLayout ):
  TableTextFormatter(),
  m_mpiLayout( mpiLayout )
{}

TableTextMpiOutput::TableTextMpiOutput( TableLayout const & tableLayout,
                                        TableMpiLayout mpiLayout ):
  TableTextFormatter( tableLayout ),
  m_mpiLayout( mpiLayout )
{}

void TableTextMpiOutput::stretchColumnsByRanks( stdVector< size_t > & columnsWidth,
                                                TableTextMpiOutput::Status const & status ) const
{
  { // we ensure we have the correct amount of columns on all ranks (for correct MPI reduction operation)
    size_t const rankColumnsCount = columnsWidth.size();
    size_t const maxRanksColumnsCount = MpiWrapper::max( rankColumnsCount );

    // TODO: contribute to the new table error system with this one
    if( status.m_isContributing )
      GEOS_ASSERT_EQ( rankColumnsCount, maxRanksColumnsCount );

    columnsWidth.resize( maxRanksColumnsCount, 0 );
  }

  // the ranks that does not contribute must not interfere in the column width computing
  if( !status.m_isContributing )
    std::fill( columnsWidth.begin(), columnsWidth.end(), 0 );

  // we keep the largest column widths so we have the same layout on all ranks
  MpiWrapper::allReduce( columnsWidth, columnsWidth, MpiWrapper::Reduction::Max );
}

void TableTextMpiOutput::convertRowsToString( CellLayoutRows const & rows,
                                              stdVector< string > & rowsConvertedInString ) const
{
  std::ostringstream rowStringStream;
  size_t rowIndex = 0;

  for( CellLayoutRow const & row : rows )
  {
    outputLine( m_tableLayout, rows, row, rowStringStream, rowIndex );
    rowIndex++;
    rowsConvertedInString.emplace_back( rowStringStream.str());
    rowStringStream.str( "" );
  }
}

stdVector< TableData::CellData > TableTextMpiOutput::reconstructRow( string_view str ) const
{
  if( str.size() < 1 )
    return stdVector< TableData::CellData >{};

  string_view rowContent( str.substr( 1 )); // skip leading '|'
  string cell;
  stdVector< TableData::CellData > reconstructedRow;

  std::string::size_type start = 0;
  std::string::size_type end = 0;

  while( (end = rowContent.find( m_verticalLine, start )) != string_view::npos )
  {
    cell =std::string( stringutilities::trimSpaces( rowContent.substr( start, end )));
    reconstructedRow.emplace_back( TableData::CellData( {CellType::Value, cell} ));
    rowContent =rowContent.substr( end + 1, rowContent.size());
  }
  return reconstructedRow;
}

void TableTextMpiOutput::gatherAndSortTableDataAcrossRanks ( TableData & gatheredTableData,
                                                             stdVector< string > & rowsAsString,
                                                             TableTextMpiOutput::Status & status ) const
{
  array1d< integer > rowsSizeAcrossAllRank( MpiWrapper::commSize());
  MpiWrapper::allGather( LvArray::integerConversion< integer >( rowsAsString.size()),
                         rowsSizeAcrossAllRank );
  integer const maxRowAcrossAllRanks = *std::max_element( rowsSizeAcrossAllRank.begin(),
                                                          rowsSizeAcrossAllRank.end());
  rowsAsString.resize( maxRowAcrossAllRanks, "" );

  for( string & row : rowsAsString )
  {
    MpiWrapper::gatherStringOnRank0( row, std::function< void(string_view) >( [&]( string_view str ){
      status.m_hasContent = true;
      gatheredTableData.addRow( reconstructRow( str ));
    } ));
  }

  std::sort( gatheredTableData.getCellsData().begin(),
             gatheredTableData.getCellsData().end(),
             *m_sortingFunctor );
}

void TableTextMpiOutput::gatherSortAndOutput( std::ostream & tableOutput,
                                              CellLayoutRows const & rows,
                                              TableTextMpiOutput::Status & status ) const
{
  stdVector< string > rowsAsString;
  convertRowsToString( rows, rowsAsString );

  TableData gatheredTableData;
  gatherAndSortTableDataAcrossRanks( gatheredTableData, rowsAsString, status );

  if( status.m_isMasterRank && status.m_hasContent )
  {
    tableOutput << toString( gatheredTableData );
  }
}

void TableTextMpiOutput::gatherAndOutputTableDataInRankOrder( std::ostream & tableOutput,
                                                              CellLayoutRows const & rows,
                                                              PreparedTableLayout const & tableLayout,
                                                              TableTextMpiOutput::Status & status ) const
{
  // master rank does the output directly to the output, other ranks will have to send it through a string.
  std::ostringstream localStringStream;
  std::ostream & rankOutput = status.m_isMasterRank ? tableOutput : localStringStream;

  if( status.m_isContributing )
  {
    if( m_mpiLayout.m_separatorBetweenRanks )
    {
      string const rankSepLine = GEOS_FMT( "{:-^{}}", m_mpiLayout.m_rankTitle, status.m_sepLine.size() - 2 );
      rankOutput << tableLayout.getIndentationStr() << m_verticalLine << rankSepLine << m_verticalLine << '\n';
    }
    outputTableData( rankOutput, tableLayout, rows );
  }
  string const rankStr = !status.m_isMasterRank && status.m_isContributing ? localStringStream.str() : "";
  stdVector< string > strsAccrossRanks;

  MpiWrapper::gatherStringOnRank0( rankStr, std::function< void(string_view) >( [&]( string_view str ){
    status.m_hasContent = true;
    strsAccrossRanks.emplace_back( str );
  } ));

  if( status.m_isMasterRank && status.m_hasContent )
  {
    for( string_view str : strsAccrossRanks )
      tableOutput << str;
  }
}

template<>
void TableTextMpiOutput::toStream< TableData >( std::ostream & tableOutput,
                                                TableData const & tableData ) const
{
  TableTextMpiOutput::Status status {
    // m_isMasterRank (only the master rank does the output of the header && bottom of the table)
    MpiWrapper::commRank() == 0,
    // m_isContributing (some ranks does not have any output to produce)
    !tableData.getCellsData().empty(),
    // m_hasContent
    false,
    // m_sepLine
    ""
  };

  CellLayoutRows headerCellsLayout;
  CellLayoutRows dataRows;
  CellLayoutRows errorRows;
  size_t tableTotalWidth = 0;

  {
    ColumnWidthModifier const columnWidthModifier = [this, status]( stdVector< size_t > & columnsWidth ) {
      stretchColumnsByRanks( columnsWidth, status );
    };
    initalizeTableGrids( m_tableLayout, tableData,
                         headerCellsLayout, dataRows, errorRows,
                         tableTotalWidth, columnWidthModifier );
    status.m_sepLine = string( tableTotalWidth, m_horizontalLine );
  }

  if( m_sortingFunctor )
  {
    gatherSortAndOutput( tableOutput, dataRows, status );
  }
  else
  {
    if( status.m_isMasterRank )
    {
      outputTableHeader( tableOutput, m_tableLayout, headerCellsLayout, status.m_sepLine );
      tableOutput.flush();
    }
    gatherAndOutputTableDataInRankOrder( tableOutput, dataRows, m_tableLayout, status );
    if( status.m_isMasterRank )
    {
      outputTableFooter( tableOutput, m_tableLayout, errorRows,
                         status.m_sepLine, status.m_hasContent );
      tableOutput.flush();
    }

  }
}

} /* namespace geos */
