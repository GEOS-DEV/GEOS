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

TableTextMpiFormatter::TableTextMpiFormatter( TableMpiLayout mpiLayout ):
  TableTextFormatter(),
  m_mpiLayout( mpiLayout )
{}

TableTextMpiFormatter::TableTextMpiFormatter( TableLayout const & tableLayout,
                                              TableMpiLayout mpiLayout ):
  TableTextFormatter( tableLayout ),
  m_mpiLayout( mpiLayout )
{}

void TableTextMpiFormatter::stretchColumnsByRanks( stdVector< size_t > & columnsWidth,
                                                   TableTextMpiFormatter::Status const & status ) const
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

stdVector< TableData::CellData > TableTextMpiFormatter::parseStringRow( string_view rowString ) const
{
  if( rowString.empty() )
    return stdVector< TableData::CellData >{};

  if( rowString.front() == '|' )
    rowString.remove_prefix( 1 );
  string_view rowContent =rowString;
  string cell;
  stdVector< TableData::CellData > dataRow;

  std::string::size_type end = 0;

  while( (end = rowContent.find( m_verticalLine )) != string_view::npos )
  {
    cell =std::string( stringutilities::trimSpaces( rowContent.substr( 0, end )));
    dataRow.emplace_back( TableData::CellData( {CellType::Value, cell} ));
    rowContent.remove_prefix( end + 1 );
  }
  return dataRow;
}

void TableTextMpiFormatter::gatherAndOutputTableDataInRankOrder( std::ostream & tableOutput,
                                                                 TableFormatter::CellLayoutRows const & rows,
                                                                 PreparedTableLayout const & tableLayout,
                                                                 TableTextMpiFormatter::Status & status ) const
{
  // master rank does the output directly to the output, other ranks will have to send it through a string.
  std::ostringstream localStringStream;
  std::ostream & rankOutput = status.m_isMasterRank ? tableOutput : localStringStream;

  if( status.m_isContributing )
  {
    if( m_mpiLayout.m_separatorBetweenRanks )
    {
      size_t const sepWidth = status.m_sepLine.size() > 2 ? status.m_sepLine.size() - 2 : 0;
      string const rankSepLine = GEOS_FMT( "{:-^{}}", m_mpiLayout.m_rankTitle, sepWidth );
      rankOutput << tableLayout.getIndentationStr() << m_verticalLine << rankSepLine << m_verticalLine << '\n';
    }
    outputTableData( rankOutput, tableLayout, rows );
  }
  string const rankStr = !status.m_isMasterRank && status.m_isContributing ? localStringStream.str() : "";
  stdVector< string > strsAccrossRanks;

  MpiWrapper::gatherStringOnRank0( rankStr, std::function< void(string_view) >( [&]( string_view str ){
    status.m_hasContent = true;
    strsAccrossRanks.emplace_back( str );
  } ) );

  if( status.m_isMasterRank && status.m_hasContent )
  {
    for( string_view str : strsAccrossRanks )
      tableOutput << str;
  }
}

TableData TableTextMpiFormatter::gatherTableDataRank0( TableData const & localTableData ) const
{
  stdVector< buffer_unit_type > serializedTableData( 0 );
  size_t totalSize = 0;

  { // allocation
    totalSize = localTableData.getSerializedSize();
    serializedTableData.reserve( totalSize );
  }

  { // Packing
    if( totalSize > 0 )
    {
      localTableData.serialize( serializedTableData );
    }
  }
  auto [globalLogRecords, counts, offsets] =
    MpiWrapper::gatherBufferRank0< stdVector< buffer_unit_type > >( serializedTableData );


  { // Unpacking
    TableData tableDataGathered;
    if( MpiWrapper::commRank() == 0 )
    {
      buffer_unit_type const * startBuff = globalLogRecords.data();
      for( size_t idxRank = 0; idxRank <  (size_t)MpiWrapper::commSize(); ++idxRank )
      {
        integer byteFromThisRank = counts[idxRank];
        buffer_unit_type const * endRowsBuff = startBuff + byteFromThisRank;
        while( startBuff < endRowsBuff )
        {
          size_t byteFromThisRow = 0;
          TableData::deserializeField( byteFromThisRow, startBuff, endRowsBuff );
          buffer_unit_type const * endRowBuff= startBuff + byteFromThisRow;
          stdVector< TableData::CellData > row;
          while( startBuff < endRowBuff )
          {
            CellType cellType;
            TableData::deserializeField( cellType, startBuff, endRowBuff );
            string cellValue;
            TableData::deserializeField( cellValue, startBuff, endRowBuff );
            row.push_back( {cellType, cellValue} );
          }
          tableDataGathered.addRow( row );
        }
      }
    }
    return tableDataGathered;
  }
}

template<>
void TableTextMpiFormatter::toStream< TableData >( std::ostream & tableOutput,
                                                   TableData const & tableData ) const
{
  TableTextMpiFormatter::Status status {
    // m_isMasterRank (only the master rank does the output of the header && bottom of the table)
    MpiWrapper::commRank() == 0,
    // m_isContributing (some ranks does not have any output to produce)
    !tableData.getCellsData().empty(),
    // m_hasContent
    false,
    // m_sepLine
    ""
  };

  if( m_sortingFunctor )
  {
    TableData tableDataGathered = gatherTableDataRank0( tableData );

    if( status.m_isMasterRank )
    {
      tableDataGathered.sort( *m_sortingFunctor );
      TableTextFormatter::toStream( tableOutput, tableDataGathered );
    }
  }
  else
  { // this version is faster (MPI cooperation) but can only be ordered by rank id
    CellLayoutRows headerRows;
    CellLayoutRows dataRows;
    CellLayoutRows errorRows;
    size_t tableTotalWidth = 0;
    { // compute layout
      ColumnWidthModifier const columnWidthModifier = [this, status]( stdVector< size_t > & columnsWidth ) {
        stretchColumnsByRanks( columnsWidth, status );
      };
      initalizeTableGrids( m_tableLayout, tableData,
                           headerRows, dataRows, errorRows,
                           tableTotalWidth, columnWidthModifier );
      status.m_sepLine = string( tableTotalWidth, m_horizontalLine );
    }

    if( status.m_isMasterRank )
    {
      outputTableHeader( tableOutput, m_tableLayout, headerRows, status.m_sepLine );
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
