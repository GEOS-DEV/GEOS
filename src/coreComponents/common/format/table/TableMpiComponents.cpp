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
#include "common/format/StringUtilities.hpp" // todo delete

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

template<>
void TableTextMpiOutput::toStream< TableData >( std::ostream & tableOutput,
                                                TableData const & tableData ) const
{
  TableTextMpiOutputStatus status {
    // m_isMasterRank (only the master rank does the output of the header && bottom of the table)
    MpiWrapper::commRank() == 0,
    // m_isContributing (some ranks does not have any output to produce)
    !tableData.getCellsData().empty(),
    // m_sepLine
    ""
  };

  CellLayoutRows headerCellsLayout;
  CellLayoutRows dataCellsLayout;
  size_t tableTotalWidth = 0;

  {
    ColumnWidthModifier const columnWidthModifier = [this, status]( stdVector< size_t > & columnsWidth ) {
      stretchColumnsByRanks( columnsWidth, status );
    };
    initalizeTableGrids( m_tableLayout, tableData,
                         headerCellsLayout, dataCellsLayout,
                         tableTotalWidth, columnWidthModifier );
  }

  if( status.m_isMasterRank )
  {
    status.m_sepLine = string( tableTotalWidth, m_horizontalLine );
    outputTableHeader( tableOutput, m_tableLayout, headerCellsLayout, status.m_sepLine );
    tableOutput.flush();
  }

  outputTableDataToRank0( tableOutput, m_tableLayout, dataCellsLayout, status );

  if( status.m_isMasterRank )
  {
    outputTableBottom( tableOutput, m_tableLayout, status.m_sepLine, !dataCellsLayout.empty() );
    tableOutput.flush();
  }
}

void TableTextMpiOutput::stretchColumnsByRanks( stdVector< size_t > & columnsWidth,
                                                TableTextMpiOutputStatus const status ) const
{
  size_t const rankColumnsCount = columnsWidth.size();
  size_t const maxRanksColumnsCount = MpiWrapper::max( rankColumnsCount );

  if( status.m_isContributing )
    GEOS_ASSERT_EQ( rankColumnsCount, maxRanksColumnsCount ); // TODO: contribute to the new table error system with this one

  columnsWidth.resize( maxRanksColumnsCount, 0 );

  MpiWrapper::allReduce( columnsWidth, columnsWidth, MpiWrapper::Reduction::Max );
}

void TableTextMpiOutput::outputTableDataToRank0( std::ostream & tableOutput,
                                                 PreparedTableLayout const & tableLayout,
                                                 CellLayoutRows const & dataCellsLayout,
                                                 TableTextMpiOutputStatus const status ) const
{
  integer const ranksCount = MpiWrapper::commSize();

  // master rank does the output directly to the output, other ranks will have to send it through a string.
  std::ostringstream localStringStream;
  std::ostream & rankOutput = status.m_isMasterRank ? tableOutput : localStringStream;

  if( status.m_isContributing )
  {
    outputTableData( rankOutput, tableLayout, dataCellsLayout );
  }

  // all other ranks than rank 0 render their output in a string and comunicate its size
  stdVector< integer > ranksStrsSizes{ ranksCount, 0 };
  string const rankStr = !status.m_isMasterRank && status.m_isContributing ? localStringStream.str() : "";
  integer const rankStrSize = rankStr.size();

  localStringStream.clear();
  MpiWrapper::allgather( &rankStrSize, 1, ranksStrsSizes.data(), 1 );

  // we compute the memory layout of the ranks strings
  stdVector< integer > ranksStrsDisps{ ranksCount, 0 };
  integer ranksStrsTotalSize = 0;
  for( integer rankId = 1; rankId < ranksCount; ++rankId )
  {
    ranksStrsDisps[rankId] = ranksStrsTotalSize;
    ranksStrsTotalSize += ranksStrsSizes[rankId];
  }

  // finally, we can send all text data to rank 0, then we output it in the output stream.
  string ranksStrs = string( ranksStrsTotalSize, '\0' );
  MpiWrapper::gatherv( &rankStr[0], rankStrSize,
                       &ranksStrs[0], ranksStrsSizes.data(), ranksStrsDisps.data(),
                       0, MPI_COMM_GEOS );
  if( status.m_isMasterRank )
  {
    for( integer rankId = 1; rankId < ranksCount; ++rankId )
    {
      if( m_mpiLayout.m_separatorBetweenRanks )
        tableOutput << tableLayout.getIndentationStr() << m_verticalLine
                    << string( status.m_sepLine.size() - 2, m_horizontalLine )
                    << m_verticalLine << '\n';

      if( ranksStrsSizes[rankId] > 0 )
      {
        tableOutput << string_view( &ranksStrs[ranksStrsDisps[rankId]], ranksStrsSizes[rankId] );
      }
    }
  }
}

} /* namespace geos */
