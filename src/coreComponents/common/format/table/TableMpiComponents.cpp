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
#include <functional>
#include <mpi.h>

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
  CellLayoutRows dataCellsLayout;
  CellLayoutRows errorCellsLayout;
  size_t tableTotalWidth = 0;

  {
    ColumnWidthModifier const columnWidthModifier = [this, status]( stdVector< size_t > & columnsWidth ) {
      stretchColumnsByRanks( columnsWidth, status );
    };
    initalizeTableGrids( m_tableLayout, tableData,
                         headerCellsLayout, dataCellsLayout, errorCellsLayout,
                         tableTotalWidth, columnWidthModifier );
    status.m_sepLine = string( tableTotalWidth, m_horizontalLine );
  }

  if( status.m_isMasterRank )
  {
    outputTableHeader( tableOutput, m_tableLayout, headerCellsLayout, status.m_sepLine );
    tableOutput.flush();
  }

  outputTableDataToRank0( tableOutput, m_tableLayout, dataCellsLayout, status );

  if( status.m_isMasterRank )
  {
    outputTableFooter( tableOutput, m_tableLayout, errorCellsLayout,
                       status.m_sepLine, status.m_hasContent );
    tableOutput.flush();
  }
}

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

void TableTextMpiOutput::outputTableDataToRank0( std::ostream & tableOutput,
                                                 PreparedTableLayout const & tableLayout,
                                                 CellLayoutRows const & dataCellsLayout,
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
    outputTableData( rankOutput, tableLayout, dataCellsLayout );
  }
  string const rankStr = !status.m_isMasterRank && status.m_isContributing ? localStringStream.str() : "";
  MpiWrapper::gatherString( rankStr, std::function< void(string_view) >( [&]( string_view str ){
    status.m_hasContent = true;
    tableOutput << str;
  } ));
}

} /* namespace geos */
