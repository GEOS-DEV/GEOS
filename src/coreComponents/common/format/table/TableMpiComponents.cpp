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

TableTextMpiOutput::TableTextMpiOutput( ParallelOutputMode parallelOutputMode ):
  TableTextFormatter(),
  m_parallelOutputMode( parallelOutputMode )
{}

TableTextMpiOutput::TableTextMpiOutput( TableLayout const & tableLayout,
                                        ParallelOutputMode parallelOutputMode ):
  TableTextFormatter( tableLayout ),
  m_parallelOutputMode( parallelOutputMode )
{}

template<>
void TableTextMpiOutput::toStream< TableData >( std::ostream & tableOutput,
                                                TableData const & tableData ) const
{
  // only the master rank does the output of the header && bottom of the table
  bool isMasterRank = MpiWrapper::commRank() == 0;
  // some ranks does not have any output to produce, they can sleep while waiting the next barrier.
  bool isContributing = !tableData.getCellsData().empty() || isMasterRank;

  CellLayoutRows headerCellsLayout;
  CellLayoutRows dataCellsLayout;
  size_t tableTotalWidth = 0;
  string sepLine;

  MPI_Barrier( MPI_COMM_GEOS );

  {
    ColumnWidthModifier const columnWidthModifier = [this]( stdVector< size_t > & columnsWidth ) {
      stretchColumnsByRanks( columnsWidth );
    };
    initalizeTableGrids( m_tableLayout, tableData,
                         headerCellsLayout, dataCellsLayout,
                         tableTotalWidth, columnWidthModifier );
  }
  MPI_Barrier( MPI_COMM_GEOS );

  if( isMasterRank )
  {
    sepLine = string( tableTotalWidth, m_horizontalLine );
    outputTableHeader( tableOutput, m_tableLayout, headerCellsLayout, sepLine );
    tableOutput.flush();
  }
  MPI_Barrier( MPI_COMM_GEOS );

  if( isContributing )
  {
    outputTableData( tableOutput, m_tableLayout, dataCellsLayout,
                     m_parallelOutputMode == ParallelOutputMode::MixedRanksRows );

    if( m_parallelOutputMode == ParallelOutputMode::InsecableRanks )
      tableOutput.flush();
  }
  MPI_Barrier( MPI_COMM_GEOS );

  if( isMasterRank )
  {
    outputTableBottom( tableOutput, m_tableLayout, sepLine, !dataCellsLayout.empty() );
    tableOutput.flush();
  }

}

void TableTextMpiOutput::stretchColumnsByRanks( stdVector< size_t > & columnsWidth ) const
{
  stdVector< size_t > oldColumnsWidth( columnsWidth );
  MpiWrapper::allReduce( oldColumnsWidth, columnsWidth, int( columnsWidth.size() ),
                         MpiWrapper::Reduction::Max, MPI_COMM_GEOS );
}

} /* namespace geos */
