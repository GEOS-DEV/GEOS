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
 * @file SolutionCheckKernel.hpp
 */

#include "physicsSolvers/fluidFlow/SolutionCheckHelpers.hpp"
#include "common/MpiWrapper.hpp"
#include "common/format/StringUtilities.hpp"
#include "common/format/table/TableFormatter.hpp"

namespace geos
{

IdReporterBuffer::IdReporterBuffer( bool enabled, IdCountType maxCollectionSize ):
  m_idsCounter( enabled ? 1 : 0 ),
  m_idsBuffer( enabled ? maxCollectionSize : 0 )
{
  if( enabled )
  {
    m_idsCounter.zero();
    m_idsBuffer.zero();
  }
}

IdReporterCollector IdReporterBuffer::createCollector( arrayView1d< globalIndex const > const & localToGlobalId ) const
{
  return IdReporterCollector( m_idsCounter, m_idsBuffer, localToGlobalId );
}

IdReporterOutput IdReporterBuffer::createOutput() const
{
  m_idsCounter.move( LvArray::MemorySpace::host, false );
  m_idsBuffer.move( LvArray::MemorySpace::host, false );
  return IdReporterOutput( *this );
}


IdReporterOutput::IdReporterOutput( IdReporterBuffer const & buffer ):
  m_buffer( buffer ),
  m_ranksSignaledIdsCount( MpiWrapper::sum( buffer.getSignaledIdsCount() ) ),
  m_ranksCollectedIdsCount( MpiWrapper::sum( buffer.getCollectedIdsCount() ) )
{}

void IdReporterOutput::outputWrongValues( string_view linesPrefix,
                                          string_view valueNaming,
                                          real64 minValue,
                                          units::Unit unit ) const
{
  if( m_buffer.enabled() )
  {
    if( m_ranksSignaledIdsCount > 0 )
    {
      string const minValueStr = GEOS_FMT( "{:.{}f} [{}]", minValue, 3, units::getSymbol( unit ) );
      GEOS_LOG_RANK_0( GEOS_FMT( "{}{} {} values encountered. Minimum value: {}.",
                                 linesPrefix, m_ranksSignaledIdsCount, valueNaming, minValueStr ) );

      if( m_ranksCollectedIdsCount > 0 )
      {
        TableLayout const layout = TableLayout().
                                     setTitle( GEOS_FMT( "Summary of {} elements", valueNaming ) ).
                                     enableLineBreak( false ).
                                     setIndentation( linesPrefix.size() ).
                                     setMargin( TableLayout::MarginValue::small ).
                                     setDefaultHeaderAlignment( TableLayout::Alignment::left );
        TableData data;
        integer const signaledCount = m_buffer.getSignaledIdsCount();
        integer const collectedCount = m_buffer.getCollectedIdsCount();

        if( signaledCount > 0 )
        {
          integer const omittedCount = signaledCount - collectedCount;
          // adding a columns for row name, each collected value, and one last if a "..." have to be added
          integer const columnsCount = MpiWrapper::max( 1 + collectedCount + integer( omittedCount > 0 ) );
          enum class Lines : integer { Title, Separator, GlobalId, Value };
          auto & cells = data.getCellsData();
          string const title = GEOS_FMT( "Rank {}, {} / {} {} values:",
                                         MpiWrapper::commRank(), collectedCount, signaledCount, valueNaming );

          data.addRow( stdVector< TableData::CellData >( columnsCount, { CellType::MergeNext, "" } ) );
          cells[integer( Lines::Title )].back() = { CellType::Header, title };

          data.addRow( stdVector< TableData::CellData >( columnsCount, { CellType::MergeNext, "" } ) );

          data.addRow( stdVector< TableData::CellData >( columnsCount, { CellType::MergeNext, "" } ) );
          cells[integer( Lines::GlobalId )].front() = { CellType::Value, "Global Id" };
          for( int i = 0; i < collectedCount; ++i )
            cells[integer( Lines::GlobalId )][i+1] = { CellType::Value, std::to_string( m_buffer[i] ) };

          data.addRow( stdVector< TableData::CellData >( columnsCount, { CellType::MergeNext, "" } ) );
          cells[integer( Lines::Value )].front() = { CellType::Value, string( units::getDescription( unit ) ) };
          for( int i = 0; i < collectedCount; ++i )
            cells[integer( Lines::Value )][i+1] = { CellType::Value, std::to_string( 0.0 ) };

          if( omittedCount > 0 )
          {
            cells[integer( Lines::GlobalId )].back() = { CellType::Value, "..." };
            cells[integer( Lines::Value )].back() = { CellType::Value, "..." };
          }
        }

        TableTextFormatter const formatter( layout );
        GEOS_LOG( formatter.toString( data ) );
      }
    }
  }
}

}   // namespace geos
