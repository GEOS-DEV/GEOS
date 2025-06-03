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
#include "common/format/table/TableMpiComponents.hpp"

namespace geos
{

ElementsReporterBuffer::ElementsReporterBuffer( bool enabled, ElementCount maxCollectionSize ):
  m_elementsCounter( enabled ? 1 : 0 ),
  m_elementsBuffer( enabled ? maxCollectionSize : 0 )
{
  if( enabled )
  {
    m_elementsCounter.zero();
    m_elementsBuffer.zero();
  }
}

ElementsReporterCollector
ElementsReporterBuffer::createCollector( arrayView1d< globalIndex const > const & localToGlobalId ) const
{
  return ElementsReporterCollector( m_elementsCounter, m_elementsBuffer, localToGlobalId );
}

ElementsReporterOutput ElementsReporterBuffer::createOutput() const
{
  m_elementsCounter.move( LvArray::MemorySpace::host, false );
  m_elementsBuffer.move( LvArray::MemorySpace::host, false );
  return ElementsReporterOutput( *this );
}


ElementsReporterOutput::ElementsReporterOutput( ElementsReporterBuffer const & buffer ):
  m_buffer( buffer ),
  m_ranksSignaledElementsCount( MpiWrapper::sum( buffer.getSignaledElementsCount() ) ),
  m_ranksCollectedElementsCount( MpiWrapper::sum( buffer.getCollectedElementsCount() ) )
{}

void ElementsReporterOutput::outputTooLowValues( string_view linesPrefix,
                                                 string_view valueNaming,
                                                 real64 minValue,
                                                 units::Unit unit ) const
{
  if( m_buffer.enabled() )
  {
    if( m_ranksSignaledElementsCount > 0 )
    {

      string const minValueStr = GEOS_FMT( "{:.{}f} [{}]", minValue, 3, units::getSymbol( unit ) );
      GEOS_LOG_RANK_0( GEOS_FMT( "{}{} {} values encountered. Minimum value: {}.",
                                 linesPrefix, m_ranksSignaledElementsCount, valueNaming, minValueStr ) );

      if( m_ranksCollectedElementsCount > 0 )
      {
        TableLayout const layout = TableLayout().
                                     setTitle( GEOS_FMT( "Summary of {} elements", valueNaming ) ).
                                     enableLineBreak( false ).
                                     setIndentation( linesPrefix.size() ).
                                     setMargin( TableLayout::MarginValue::small ).
                                     setDefaultHeaderAlignment( TableLayout::Alignment::left );
        TableData data;
        integer const signaledCount = m_buffer.getSignaledElementsCount();
        integer const collectedCount = m_buffer.getCollectedElementsCount();
        integer const omittedCount = signaledCount - collectedCount;
        integer const tableColumnsCount = MpiWrapper::max( 1 + collectedCount + integer( omittedCount > 0 ) );

        if( signaledCount > 0 )
        {
          // adding a columns for row name, each collected value, and one last if a "..." have to be added
          auto & cells = data.getCellsData();
          string const title = GEOS_FMT( "Rank {}, {} / {} {} values:",
                                         MpiWrapper::commRank(),
                                         collectedCount,
                                         signaledCount,
                                         valueNaming );
          static constexpr integer titleLine = 0;
          static constexpr integer globalIdLine = 2;
          static constexpr integer valueLine = 3;

          data.addRow( stdVector< TableData::CellData >( tableColumnsCount, { CellType::MergeNext, "" } ) );
          cells[titleLine].back() = { CellType::Header, title };

          data.addRow( stdVector< TableData::CellData >( tableColumnsCount, { CellType::MergeNext, "" } ) );

          data.addRow( stdVector< TableData::CellData >( tableColumnsCount, { CellType::MergeNext, "" } ) );
          cells[globalIdLine].front() = { CellType::Value, "Global Id" };
          for( int i = 0; i < collectedCount; ++i )
            cells[globalIdLine][i+1] = { CellType::Value, std::to_string( m_buffer[i].m_id ) };

          data.addRow( stdVector< TableData::CellData >( tableColumnsCount, { CellType::MergeNext, "" } ) );
          cells[valueLine].front() = { CellType::Value, string( units::getDescription( unit ) ) };
          for( int i = 0; i < collectedCount; ++i )
            cells[valueLine][i+1] = { CellType::Value, std::to_string( m_buffer[i].m_value ) };

          if( omittedCount > 0 )
          {
            cells[globalIdLine].back() = { CellType::Value, "..." };
            cells[valueLine].back() = { CellType::Value, "..." };
          }
        }

        TableTextMpiOutput const formatter = TableTextMpiOutput( layout, TableMpiLayout{ true } );
        formatter.toStream( std::cout, data );
        GEOS_LOG_RANK_0( '\n' );
      }
      else
      {
        GEOS_LOG( GEOS_FMT( "{}Increase the log level to enable a reporting of the {} values.",
                            string( linesPrefix.size(), ' ' ), valueNaming ) );
      }
    }
  }
}

}   // namespace geos
