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
 * @file SolutionCheckKernel.cpp
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
                                     addColumns( { "Global Id", units::getDescription( unit ) } ).
                                     enableLineBreak( false ).
                                     setIndentation( linesPrefix.size() ).
                                     setDefaultHeaderAlignment( TableLayout::Alignment::left );
        TableData data;
        integer const signaledCount = m_buffer.getSignaledElementsCount();
        integer const collectedCount = m_buffer.getCollectedElementsCount();
        integer const omittedCount = signaledCount - collectedCount;

        TableMpiLayout mpiLayout;
        mpiLayout.m_separatorBetweenRanks = true;

        if( signaledCount > 0 )
        {
          mpiLayout.m_rankTitle = GEOS_FMT( "Rank {}, {} / {} values",
                                            MpiWrapper::commRank(), collectedCount, signaledCount );

          for( ElementReport const & report : m_buffer )
          {
            data.addRow( report.m_id, report.m_value );
          }

          // adding one last line for signaling partial data & readability
          if( omittedCount > 0 )
          {
            data.addRow( "...", "..." );
          }
        }

        TableTextMpiFormatter const formatter = TableTextMpiFormatter( layout, mpiLayout );
        formatter.toStream( std::cout, data );
        GEOS_LOG_RANK_0( '\n' );
      }
      else
      {
        GEOS_LOG_RANK_0( GEOS_FMT( "{}Increase the log level to enable a reporting of the {} values.",
                                   string( linesPrefix.size(), ' ' ), valueNaming ) );
      }
    }
  }
}

}   // namespace geos
