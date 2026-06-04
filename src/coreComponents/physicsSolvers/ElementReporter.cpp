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
 * @file ElementReporter.cpp
 */

#include "physicsSolvers/ElementReporter.hpp"
#include "common/MpiWrapper.hpp"
#include "common/format/StringUtilities.hpp"
#include "common/format/table/TableMpiComponents.hpp"

namespace geos
{

ElementReporterBuffer::ElementReporterBuffer( bool enabled, ElementCount maxCollectionSize ):
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
ElementReporterBuffer::createCollector( arrayView1d< globalIndex const > const & localToGlobalId ) const
{
  return ElementsReporterCollector( m_elementsCounter, m_elementsBuffer, localToGlobalId );
}

void ElementReporterBuffer::finalizeCollection()
{
  m_elementsCounter.move( LvArray::MemorySpace::host, false );
  m_elementsBuffer.move( LvArray::MemorySpace::host, false );
}

void ElementReporterBuffer::clear()
{
  if( enabled() )
    m_elementsCounter[0] = 0;
}


ElementReporterOutput::ElementReporterOutput():
  m_buffer( false, 0 )
{}

ElementReporterOutput & ElementReporterOutput::setValueMetadata( string_view valuesNaming,
                                                                 units::Unit valuesUnit )
{
  m_valuesNaming = valuesNaming;
  m_valuesUnit = valuesUnit;
  return *this;
}

ElementReporterOutput & ElementReporterOutput::setLogPrefix( string_view logPrefix )
{
  m_logPrefix = logPrefix;
  return *this;
}

ElementReporterOutput & ElementReporterOutput::setRanges( std::optional< real64 > minValue,
                                                          std::optional< real64 > maxValue )
{
  if( minValue )
    m_minEncounteredValue = *minValue;

  if( maxValue )
    m_maxEncounteredValue = *maxValue;

  return *this;
}

ElementReporterOutput & ElementReporterOutput::setValues( ElementReporterBuffer && buffer,
                                                          std::optional< real64 > rankMinValue,
                                                          std::optional< real64 > rankMaxValue )
{
  m_ranksSignaledElementsCount = MpiWrapper::sum( buffer.getSignaledElementsCount() );
  m_ranksCollectedElementsCount = MpiWrapper::sum( buffer.getCollectedElementsCount() );

  if( rankMinValue )
    m_minEncounteredValue = MpiWrapper::min( *rankMinValue );

  if( rankMaxValue )
    m_maxEncounteredValue = MpiWrapper::max( *rankMaxValue );

  buffer.finalizeCollection();
  m_buffer = std::move( buffer );

  return *this;
}

void ElementReporterOutput::outputReport( string_view eventTitle ) const
{
  if( m_buffer.enabled() )
  {
    if( m_ranksSignaledElementsCount > 0 )
    {
      { // report top description
        auto const unitStr = units::getSymbol( m_valuesUnit );
        std::ostringstream descOSS;
        descOSS << GEOS_FMT( "{}{} {} values encountered",
                             m_logPrefix, m_ranksSignaledElementsCount, m_valuesNaming );
        if( m_minEncounteredValue || m_maxEncounteredValue )
        {
          descOSS << ", ";
          descOSS << ( m_minEncounteredValue ? GEOS_FMT( " from {} [{}]", *m_minEncounteredValue, unitStr ) : "" );
          descOSS << ( m_maxEncounteredValue ? GEOS_FMT( " to {} [{}]", *m_maxEncounteredValue, unitStr ) : "" );
          if( m_minAcceptableValue || m_maxAcceptableValue )
          {
            descOSS << '\n' << string( m_logPrefix.size(), ' ' ) << "Acceptable range ";
            descOSS << ( m_minAcceptableValue ? GEOS_FMT( " from {} [{}]", *m_minAcceptableValue, unitStr ) : "" );
            descOSS << ( m_maxAcceptableValue ? GEOS_FMT( " to {} [{}]", *m_minAcceptableValue, unitStr ) : "" );
          }
        }
        string const rangesJoinedStr = descOSS.str();
        GEOS_LOG_RANK_0( GEOS_FMT( "{}{} {} values encountered",
                                   m_logPrefix, m_ranksSignaledElementsCount, m_valuesNaming ) );
      }

      if( m_ranksCollectedElementsCount > 0 )
      {
        TableLayout const layout = TableLayout().
                                     setTitle( GEOS_FMT( "{}\nSummary of {} elements", eventTitle, m_valuesNaming ) ).
                                     addColumns( { "Global Id", units::getDescription( m_valuesUnit ) } ).
                                     enableLineBreak( false ).
                                     setIndentation( m_logPrefix.size() ).
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

        TableTextMpiOutput const formatter = TableTextMpiOutput( layout, mpiLayout );
        formatter.toStream( std::cout, data );
        GEOS_LOG_RANK_0( '\n' );
      }
      else
      {
        GEOS_LOG_RANK_0( GEOS_FMT( "{}Increase the log level to enable a reporting of the {} values.",
                                   string( m_logPrefix.size(), ' ' ), m_valuesNaming ) );
      }
    }
  }
}

void ElementReporterOutput::clear()
{
  m_buffer.clear();

  m_ranksSignaledElementsCount = 0;
  m_ranksCollectedElementsCount = 0;
  m_minEncounteredValue.reset();
  m_maxEncounteredValue.reset();
}

}   // namespace geos
