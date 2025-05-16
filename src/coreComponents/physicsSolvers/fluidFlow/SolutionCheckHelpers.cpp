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

IdReporterCollector IdReporterBuffer::createCollector( arrayView1d< globalIndex > const & localToGlobalId ) const
{
  return IdReporterCollector( m_idsCounter, m_idsBuffer, localToGlobalId );
}

IdReporterOutput IdReporterBuffer::createOutput() const
{
  return IdReporterOutput( *this );
}


IdReporterOutput::IdReporterOutput( IdReporterBuffer const & buffer ):
  m_buffer( buffer )
{}

void IdReporterOutput::outputWrongValues( string_view linesPrefix,
                                          string_view valueNaming,
                                          real64 minValue,
                                          units::Unit unit ) const
{
  if( m_buffer.enabled() )
  {
    integer numNegativeValues = MpiWrapper::sum( m_buffer.getSignaledIdsCount() );
    if( numNegativeValues > 0 )
    {
      string const minValueStr = GEOS_FMT( "{:.{}f} [{}]", minValue, 3, units::getSymbol( unit ) );
      GEOS_LOG_RANK_0( GEOS_FMT( "{}{} {} values encountered. Minimum value: {}.",
                                 linesPrefix, numNegativeValues, valueNaming, minValueStr ) );
      GEOS_LOG_RANK_0( GEOS_FMT( "{}{} element ids:",
                                 linesPrefix, valueNaming ) );

      MpiWrapper::barrier();
      if( m_buffer.getCollectedIdsCount() > 0 )
      {
        GEOS_LOG( GEOS_FMT( "{}- rank {} ({} values): {} {}",
                            string( ' ', linesPrefix.size() ),
                            MpiWrapper::commRank(),
                            m_buffer.getSignaledIdsCount(),
                            stringutilities::join( m_buffer, ", " ),
                            ( m_buffer.isComplete() ? "..." : "." ) ) );
      }
    }
  }
}

} // namespace geos
