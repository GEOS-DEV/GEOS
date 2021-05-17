/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "TimeHistoryOutput.hpp"

namespace geosx
{
TimeHistoryOutput::TimeHistoryOutput( string const & name,
                                      Group * const parent ):
  OutputBase( name, parent ),
  m_collectorPaths( ),
  m_format( ),
  m_filename( ),
  m_recordCount( 0 ),
  m_io( )
{
  registerWrapper( viewKeys::timeHistoryOutputTarget, &m_collectorPaths ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "A list of collectors from which to collect and output time history information." );

  registerWrapper( viewKeys::timeHistoryOutputFilename, &m_filename ).
    setApplyDefaultValue( "TimeHistory" ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The filename to which to write time history output." );

  registerWrapper( viewKeys::timeHistoryOutputFormat, &m_format ).
    setApplyDefaultValue( "hdf" ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The output file format for time history output." );

  registerWrapper( viewKeys::timeHistoryRestart, &m_recordCount ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "The current history record to be written, on restart from an earlier time allows use to remove invalid future history." );

}

void TimeHistoryOutput::initCollectorParallel( DomainPartition & domain, HistoryCollection & collector )
{
  bool const freshInit = ( m_recordCount == 0 );

  string const outputDirectory = getOutputDirectory();
  string const outputFile = joinPath( outputDirectory, m_filename );

  // rank == 0 do time output for the collector
  for( localIndex ii = 0; ii < collector.getCollectionCount( ); ++ii )
  {
    HistoryMetadata metadata = collector.getMetadata( domain, ii );
    m_io.emplace_back( std::make_unique< HDFHistIO >( outputFile, metadata, m_recordCount ) );
    collector.registerBufferCall( ii, [this, ii, &pm, &collector]()
    { 
      HistoryMetadata metadata = collector.getMetadata( pm, ii );
      m_io[ii]->updateCollectingCount( metadata.getDims( )[0] );
      return m_io[ii]->getBufferHead( ); 
    } );
    m_io.back()->init( !freshInit );
  }

  // rank == 0 does time output for the collector
  if( MpiWrapper::commRank( MPI_COMM_GEOSX ) == 0 )
  {
    HistoryMetadata timeMetadata = collector.getTimeMetadata( );
    m_io.emplace_back( std::make_unique< HDFHistIO >( outputFile, timeMetadata, m_recordCount, 2, 2, MPI_COMM_SELF ) );
    collector.registerTimeBufferCall( [this]() { return m_io.back()->getBufferHead( ); } );
    m_io.back()->init( !freshInit );
  }
  if( freshInit )
  {
    GEOSX_MARK_SCOPE( "Fresh Time-hist initialization." );
    // do any 1-time metadata output
    localIndex metaCollectorCount = collector.getNumMetaCollectors( );
    for( localIndex metaIdx = 0; metaIdx < metaCollectorCount; ++metaIdx )
    {
      std::unique_ptr< HistoryCollection > metaCollector = collector.getMetaCollector( domain, metaIdx );
      std::vector< std::unique_ptr< HDFHistIO > > metaIOs( metaCollector->getCollectionCount( ) );
      for( localIndex ii = 0; ii < metaCollector->getCollectionCount( ); ++ii )
      {
        HistoryMetadata metaMetadata = metaCollector->getMetadata( domain, ii );
        metaMetadata.setName( collector.getTargetName() + " " + metaMetadata.getName( ) );
        metaIOs[ii] = std::make_unique< HDFHistIO >( outputFile, metaMetadata, 0, 1 );
        metaCollector->registerBufferCall( ii, [&metaIOs, ii] () { return metaIOs[ii]->getBufferHead( ); } );
        metaIOs[ii]->init( false );
      }

      metaCollector->execute( 0.0, 0.0, 0, 0, 0, domain );
      for( localIndex ii = 0; ii < metaCollector->getCollectionCount( ); ++ii )
      {
        metaIOs[ii]->write( );
      }
    }
  }
  MpiWrapper::barrier( MPI_COMM_GEOSX );
}

void TimeHistoryOutput::initializePostSubGroups()
{
  {
    // check whether to truncate or append to the file up front so we don't have to bother during later accesses
    string const outputDirectory = getOutputDirectory();
    if( MpiWrapper::commRank( MPI_COMM_GEOSX ) == 0 )
    {
      makeDirsForPath( outputDirectory );
    }
    MpiWrapper::barrier( MPI_COMM_GEOSX );
    string const outputFile = joinPath( outputDirectory, m_filename );
    HDFFile( outputFile, (m_recordCount == 0), true, MPI_COMM_GEOSX );
  }

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  for( auto collector_path : m_collectorPaths )
  {
    HistoryCollection & collector = this->getGroupByPath< HistoryCollection >( collector_path );
    collector.initializePostSubGroups();
    initCollectorParallel( domain, collector );
  }
}

bool TimeHistoryOutput::execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                 real64 const GEOSX_UNUSED_PARAM( dt ),
                                 integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                 integer const GEOSX_UNUSED_PARAM( eventCounter ),
                                 real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                                 DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  GEOSX_MARK_FUNCTION;
  localIndex newBuffered = m_io.front()->getBufferedCount( );
  for( auto & th_io : m_io )
  {
    GEOSX_ERROR_IF( newBuffered != th_io->getBufferedCount( ), "Inconsistent buffered time history count from single collector." );
    th_io->write( );
  }
  m_recordCount += newBuffered;

  return false;
}

void TimeHistoryOutput::cleanup( real64 const time_n,
                                 integer const cycleNumber,
                                 integer const eventCounter,
                                 real64 const eventProgress,
                                 DomainPartition & domain )
{
  execute( time_n, 0.0, cycleNumber, eventCounter, eventProgress, domain );
  // remove any unused trailing space reserved to write additional histories
  for( auto & th_io : m_io )
  {
    th_io->compressInFile();
  }
}

REGISTER_CATALOG_ENTRY( OutputBase, TimeHistoryOutput, string const &, Group * const )
}
