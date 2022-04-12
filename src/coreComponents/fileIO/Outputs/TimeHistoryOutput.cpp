/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "TimeHistoryOutput.hpp"

#include "fileIO/timeHistory/TestingBuffer.hpp"

#if defined(GEOSX_USE_PYGEOSX)
#include "fileIO/python/PyHistoryOutputType.hpp"
#endif

namespace geosx
{
TimeHistoryOutput::TimeHistoryOutput( string const & name,
                                      Group * const parent ):
  OutputBase( name, parent ),
  m_collectorPaths( ),
  m_format( ),
  m_filename( ),
  m_hasRecordStarted( false ),
  m_io( )
{
  registerWrapper( viewKeys::timeHistoryOutputTargetString(), &m_collectorPaths ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "A list of collectors from which to collect and output time history information." );

  registerWrapper( viewKeys::timeHistoryOutputFilenameString(), &m_filename ).
    setApplyDefaultValue( "TimeHistory" ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The filename to which to write time history output." );

  registerWrapper( viewKeys::timeHistoryOutputFormatString(), &m_format ).
    setApplyDefaultValue( "hdf" ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The output file format for time history output." );

  registerWrapper( viewKeys::timeHistoryRestartString(), &m_hasRecordStarted ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "The current history record to be written, on restart from an earlier time allows use to remove invalid future history." );

}

void TimeHistoryOutput::initCollectorParallel( DomainPartition const & domain, HistoryCollection & collector )
{
  GEOSX_ASSERT( m_io.empty() );

  string const outputDirectory = getOutputDirectory();
  string const outputFile = joinPath( outputDirectory, m_filename );

  // TODO Is HistoryCollection a Collection or a Collector?
  auto registerBufferCalls = [&]( HistoryCollection & hc, string prefix = "" )
  {
    for( localIndex collectorIdx = 0; collectorIdx < hc.numCollectors(); ++collectorIdx )
    {
      HistoryMetadata metadata = hc.getMetaData( domain, collectorIdx );

      // TODO Deal with this in PackCollection::buildMetaDataCollectors?
      if( !prefix.empty() )
      { metadata.setName( prefix + metadata.getName() ); }

      m_io.emplace_back( std::make_unique< HDFHistIO >( outputFile, metadata ) );
//      m_io.emplace_back( std::make_unique< TestingBuffer >( metadata.getName(), metadata.getDims() ) );
      hc.registerBufferProvider( collectorIdx, [this, idx = m_io.size() - 1]( localIndex count )
      {
        m_io[idx]->updateCollectingCount( count );
        return m_io[idx]->getBufferHead();
      } );
      m_io.back()->init( m_hasRecordStarted );
    }
  };

  // FIXME Why stop (pseudo) recursion at one single level?
  registerBufferCalls( collector );

  for( localIndex metaIdx = 0; metaIdx < collector.numMetaDataCollectors(); ++metaIdx )
  {
    registerBufferCalls( collector.getMetaDataCollector( metaIdx ), collector.getTargetName() + " " );
  }

  // do the time output last so its at the end of the m_io list, since writes are parallel we need
  //  the rest of the collectors to share position in the list across the world comm

  // rank == 0 does time output for the collector
  if( MpiWrapper::commRank() == 0 )
  {
    HistoryMetadata timeMetadata = collector.getTimeMetaData();
    m_io.emplace_back( std::make_unique< HDFHistIO >( outputFile, timeMetadata, MPI_COMM_SELF ) );
//    m_io.emplace_back( std::make_unique< TestingBuffer >( timeMetadata.getName(), timeMetadata.getDims() ) );
    // We copy the back `idx` not to rely on possible future appends to `m_io`.
    collector.registerTimeBufferProvider( [this, idx = m_io.size() - 1]()
                                          { return m_io[idx]->getBufferHead(); } );
    m_io.back()->init( m_hasRecordStarted );
  }

  MpiWrapper::barrier( MPI_COMM_GEOSX );
}

void TimeHistoryOutput::initializePostInitialConditionsPostSubGroups()
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
    HDFFile( outputFile, !m_hasRecordStarted, true, MPI_COMM_GEOSX );
  }

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  for( auto collectorPath : m_collectorPaths )
  {
    HistoryCollection & collector = this->getGroupByPath< HistoryCollection >( collectorPath );
    collector.initializePostSubGroups();
    initCollectorParallel( domain, collector );
  }
}

void TimeHistoryOutput::setFileName( string const & root )
{
  m_filename = root;
}

void TimeHistoryOutput::reinit()
{
  m_hasRecordStarted = false;
  m_io.clear();
  initializePostInitialConditionsPostSubGroups();
}

bool TimeHistoryOutput::execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                 real64 const GEOSX_UNUSED_PARAM( dt ),
                                 integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                 integer const GEOSX_UNUSED_PARAM( eventCounter ),
                                 real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                                 DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  GEOSX_MARK_FUNCTION;
  for( auto & th_io : m_io )
  {
    th_io->write( );
  }
  m_hasRecordStarted = true;
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

#if defined(GEOSX_USE_PYGEOSX)
PyTypeObject * TimeHistoryOutput::getPythonType() const
{ return python::getPyHistoryOutputType(); }
#endif

REGISTER_CATALOG_ENTRY( OutputBase, TimeHistoryOutput, string const &, Group * const )
}
