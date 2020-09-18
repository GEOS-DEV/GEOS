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
  registerWrapper( viewKeys::timeHistoryOutputTarget, &m_collectorPaths )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "A list of collectors from which to collect and output time history information." );

  registerWrapper( viewKeys::timeHistoryOutputFilename, &m_filename )->
    setApplyDefaultValue( "TimeHistory" )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "The filename to which to write time history output." );

  registerWrapper( viewKeys::timeHistoryOutputFormat, &m_format )->
    setApplyDefaultValue( "hdf" )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "The output file format for time history output." );

  registerWrapper( viewKeys::timeHistoryRestart, &m_recordCount )->
    setApplyDefaultValue( 0 )->
    setInputFlag( InputFlags::FALSE )->
    setRestartFlags( RestartFlags::WRITE_AND_READ )->
    setDescription( "The current history record to be written, on restart from an earlier time allows use to remove invalid future history." );

}

void TimeHistoryOutput::initCollectorParallel( ProblemManager & pm, HistoryCollection * collector )
{
  bool freshInit = ( m_recordCount == 0 );
  // rank == 0 do time output for the collector
  for( localIndex ii = 0; ii < collector->getCollectionCount( ); ++ii )
  {
    HistoryMetadata metadata = collector->getMetadata( pm, ii );
    m_io.emplace_back( std::make_unique< HDFHistIO >( m_filename, metadata, m_recordCount ) );
    collector->registerBufferCall( ii, [this, ii]() { return m_io[ii]->getBufferHead( ); } );
    m_io.back()->init( !freshInit );
  }
  int rnk = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  if( rnk == 0 )
  {
    HistoryMetadata timeMetadata = collector->getTimeMetadata( );
    m_io.emplace_back( std::make_unique< HDFHistIO >( m_filename, timeMetadata, m_recordCount, 2, 2, MPI_COMM_SELF ) );
    collector->registerTimeBufferCall( [this]() { return m_io.back()->getBufferHead( ); } );
    m_io.back()->init( !freshInit );
  }
  if( freshInit )
  {
    // do any 1-time metadata output
    localIndex metaCollectorCount = collector->getNumMetaCollectors( );
    Group * domainGroup = dynamicCast< Group * >( pm.getDomainPartition( ) );
    for( localIndex metaIdx = 0; metaIdx < metaCollectorCount; ++metaIdx )
    {
      std::unique_ptr< HistoryCollection > metaCollector = collector->getMetaCollector( pm, metaIdx );
      std::vector< std::unique_ptr< HDFHistIO > > metaIOs( metaCollector->getCollectionCount( ) );
      for( localIndex ii = 0; ii < metaCollector->getCollectionCount( ); ++ii )
      {
        HistoryMetadata metaMetadata = metaCollector->getMetadata( pm, ii );
        metaMetadata.setName( collector->getTargetName() + " " + metaMetadata.getName( ) );
        metaIOs[ii] = std::make_unique< HDFHistIO >( m_filename, metaMetadata, 0, 1 );
        metaCollector->registerBufferCall( ii, [&metaIOs, ii] () { return metaIOs[ii]->getBufferHead( ); } );
        metaIOs[ii]->init( false );
      }
      metaCollector->Execute( 0.0, 0.0, 0, 0, 0, domainGroup );
      for( localIndex ii = 0; ii < metaCollector->getCollectionCount( ); ++ii )
      {
        metaIOs[ii]->write( );
      }
    }
  }
  MpiWrapper::Barrier( MPI_COMM_GEOSX );
}

void TimeHistoryOutput::InitializePostSubGroups( Group * const group )
{
  {
    // check whether to truncate or append to the file up front so we don't have to bother during later accesses
    HDFFile( m_filename, (m_recordCount == 0), true, MPI_COMM_GEOSX );
  }
  ProblemManager & pm = dynamicCast< ProblemManager & >( *group );
  for( auto collector_path : m_collectorPaths )
  {
    Group * tmp = this->GetGroupByPath( collector_path );
    HistoryCollection * collector = Group::group_cast< HistoryCollection * >( tmp );
    GEOSX_ERROR_IF( collector == nullptr, "The target of a time history output event must be a collector! " << collector_path );
    collector->InitializePostSubGroups( group );
    initCollectorParallel( pm, collector );
  }
}

void TimeHistoryOutput::Execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                 real64 const GEOSX_UNUSED_PARAM( dt ),
                                 integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                 integer const GEOSX_UNUSED_PARAM( eventCounter ),
                                 real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                                 dataRepository::Group * GEOSX_UNUSED_PARAM( domain ) )
{
  GEOSX_MARK_FUNCTION;
  localIndex newBuffered = m_io.front()->getBufferedCount( );
  for( auto & th_io : m_io )
  {
    GEOSX_ERROR_IF( newBuffered != th_io->getBufferedCount( ), "Inconsistent buffered time history count from single collector." );
    th_io->write( );
  }
  m_recordCount += newBuffered;
}

void TimeHistoryOutput::Cleanup( real64 const time_n,
                                 integer const cycleNumber,
                                 integer const eventCounter,
                                 real64 const eventProgress,
                                 dataRepository::Group * domain )
{
  Execute( time_n, 0.0, cycleNumber, eventCounter, eventProgress, domain );
  // remove any unused trailing space reserved to write additional histories
  for( auto & th_io : m_io )
  {
    th_io->compressInFile();
  }
}

REGISTER_CATALOG_ENTRY( OutputBase, TimeHistoryOutput, std::string const &, Group * const )
}
