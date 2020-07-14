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

void TimeHistoryOutput::initCollectorSerial( ProblemManager & pm, HistoryCollection * collector )
{
  HistoryMetadata metadata = collector->getMetadata( pm );
  HistoryMetadata time_metadata = collector->getTimeMetadata( );
  time_metadata.setName( metadata.getName( ) + " " + time_metadata.getName( ) );
  m_io.emplace_back( std::make_pair( std::make_unique< HDFSerialHistIO >( m_filename, metadata, m_recordCount ),
                                     std::make_unique< HDFSerialHistIO >( m_filename, time_metadata, m_recordCount ) ) );
  collector->registerTimeBufferCall( [this]() { return this->m_io.back().second->getBufferHead( ); } );
  m_io.back().second->Init( ( m_recordCount > 0 ) );
  collector->registerBufferCall( [this]() { return this->m_io.back().first->getBufferHead( ); } );
  m_io.back().first->Init( ( m_recordCount > 0 ) );
  if( m_recordCount == 0 )
  {
    // do any 1-time metadata output
    localIndex meta_collector_count = collector->GetNumMetaCollectors( );
    Group * domain_group = dynamicCast< Group * >( dynamicCast< ProblemManager * >( group )->getDomainPartition( ) );
    for( localIndex meta_idx = 0; meta_idx < meta_collector_count; ++meta_idx )
    {
      std::unique_ptr< HistoryCollection > meta_collector = collector->getMetaCollector( group, meta_idx, 0 );
      HistoryMetadata meta_metadata = meta_collector->getMetadata( pm );
      meta_metadata.setName( metadata.getName() + " " + meta_metadata.getName());
      std::unique_ptr< HDFSerialHistIO > meta_io = std::make_unique< HDFSerialHistIO >( m_filename, meta_metadata, 0, 1 );
      meta_collector->registerBufferCall( [&meta_io] () { return meta_io->getBufferHead( ); } );
      meta_io->init( false );
      meta_collector->Execute( 0.0, 0.0, 0, 0, 0, domain_group );
      meta_io->write( );
    }
  }
  MpiWrapper::Barrier( MPI_COMM_GEOSX );
}

void TimeHistoryOutput::initCollectorParallel( ProblemManager & pm, HistoryCollection * collector )
{
  HistoryMetadata metadata = collector->getMetadata( pm );
  int rnk = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  if( rnk == 0 )
  {
    HistoryMetadata time_metadata = collector->getTimeMetadata( );
    time_metadata.setName( metadata.getName( ) + " " + time_metadata.getName( ) );
    m_io.emplace_back( std::make_pair( std::make_unique< HDFHistIO >( m_filename, metadata, m_recordCount ),
                                       std::make_unique< HDFHistIO >( m_filename, time_metadata, m_recordCount ) ) );
    collector->registerTimeBufferCall( [this]() { return this->m_io.back().second->getBufferHead( ); } );
    m_io.back().second->init( ( m_recordCount > 0 ) );
  }
  else
  {
    m_io.emplace_back( std::make_pair( std::make_unique< HDFHistIO >( m_filename, metadata, m_recordCount ),
                                       std::unique_ptr< HDFHistIO >( nullptr ) ) );
  }
  collector->registerBufferCall( [this]() { return this->m_io.back().first->getBufferHead( ); } );
  m_io.back().first->init( ( m_recordCount > 0 ) );
  if( m_recordCount == 0 )
  {
    // do any 1-time metadata output
    globalIndex global_rank_offset = m_io.back().first->getRankOffset( );
    localIndex meta_collector_count = collector->getNumMetaCollectors( );
    Group * domain_group = dynamicCast< Group * >( pm->getDomainPartition( ) );
    for( localIndex meta_idx = 0; meta_idx < meta_collector_count; ++meta_idx )
    {
      std::unique_ptr< HistoryCollection > meta_collector = collector->getMetaCollector( pm, meta_idx, global_rank_offset );
      HistoryMetadata meta_metadata = meta_collector->getMetadata( pm );
      meta_metadata.setName( metadata.getName() + " " + meta_metadata.getName());
      std::unique_ptr< HDFHistIO > meta_io = std::make_unique< HDFHistIO >( m_filename, meta_metadata, 0, 1 );
      meta_collector->RegisterBufferCall( [&meta_io] () { return meta_io->GetBufferHead( ); } );
      meta_io->init( false );
      meta_collector->Execute( 0.0, 0.0, 0, 0, 0, domain_group );
      meta_io->write( );
    }
  }
  MpiWrapper::Barrier( MPI_COMM_GEOSX );
}

void TimeHistoryOutput::InitializePostSubGroups( Group * const group )
{
  {
    // check whether to truncate or append to the file up front so we don't have to bother during later accesses
    HDFFile( m_filename, (m_recordCount == 0), false, MPI_COMM_GEOSX );
    //HDFFile( m_filename, (m_recordCount == 0), true, MPI_COMM_GEOSX );
  }
  ProblemManager & pm = dynamicCast< ProblemManager & >( *group );
  for( auto collector_path : m_collectorPaths )
  {
    Group * tmp = this->GetGroupByPath( collector_path );
    HistoryCollection * collector = Group::group_cast< HistoryCollection * >( tmp );
    GEOSX_ERROR_IF( collector == nullptr, "The target of a time history output event must be a collector! " << collector_path );
    initCollectorSerial( pm, collector );
    //initCollectorParallel( pm, collector );
  }
}

void TimeHistoryOutput::Execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                 real64 const GEOSX_UNUSED_PARAM( dt ),
                                 integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                 integer const GEOSX_UNUSED_PARAM( eventCounter ),
                                 real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                                 dataRepository::Group * GEOSX_UNUSED_PARAM( domain ) )
{
  for( auto & io_pair : m_io )
  {
    m_recordCount += io_pair.first->getBufferedCount( );
    io_pair.first->write( );
    if( io_pair.second )
    {
      io_pair.second->write( );
    }
  }
}

void TimeHistoryOutput::Cleanup( real64 const time_n,
                                 integer const cycleNumber,
                                 integer const eventCounter,
                                 real64 const eventProgress,
                                 dataRepository::Group * domain )
{
  Execute( time_n, 0.0, cycleNumber, eventCounter, eventProgress, domain );
  // remove any unused trailing space reserved to write additional histories
  for( auto & io_pair : m_io )
  {
    io_pair.first->compressInFile();
    if( io_pair.second )
    {
      io_pair.second->compressInFile();
    }
  }
}

REGISTER_CATALOG_ENTRY( OutputBase, TimeHistoryOutput, std::string const &, Group * const )
}
