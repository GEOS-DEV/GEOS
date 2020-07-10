#include "TimeHistoryOutput.hpp"
namespace geosx
{
TimeHistoryOutput::TimeHistoryOutput( string const & name,
                                      Group * const parent ):
  OutputBase( name, parent ),
  m_collector_paths( ),
  m_format( ),
  m_filename( ),
  m_record_count( 0 ),
  m_io( )
{
  registerWrapper( viewKeys::timeHistoryOutputTarget, &m_collector_paths )->
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

  registerWrapper( viewKeys::timeHistoryRestart, &m_record_count )->
    setApplyDefaultValue( 0 )->
    setInputFlag( InputFlags::FALSE )->
    setRestartFlags( RestartFlags::WRITE_AND_READ )->
    setDescription( "The current history record to be written, on restart from an earlier time allows use to remove invalid future history." );

}

void TimeHistoryOutput::InitializePostSubGroups( Group * const group )
{
  {
    // check whether to truncate or append to the file up front so we don't have to bother during later accesses
    HDFFile( m_filename, (m_record_count == 0), false, MPI_COMM_GEOSX );
  }
  for( auto collector_path : m_collector_paths )
  {
    Group * tmp = this->GetGroupByPath( collector_path );
    HistoryCollection * collector = Group::group_cast< HistoryCollection * >( tmp );
    GEOSX_ERROR_IF( collector == nullptr, "The target of a time history output event must be a collector! " << collector_path );
    // todo: switch based on m_format, always hdf for now
    HistoryMetadata metadata = collector->GetMetadata( group );
    int rank = MpiWrapper::Comm_rank();
    if( rank == 0 )
    {
      HistoryMetadata time_metadata = collector->GetTimeMetadata( );
      time_metadata.setName( metadata.getName() + string( " " ) + time_metadata.getName());
      m_io.emplace_back( std::make_pair( std::make_unique< HDFSerialHistIO >( m_filename, metadata, m_record_count ),
                                         std::make_unique< HDFSerialHistIO >( m_filename, time_metadata, m_record_count, 4, MPI_COMM_SELF ) ) );
      collector->RegisterTimeBufferCall( [this]() { return this->m_io.back().second->GetBufferHead( ); } );
      m_io.back().second->Init( ( m_record_count > 0 ) );
    }
    else
    {
      m_io.emplace_back( std::make_pair( std::make_unique< HDFSerialHistIO >( m_filename, metadata, m_record_count ), std::unique_ptr< HDFSerialHistIO >( nullptr ) ) );
    }
    collector->RegisterBufferCall( [this]() { return this->m_io.back().first->GetBufferHead( ); } );
    m_io.back().first->Init( ( m_record_count > 0 ) );
    MpiWrapper::Barrier( MPI_COMM_GEOSX );
    if( m_record_count == 0 )
    {
      // do any 1-time metadata output
      globalIndex global_rank_offset = m_io.back().first->GetRankOffset( );
      localIndex meta_collector_count = collector->GetNumMetaCollectors( );
      Group * domain_group = dynamicCast< Group * >( dynamicCast< ProblemManager * >( group )->getDomainPartition( ) );
      for( localIndex meta_idx = 0; meta_idx < meta_collector_count; ++meta_idx )
      {
        std::unique_ptr< HistoryCollection > meta_collector = collector->GetMetaCollector( group, meta_idx, global_rank_offset );
        HistoryMetadata meta_metadata = meta_collector->GetMetadata( group );
        meta_metadata.setName( metadata.getName() + " " + meta_metadata.getName());
        std::unique_ptr< HDFSerialHistIO > meta_io = std::make_unique< HDFSerialHistIO >( m_filename, meta_metadata, 0, 1 );
        meta_collector->RegisterBufferCall( [&meta_io] () { return meta_io->GetBufferHead( ); } );
        meta_io->Init( false );
        meta_collector->Execute( 0.0, 0.0, 0, 0, 0, domain_group );
        meta_io->Write( );
      }
    }
    MpiWrapper::Barrier( MPI_COMM_GEOSX );
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
    m_record_count += io_pair.first->GetBufferedCount( );
    io_pair.first->Write( );
    if( io_pair.second )
    {
      io_pair.second->Write( );
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
    io_pair.first->CompressInFile();
    if( io_pair.second )
    {
      io_pair.second->CompressInFile();
    }
  }
}

REGISTER_CATALOG_ENTRY( OutputBase, TimeHistoryOutput, std::string const &, Group * const )
}
