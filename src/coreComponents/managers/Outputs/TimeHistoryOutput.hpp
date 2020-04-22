#ifndef GEOSX_TIME_HISTORY_OUTPUT_HPP_
#define GEOSX_TIME_HISTORY_OUTPUT_HPP_

#include "OutputBase.hpp"
#include "TimeHistoryCollector.hpp"
#include "TimeHistory.hpp"
#include "managers/TimeHistory/HistoryIO.hpp"
#include "fileIO/hdf/HDFFile.hpp"

#include "cxx-utilities/src/Array.hpp" // just for collector

namespace geosx
{
  class TimeHistoryOutput : public OutputBase
  {
  public:
    TimeHistoryOutput( string const & name,
                       Group * const parent ):
      OutputBase(name,parent),
      m_collector_path( ),
      m_format( ),
      m_filename( ),
      m_record_count(false),
      m_time_collection( "TimeCollection", this ),
      m_io(nullptr),
      m_t_io(nullptr)
    {
      registerWrapper(viewKeys::timeHistoryOutputTarget, &m_collector_path, false)->
        setInputFlag(InputFlags::REQUIRED)->
        setDescription("A collector from which to collect and output time history information.");

      registerWrapper(viewKeys::timeHistoryOutputFilename, &m_filename, false)->
        setApplyDefaultValue("TimeHistory")->
        setInputFlag(InputFlags::OPTIONAL)->
        setDescription("The filename to which to write time history output.");

      registerWrapper(viewKeys::timeHistoryOutputFormat, &m_format, false)->
        setApplyDefaultValue("hdf")->
        setInputFlag(InputFlags::OPTIONAL)->
        setDescription("The output file format for time history output.");

      registerWrapper(viewKeys::timeHistoryRestart, &m_record_count, false)->
        setApplyDefaultValue(0)->
        setInputFlag(InputFlags::FALSE)->
        setRestartFlags(RestartFlags::WRITE_AND_READ)->
        setDescription("The current history record to be written, on restart from an earlier time allows use to remove invalid future history.");
    }

    virtual ~TimeHistoryOutput() override
    { }

    virtual void SetupDirectoryStructure() override
    {
      Group * tmp = this->GetGroupByPath(m_collector_path);
      HistoryCollection * collector = Group::group_cast<HistoryCollection*>(tmp);
      GEOSX_ERROR_IF(collector == nullptr, "The target of a time history output event must be a collector! " << m_collector_path);
      // switch based on m_format, always hdf for now
      m_io = std::make_unique<HDFHistIO>( m_filename, collector->GetMetadata( ), m_record_count );
      collector->RegisterBufferCall([this]() { return this->m_io->GetBufferHead( ); });
      m_io->Init( m_record_count > 0 );
      int rank = MpiWrapper::Comm_rank();
      if( rank == 0 )
      {
        m_t_io = std::make_unique<HDFHistIO>( m_filename, m_time_collection.GetMetadata( ), m_record_count, 4, MPI_COMM_SELF );
        m_time_collection.RegisterBufferCall([this]() { return this->m_t_io->GetBufferHead( ); });
        m_t_io->Init( m_record_count > 0 );
      }
    }

    /// This method will be called by the event manager if triggered
    virtual void Execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                          real64 const GEOSX_UNUSED_PARAM( dt ),
                          integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                          integer const GEOSX_UNUSED_PARAM( eventCounter ),
                          real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                          dataRepository::Group * GEOSX_UNUSED_PARAM( domain ) ) override
    {
      m_record_count += m_io->GetBufferedCount( );
      m_io->Write( );
      int rank = MpiWrapper::Comm_rank();
      if( rank == 0 )
      {
        m_t_io->Write( );
      }
    }

    /// Write one final output as the code exits
    virtual void Cleanup( real64 const time_n,
                          integer const cycleNumber,
                          integer const eventCounter,
                          real64 const eventProgress,
                          dataRepository::Group * domain ) override
    {
      Execute(time_n,0.0,cycleNumber,eventCounter,eventProgress,domain);
    }

    static string CatalogName() { return "TimeHistoryOutput"; }

    struct viewKeys
    {
      static constexpr auto timeHistoryOutputTarget = "source";
      static constexpr auto timeHistoryOutputFilename = "filename";
      static constexpr auto timeHistoryOutputFormat = "format";
      static constexpr auto timeHistoryWriteHead = "write_head";
      static constexpr auto timeHistoryRestart = "restart";
    } timeHistoryOutputViewKeys;


    private:
      string m_collector_path;
      string m_format;
      string m_filename;

      integer m_record_count;

      TimeCollection m_time_collection;

      std::unique_ptr<BufferedHistoryIO> m_io;
      std::unique_ptr<BufferedHistoryIO> m_t_io;
  };
}

#endif
