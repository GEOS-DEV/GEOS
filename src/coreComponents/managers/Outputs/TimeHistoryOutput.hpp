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
      m_restart(false),
      m_io(nullptr),
      m_collector(nullptr)
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

      registerWrapper(viewKeys::timeHistoryRestart, &m_restart, false)->
        setApplyDefaultValue(false)->
        setInputFlag(InputFlags::FALSE)->
        setRestartFlags(RestartFlags::WRITE_AND_READ)->
        setDescription("The current history record to be written, on restart from an earlier time allows use to remove invalid future history.");
    }

    virtual ~TimeHistoryOutput() override
    { }

    virtual void SetupDirectoryStructure() override
    {
      // todo: right now if the collector isn't associated with an output operation it will not be able to retrieve a valid
      //       write buffer, and so will be unable to actually collect
      Group * tmp = this->GetGroupByPath(m_collector_path);
      m_collector = Group::group_cast<TimeHistoryCollector*>(tmp);
      GEOSX_ERROR_IF(m_collector == nullptr, "The target of a time history output event must be a collector! " << m_collector_path);
      m_collector->RegisterBufferCall([this]() { return this->m_io->GetBufferHead( ); });
      m_io = new HDFHistIO( m_filename,  m_collector->GetMetadata( ) );
      // its okay for the file to already have this data specificaiton inside it if we are restarting
      //  (write_head > 0) isn't an exact 1-to-1 on that though, as we *could* init the file and file prior to
      //  doing any writes on the file, an IS_RESTART flag of some sort would be nice
      m_io->Init( m_restart );
      m_restart = true;
    }

    /// This method will be called by the event manager if triggered
    virtual void Execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                          real64 const GEOSX_UNUSED_PARAM( dt ),
                          integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                          integer const GEOSX_UNUSED_PARAM( eventCounter ),
                          real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                          dataRepository::Group * GEOSX_UNUSED_PARAM( domain ) ) override
    {
      m_io->Write( );
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

      integer m_restart;

      BufferedHistoryIO * m_io;
      TimeHistoryCollector * m_collector;
  };
}

#endif