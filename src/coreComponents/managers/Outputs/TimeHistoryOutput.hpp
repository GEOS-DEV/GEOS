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
      //m_time_history(nullptr),
      m_time_history_filename( ),
      m_time_history_path( )
    {
      // the filename here and the write_head in hdftableio/hdfdataio should be paired in some way
      // as the write head doesn't make sense if the filename changes

      registerWrapper(viewKeysStruct::timeHistoryOutputFilename, &m_time_history_filename, false)->
        setApplyDefaultValue("TimeHistory")->
        setInputFlag(InputFlags::OPTIONAL)->
        setDescription("The filename to which to write time history output.");

      // it would be best to allow multiple targets, but
      // i don't know how required_nonunique works and it isn't used anywhere else
      registerWrapper(viewKeysStruct::timeHistoryOutputTarget, &m_time_history_path, false)->
        setInputFlag(InputFlags::REQUIRED)->
        setDescription("A time history to output to the history file.");
    }

    virtual ~TimeHistoryOutput() override
    { }

    virtual void SetupDirectoryStructure() override
    {
      GetTimeHistoryTarget( );
      m_time_history->Init( m_time_history_filename );
    }

    /// This method will be called by the event manager if triggered
    virtual void Execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                          real64 const GEOSX_UNUSED_PARAM( dt ),
                          integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                          integer const GEOSX_UNUSED_PARAM( eventCounter ),
                          real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                          dataRepository::Group * GEOSX_UNUSED_PARAM( domain ) ) override
    {
      m_time_history->Write( m_time_history_filename );
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

    inline void GetTimeHistoryTarget ( )
    {
      Group * tmp = this->GetGroupByPath(m_time_history_path);
      m_time_history = Group::group_cast<TimeHistory*>(tmp);
      GEOSX_ERROR_IF(m_time_history == nullptr, "The target of a time history output event must be a time history! " << m_time_history_path);
    }

    static string CatalogName() { return "TimeHistoryOutput"; }

    struct viewKeysStruct
    {
      static constexpr auto timeHistoryOutputFilename = "filename";
      static constexpr auto timeHistoryOutputTarget = "target";
    } timeHistoryOutputViewKeys;

    // struct groupKeysStruct
    // {
    //   static constexpr auto timeHistoryOutputHistory
    // } timeHistoryOutputGroups;
    private:
      TimeHistory * m_time_history;
      string m_time_history_filename;
      string m_time_history_path;
  };
}

#endif