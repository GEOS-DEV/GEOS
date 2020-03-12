#ifndef GEOSX_TIME_HISTORY_HPP
#define GEOSX_TIME_HISTORY_HPP

#include "TimeHistoryCollector.hpp"
#include "TimeHistoryOutput.hpp"

namespace geosx
{
  class TimeHistory : public ExecutableGroup
  {
  public:
    TimeHistory( string const & name, Group * parent )
      : ExecutableGroup( name, parent )
      , m_title()
      , m_id()
      , m_data_spec(new HDFTable)
      , m_collector(nullptr)
      , m_io(nullptr)
    {
      registerWrapper(viewKeys::timeHistoryFilename, &m_title, false)->
        setApplyDefaultValue("TimeHistory")->
        setInputFlag(InputFlags::REQUIRED)->
        setDescription("The title of the data set containing the time history.");

      // it would be best to allow multiple targets, but
      // i don't know how required_nonunique works and it isn't used anywhere else
      registerWrapper(viewKeys::timeHistoryTarget, &m_tid, false)->
       setApplyDefaultValue("time_hist")
        setInputFlag(InputFlags::REQUIRED)->
        setDescription("A short unique identifier for the data related to the time history.");

      // RegisterGroup<TimeHistoryCollector>(keys::collectorKey);
      // RegisterGroup<BufferedHistoryIO>(keys::ioKey)
      registerWrapper(groupKeys::dataSpec, m_data_spec);
    }

    virtual void Init( string const & target_name )
    {
      GetGroups( );
      m_data_spec->SetTitleID( m_title, m_id );
      m_collector->AddToSpec( *m_data_spec );
      m_data_spec->Finalize();
      m_io->Init( target_name, *m_data_spec );
    }
    virtual void Update( real64 const time_n, real64 const dt )
    {
      m_collector->Collect( time_n, dt, m_io->GetRowHead( DataSpec const * spec ) );
    }
    virtual void Write( string const & target_name )
    {
      m_io.Write( target_name, *m_data_spec );
    }

    inline void GetGroups()
    {
      Group * tmp = Group::GetGroup(keys::collectorKey);
      m_collector = Group::group_cast<TimeHistoryCollector*>(tmp);
      // GEOSX_ERROR_IF(m_time_history == nullptr, "The target of a time history output event must be a time history! " << m_time_history_path);
      tmp = Group::GetGroup(keys::IO);
      m_io = Group::group_cast<BufferedHistoryIO*>(tmp);
      // GEOSX_ERROR_IF(m_time_history == nullptr, "The target of a time history output event must be a time history! " << m_time_history_path);
    }

    inline DataSpec const * GetDataSpec()
    {
      return m_data_spec;
    }

    static string CatalogName() { return "TimeHistory"; }

    struct viewKeys
    {
      static constexpr auto title = "title";
      static constexpr auto id = "id";
    } timeHistoryViewKeys;

    struct groupKeys
    {
      static constexpr auto dataSpec = "DataSpec";
      static constexpr auto collectorKey = "Collector";
      static constexpr auto ioKey = "IO";
    } timeHistoryGroupKeys;
  protected:
    string m_title;
    string m_id;
    DataSpec * m_data_spec;
    TimeHistoryCollector * m_collector;
    BufferedHistoryIO * m_io;
  };
}

#endif