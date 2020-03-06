#include "OutputBase.hpp"
#include "fileIO/hdf/HDFFile.hpp"

#include "cxx-utilities/src/Array.hpp" // just for collector

namespace geosx
{
  using namespace traits;

  class TimeHistoryCollector
  {
  public:
    virtual HDFTable GenerateHistorySpec( const string & title, const string & id ) = 0;
    virtual void Collect( real64 const time_n, real64 const dt ) = 0;
    virtual buffer_unit_type * Provide( ) = 0;
  };

  template < typename ARRAY_T, typename ENABLE = void >
  class ArrayTimeHistoryCollector;

  template <typename ARRAY_T >
  class ArrayTimeHistoryCollector< ARRAY_T, typename std::enable_if< is_array< ARRAY_T > >::type > : public TimeHistoryCollector
  {
  public:
    ArrayTimeHistoryCollector ( ARRAY_T const & array ) :
      TimeHistoryCollector(),
      m_arr(array),
      m_offset(0),
      m_collection(array.size() + (ARRAY_T::ndim * sizeof(typename ARRAY_T::index_type)) )
    {}

    virtual HDFTable GenerateHistorySpec( const string & title, const string & id ) override
    {
      HDFTable spec(title,id);
      SpecFromArray(spec,m_arr);
      spec.Finalize();
      return spec;
    }

    virtual void Collect ( real64 const GEOSX_UNUSED_PARAM(time_n), real64 const GEOSX_UNUSED_PARAM(dt) ) override
    {
      buffer_unit_type * buf_head = NULL;
      size_t buffer_size = bufferOps::Pack<false>(buf_head,m_arr.toView());
      m_collection.resize(buffer_size);
      buf_head = &m_collection[0];
      bufferOps::Pack<true>(buf_head,m_arr.toView());
      m_offset = ARRAY_T::ndim * sizeof(typename ARRAY_T::index_type);
    }

    virtual buffer_unit_type * Provide( ) override
    {
      return &m_collection[m_offset];
    }

  private:
    ARRAY_T const & m_arr;
    size_t m_offset;
    std::vector<buffer_unit_type> m_collection;
  };

  class TimeHistory
  {
  public:
    void AddHistory( string const & id, string const & title, TimeHistoryCollector * collector )
    {
      HDFTable spec = collector->GenerateHistorySpec( title, id );
      m_collectors[id] = collector;
      m_time_series.insert(std::make_pair(id,TimeSeries( spec )));
      m_hist_series.insert(std::make_pair(id,HDFTableIO( spec )));
    }
    void UpdateHistories( real64 const time_n, real64 const dt )
    {
      forEachCollector([&](TimeHistoryCollector * collector)
        {
          collector->Collect( time_n, dt );
        }
      );
      forEachTimeSeries([&](HDFTableIO & time_series)
        {
          time_series.BufferRow( reinterpret_cast<buffer_unit_type const *>(&time_n) );
        }
      );
      for( auto kv : m_hist_series )
      {
        kv.second.BufferRow( m_collectors[kv.first]->Provide( ) );
      }

    }
    template < typename LAMBDA >
    inline void forEachHistory( LAMBDA && lambda )
    {
      forEachTimeSeries( lambda );
      forEachHistorySeries( lambda );
    }
    template < typename LAMBDA >
    void forHistory( string const & id, LAMBDA && lambda )
    {
      lambda( m_time_series[id] );
      lambda( m_hist_series[id] );
    }
    template < typename LAMBDA >
    void forEachCollector( LAMBDA && lambda )
    {
      for ( auto kv : m_collectors )
      {
        lambda( kv.second );
      }
    }
    template < typename LAMBDA >
    void forCollector( string const & id, LAMBDA && lambda )
    {
      lambda( m_collectors[id] );
    }
  private:
    template < typename LAMBDA >
    inline void forEachTimeSeries( LAMBDA && lambda )
    {
      for( auto kv : m_time_series )
      {
        lambda( kv.second );
      }
    }
    template < typename LAMBDA >
    inline void forEachHistorySeries( LAMBDA && lambda )
    {
      for( auto kv : m_hist_series )
      {
        lambda( kv.second );
      }
    }

    map<string,TimeHistoryCollector*> m_collectors;
    map<string,HDFTableIO> m_time_series;
    map<string,HDFTableIO> m_hist_series;
  };


// who should own time history? presumably update
  // and output should retrieve and reference it, but there may be
  // multiple time history update events
  // or rather the time history update event may be called many different times for
  // different time histories

  class TimeHistoryUpdate : public OutputBase
  {
  public:
    TimeHistoryUpdate( string const & name, Group * const parent ):
      OutputBase(name,parent),
      m_time_hist()
    {
      // add to data repo
      // m_time_hist
    }

    /// This method will be called by the event manager if triggered
    virtual void Execute( real64 const time_n,
                          real64 const dt,
                          integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                          integer const GEOSX_UNUSED_PARAM( eventCounter ),
                          real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                          dataRepository::Group * GEOSX_UNUSED_PARAM( domain ) ) override
    {
      m_time_hist.UpdateHistories( time_n, dt );
    }

    inline TimeHistory & getTimeHistory( ) { return m_time_hist; }

  private:
    TimeHistory m_time_hist;
  };

  class TimeHistoryOutput : public OutputBase
  {
  public:
    TimeHistoryOutput( string const & hist_filename,
                       string const & name,
                       Group * const parent ):
      OutputBase(name,parent),
      m_thist_filename( hist_filename ),
      m_table_names()
    { }

    virtual ~TimeHistoryOutput() override
    { }

    static string CatalogName() { return "TimeHistoryOutput"; }

    virtual void SetupDirectoryStructure() override
    {

    }

    /// This method will be called by the event manager if triggered
    virtual void Execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                          real64 const GEOSX_UNUSED_PARAM( dt ),
                          integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                          integer const GEOSX_UNUSED_PARAM( eventCounter ),
                          real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                          dataRepository::Group * GEOSX_UNUSED_PARAM( domain ) ) override
    {
      TimeHistory /*&*/ time_hist; // = Wrapper::getRefernce(...)
      HDFFile out(m_thist_filename);
      time_hist.forEachHistory([&out](HDFTableIO & history)
        {
          history.WriteBuffered( out );
        }
      );
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

    void InitHistoryFile()
    {
      TimeHistory /*&*/ time_hist; // = Wrapper::getRefernce(...)
      HDFFile out(m_thist_filename);
      time_hist.forEachHistory( [&](HDFTableIO & history)
        {
          history.CreateInTarget( out );
        }
      );
    }

    struct viewKeysStruct
    {
      static constexpr auto xxxString = "xxxString";
    } timeHistoryOutputViewKeys;

    private:
      string m_thist_filename;
      std::vector<std::string> m_table_names;
  };
}