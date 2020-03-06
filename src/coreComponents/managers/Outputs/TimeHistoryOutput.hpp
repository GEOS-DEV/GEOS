#include "OutputBase.hpp"

namespace geosx
{

  class TimeHistoryCollector
  {
  public:
    template < typename ARRAY_T >
    typename std::enable_if< is_array< ARRAY_T >, void >::type
    InitCollectionFrom
  };

  class TimeHistory
  {
  public:
    void AddTimeHistory(string const & name, HDFTable table)
    {
      m_tables[name] = HDFTableIO(table);
    }
    HDFTableIO & getTable(string const & name)
    {
      // check the key exists
      return m_tables[name];
    }
  private:
    map<string,HDFTableIO> m_tables;
  };

  class TimeHistoryOutput : public OutputBase
  {
    TimeHistoryOutput( std::string const & name, Group * const parent )
    : m_thist_filename( name )
    , m_tables()
    { }

    virtual ~TimeHistoryOutput() override;

    static string CatalogName() { return "TimeHistoryOutput"; }
    virtual void SetupDirectoryStructure();

    /// This method will be called by the event manager if triggered
    virtual void Execute( real64 const time_n,
                          real64 const dt,
                          integer const cycleNumber,
                          integer const eventCounter,
                          real64 const eventProgress,
                          dataRepository::Group * domain ) override
    {
      // get time history from the catalog
      TimeHistory * t_hist = NULL;
      HDFFile out(m_thist_filename);
      for( auto hist_name : m_table_names )
      {
        HDFTableIO & table = t_hist->getTable( hist_name )
        table.Open( file );
        table.WriteBuffered( );
        table.Close( );
      }
    }

    /// Write one final output as the code exits
    virtual void Cleanup( real64 const time_n,
                          integer const cycleNumber,
                          integer const eventCounter,
                          real64 const eventProgress,
                          dataRepository::Group * domain ) override;

    struct viewKeysStruct
    {
      static constexpr auto xxxString = "xxxString";
    } timeHistoryOutputViewKeys;
    private:
      std::string m_thist_filename;
      vector<std::string> m_table_names;
  };
}