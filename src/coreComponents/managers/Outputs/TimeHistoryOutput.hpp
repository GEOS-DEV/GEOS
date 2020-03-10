#include "OutputBase.hpp"
#include "fileIO/hdf/HDFFile.hpp"

#include "cxx-utilities/src/Array.hpp" // just for collector

namespace geosx
{
  using namespace traits;

  class TimeHistoryCollector
  {
  public:
    virtual void AddToSpec( HDFTable table ) const = 0;
    virtual void Collect( real64 const time_n, real64 const dt ) = 0;
    virtual buffer_unit_type * Provide( ) = 0;
    virtual size_t size() const = 0;
  };

  class CollectorSet : public TimeHistoryCollector
  {
  public:
    void AddCollector( TimeHistoryCollector * to_add )
    {
      m_collectors.insert(to_add);
    }

    virtual void AddToSpec( HDFTable spec ) const override
    {
      forEachCollector([&spec](const TimeHistoryCollector * coll)
        {
          coll->AddToSpec(spec);
        });
    }

    virtual void Collect( real64 const time_n, real64 const dt ) override
    {
      size_t total_size = 0;
      forEachCollector([&time_n,&dt,&total_size](TimeHistoryCollector * coll)
        {
          coll->Collect(time_n,dt);
          total_size += coll->size( );
        });
      m_buffer.resize(total_size);
      size_t offset = 0;
      buffer_unit_type * buffer_head = &m_buffer[0];
      forEachCollector([&offset,&buffer_head](TimeHistoryCollector * coll)
        {
          size_t size = coll->size( );
          memcpy(buffer_head + offset,coll->Provide( ),size);
          offset += size;
        }) ;
    }

    virtual buffer_unit_type * Provide( ) override
    {
      return &m_buffer[0];
    }

    virtual size_t size( ) const override
    {
      size_t total_size = 0;
      forEachCollector([&total_size](const TimeHistoryCollector * coll)
        {
          total_size += coll->size();
        });
      return total_size;
    }

  private:

    template < typename LAMBDA >
    void forEachCollector( LAMBDA && lambda )
    {
      for( TimeHistoryCollector * coll : m_collectors )
      {
        lambda(coll);
      }
    }

    template < typename LAMBDA >
    void forEachCollector( LAMBDA && lambda ) const
    {
      for( const TimeHistoryCollector * coll : m_collectors )
      {
        lambda(coll);
      }
    }


    set<TimeHistoryCollector*> m_collectors;
    std::vector<buffer_unit_type> m_buffer;
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
    {
      // this points past the packing metadata
      // todo: move this calc into bufferOps somewhere
      m_offset = ARRAY_T::ndim * sizeof(typename ARRAY_T::index_type);
    }

    virtual void AddToSpec( HDFTable spec ) const override
    {
      spec.AddArrayCol(m_arr);
    }

    virtual void Collect ( real64 const GEOSX_UNUSED_PARAM(time_n), real64 const GEOSX_UNUSED_PARAM(dt) ) override
    {
      buffer_unit_type * buf_head = NULL;
      size_t buffer_size = bufferOps::PackDevice<false>(buf_head,m_arr.toView());
      m_collection.resize(buffer_size);
      buf_head = &m_collection[0];
      bufferOps::PackDevice<true>(buf_head,m_arr.toView());
    }

    virtual buffer_unit_type * Provide( ) override
    {
      return &m_collection[m_offset];
    }

    virtual size_t size( ) const override
    {
      return m_collection.size() - m_offset;
    }

  private:
    ARRAY_T const & m_arr;
    size_t m_offset;
    std::vector<buffer_unit_type> m_collection;
  };

  template < typename ARRAY_T, typename INDEX_ARRAY_T, typename ENABLE = void >
  class ArrayIndexedTimeHistoryCollector;

  template <typename ARRAY_T, typename INDEX_ARRAY_T >
  class ArrayIndexedTimeHistoryCollector< ARRAY_T, INDEX_ARRAY_T, typename std::enable_if< is_array< ARRAY_T > &&
                                                                                           is_array< INDEX_ARRAY_T > &&
                                                                                           LvArray::is_integer< typename INDEX_ARRAY_T::value_type >::value >::type > : public TimeHistoryCollector
  {
  public:
    ArrayIndexedTimeHistoryCollector ( ARRAY_T const & array, INDEX_ARRAY_T const & idx_arr ) :
      TimeHistoryCollector(),
      m_arr(array),
      m_idxs(idx_arr),
      m_offset(0),
      m_collection(array.size() + (ARRAY_T::ndim * sizeof(typename ARRAY_T::index_type)) )
    {}

    virtual void AddToSpec ( HDFTable spec ) const override
    {
      spec.AddArrayIndicesCol(m_arr,m_idxs.size( ));
    }

    virtual void Collect ( real64 const GEOSX_UNUSED_PARAM(time_n), real64 const GEOSX_UNUSED_PARAM(dt) ) override
    {
      buffer_unit_type * buf_head = NULL;
      size_t buffer_size = bufferOps::PackByIndexDevice<false>(buf_head,m_arr.toView(),m_idxs);
      m_collection.resize(buffer_size);
      buf_head = &m_collection[0];
      bufferOps::PackByIndexDevice<true>(buf_head,m_arr.toView(),m_idxs);
      // this points past the packing metadata
      // todo: move this calc into bufferOps somewhere
      m_offset = ARRAY_T::ndim * sizeof(typename ARRAY_T::index_type);
    }

    virtual buffer_unit_type * Provide ( ) override
    {
      return &m_collection[m_offset];
    }

    virtual size_t size( ) const override
    {
      return m_collection.size() - m_offset;
    }

  private:
    ARRAY_T const & m_arr;
    INDEX_ARRAY_T const & m_idxs;
    size_t m_offset;
    std::vector<buffer_unit_type> m_collection;
  };

  // move into hdffile once the collectors are in a seperate header
  HDFTableIO InitTimeHistoryIO( string const & name, string const & id, TimeHistoryCollector * coll )
  {
    HDFTable spec = InitHistoryTable( name,id );
    coll->AddToSpec( spec );
    spec.Finalize();
    return HDFTableIO( spec );
  }

  class TimeHistory
  {
  public:
    TimeHistory( string const & name, string const & id, TimeHistoryCollector * coll )
      : m_collector(coll)
      , m_hist_io(InitTimeHistoryIO(name,id,coll))
    { }
    void InitFile( string const & filename )
    {
      HDFFile out_file( filename );
      m_hist_io.CreateInTarget( out_file );
    }
    void Update( real64 const time_n, real64 const dt )
    {
      m_collector->Collect( time_n, dt );
      m_hist_io.BufferRow( m_collector->Provide() );
    }
    void WriteToFile( string const & filename )
    {
      HDFFile out_file( filename );
      m_hist_io.WriteBuffered( out_file );
    }
  private:
    TimeHistoryCollector * m_collector;
    HDFTableIO m_hist_io;
  };


// who should own time history? presumably update
  // and output should retrieve and reference it, but there may be
  // multiple time history update events
  // or rather the time history update event may be called many different times for
  // different time histories

  class TimeHistoryUpdate : public OutputBase
  {
  public:
    TimeHistoryUpdate( string const & target, string const & name, Group * const parent ):
      OutputBase(name,parent),
      m_target(target)
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
      TimeHistory /*&*/ * target = nullptr;
      target->Update( time_n, dt );
    }

    //inline TimeHistory & getTimeHistoryTarget( ) { return m_time_hist; }

  private:
    string m_target;
  };

  class TimeHistoryOutput : public OutputBase
  {
  public:
    TimeHistoryOutput( string const & hist_filename,
                       string const & target,
                       string const & name,
                       Group * const parent ):
      OutputBase(name,parent),
      m_thist_filename( hist_filename ),
      m_target(target)
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
      TimeHistory /*&*/ * time_hist = nullptr; // = Wrapper::getReference(m_target...t...)
      time_hist->WriteToFile( m_thist_filename );
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
      TimeHistory /*&*/ * time_hist = nullptr; // = Wrapper::getRefernce(...)
      time_hist->InitFile( m_thist_filename );
    }

    struct viewKeysStruct
    {
      static constexpr auto xxxString = "xxxString";
    } timeHistoryOutputViewKeys;

    private:
      string m_thist_filename;
      string m_target;
  };
}