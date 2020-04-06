#ifndef GEOSX_TIMEHISTORYCOLLECTOR_HPP_
#define GEOSX_TIMEHISTORYCOLLECTOR_HPP_

#include "dataRepository/Group.hpp"
#include "managers/TimeHistory/HistoryDataSpec.hpp"
#include "dataRepository/BufferOpsDevice.hpp"

#include <functional>

namespace geosx
{
  using namespace dataRepository;

  // todo make this a group, required
  class TimeHistoryCollector : public ExecutableGroup
  {
  public:
    TimeHistoryCollector( string const & name, Group * parent );

    virtual HistoryMetadata GetMetadata( ) const = 0;
    virtual void Collect( real64 const time_n, real64 const dt, buffer_unit_type *& buffer ) = 0;

    virtual void Execute( real64 const time_n,
                          real64 const dt,
                          integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                          integer const GEOSX_UNUSED_PARAM( eventCounter ),
                          real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                          dataRepository::Group * GEOSX_UNUSED_PARAM( domain ) ) override
    {
      buffer_unit_type * buffer = m_buffer_call();
      Collect( time_n, dt, buffer );
    }

    void RegisterBufferCall( std::function<buffer_unit_type*()> buffer_call )
    {
      m_buffer_call = buffer_call;
    }

    static string CatalogName() { return "TimeHistoryCollector"; }
    struct viewKeysStruct
    {
      static constexpr auto dataTitleKey = "title";
    } timeHistoryCollectorKeys;
    struct groupKeysStruct
    {
      // static constexpr auto timeCollectorKey = "TimeCollector";
    } timeHistoryGroupKeys;
  protected:
    string m_data_title;
  private:
    void collect_with_time( real64 time_n, real64 dt, buffer_unit_type *& buffer )
    {
      bufferOps::PackPointerDevice< true >( buffer, &time_n, 1 );
      Collect(time_n, dt, buffer);
    }
    std::function<buffer_unit_type*()> m_buffer_call;
  };

  template < typename VALUE_TYPE, typename ENABLE = void >
  class ScalarCollector;

  template < typename VALUE_T >
  class ScalarCollector< VALUE_T, typename std::enable_if< can_history_io< VALUE_T > >::type > : public TimeHistoryCollector
  {
  public:
    ScalarCollector( string const & name, Group * parent )
      : TimeHistoryCollector(name,parent)
      , m_scalar_path("")
    {
      registerWrapper(ScalarCollector<VALUE_T>::viewKeysStruct::scalarPath, &m_scalar_path, false)->
        setInputFlag(InputFlags::REQUIRED)->
        setDescription("A scalar value to collect for time history output.");
    }

    virtual HistoryMetadata GetMetadata( ) const override
    {
      return HistoryMetadata(m_data_title,1,std::type_index(typeid(VALUE_T)));
    }
    virtual void Collect( real64 const GEOSX_UNUSED_PARAM(time_n), real64 const GEOSX_UNUSED_PARAM(dt), buffer_unit_type *& buffer ) override
    {
      VALUE_T & value = Group::getReference( m_scalar_path );
      bufferOps::PackPointerDevice< true >( buffer, &value, 1 );
    }
    static string CatalogName() { return "ScalarCollector"; }

    struct viewKeysStruct
    {
      static constexpr auto scalarPath = "target";
    } scalarCollectorKeys;

  private:
    string m_scalar_path;
  };


  class ArrayCollector : public TimeHistoryCollector
  {
    public:
    ArrayCollector ( string const & name, Group * parent )
    : TimeHistoryCollector( name, parent )
    {
      registerWrapper(ArrayCollector<ARRAY_T>::viewKeysStruct::m_array_path, &m_array_path, false)->
       setInputFlag(InputFlags::REQUIRED)->
       setDescription("The path of an array to collect for time history output.");
    }

    static string CatalogName() { return "ArrayCollector"; }

    struct viewKeysStruct
    {
      static constexpr auto arrayPath = "target";
    } arrayCollectorKeys;

  protected:
    string m_array_path;
  };

  template < typename ARRAY_T, typename ENABLE = void >
  class ArrayCollectorT;

  template <typename ARRAY_T >
  class ArrayCollectorT< ARRAY_T, typename std::enable_if< is_array< ARRAY_T > && can_history_io< ARRAY_T::value_type> >::type > : public ArrayCollector
  {
  public:
    ArrayCollectorT ( string const & name, Group * parent )
    : ArrayCollector( name, parent )
    { }

    virtual HistoryMetadata GetMetadata( ) const override
    {
      typename ARRAY_T::view_type & av = Group::getReference(m_array_path);
      return ArrayMetadata(m_data_title,av);
    }

    virtual void Collect ( real64 const GEOSX_UNUSED_PARAM(time_n), real64 const GEOSX_UNUSED_PARAM(dt), buffer_unit_type *& buffer ) override
    {
      typename ARRAY_T::view_type & av = Group::getReference(m_array_path);
      bufferOps::PackDevice<true>(buffer,av);
    }
  };



  class IndexedArrayCollector : public TimeHistoryCollector
  {
  public:
    IndexedArrayCollector ( string const & name, Group * parent ) :
      TimeHistoryCollector( name, parent )
    {
      registerWrapper(ArrayCollector<ARRAY_T>::viewKeysStruct::m_array_path, &m_array_path, false)->
       setInputFlag(InputFlags::REQUIRED)->
       setDescription("The path of an array to collect for time history output.");
      registerWrapper(ArrayCollector<ARRAY_T>::viewKeysStruct::m_probe_path, &m_probe_path, false)->
       setInputFlag(InputFlags::REQUIRED)->
       setDescription("The path of a probe to supply indices for time history collection.");
    }

    static string CatalogName() { return "ArrayIndexedCollector"; }

    struct viewKeysStruct
    {
      static constexpr auto arrayPath = "target";
    } arrayCollectorKeys;

  protected:
    string m_probe_path;
    string m_array_path;
  };


  template < typename ARRAY_T, typename INDEX_ARRAY_T, typename ENABLE = void >
  class IndexedArrayCollectorT;

  template <typename ARRAY_T, typename INDEX_ARRAY_T >
  class IndexedArrayCollectorT< ARRAY_T, INDEX_ARRAY_T, typename std::enable_if< is_array< ARRAY_T > &&
                                                                                           can_history_io< ARRAY_T::value_type> &&
                                                                                           is_array< INDEX_ARRAY_T > &&
                                                                                           LvArray::is_integer< typename INDEX_ARRAY_T::value_type >::value >::type > : public IndexedArrayCollector
  {
  public:
    IndexedArrayCollectorT ( string const & name, Group * parent )
      : IndexedArrayCollector( name, parent )
    { }

    virtual HistoryMetadata GetMetadata( ) const override
    {
      // Probe<INDEX_ARRAY_T> & probe = Group::getReference(m_probe_path);
      typename INDEX_ARRAY_T::view_type /*&*/ idxv;// = probe.getReference(Probe<INDEX_ARRAY_T>::viewKeysStruct::indicesKey);
      typename ARRAY_T::view_type & av = Group::getReference(m_array_path);
      return ArrayIndicesMetadata(m_data_title,av,idxv.size( ));
    }

    virtual void Collect ( real64 const GEOSX_UNUSED_PARAM(time_n), real64 const GEOSX_UNUSED_PARAM(dt), buffer_unit_type *& buffer ) override
    {
      // Probe<INDEX_ARRAY_T> & probe = Group::getReference(m_probe_path);
      typename INDEX_ARRAY_T::view_type /*&*/ idxv;// = probe.getReference(Probe<INDEX_ARRAY_T>::viewKeysStruct::indicesKey);
      typename ARRAY_T::view_type & av = Group::getReference(m_array_path);
      bufferOps::PackByIndexDevice<true>(buffer,av,idxv.size( ));
    }

  };
}

#endif