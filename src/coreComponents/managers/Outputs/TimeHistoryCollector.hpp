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

    virtual void AddToSpec( DataSpec & data_spec ) const = 0;
    virtual void Collect( real64 const time_n, real64 const dt, buffer_unit_type *& buffer ) = 0;

    void InitSpec( DataSpec & data_spec ) const
    {
      data_spec.Append<real64>(1,1,"Time");
      AddToSpec(data_spec);
      data_spec.Finalize( );
    }

    virtual void Execute( real64 const time_n,
                          real64 const dt,
                          integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                          integer const GEOSX_UNUSED_PARAM( eventCounter ),
                          real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                          dataRepository::Group * GEOSX_UNUSED_PARAM( domain ) ) override
    {
      buffer_unit_type * buffer = m_buffer_call();
      collect_with_time( time_n, dt, buffer );
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
      // it would be best to allow multiple targets, but
      // i don't know how required_nonunique works and it isn't used anywhere else
      registerWrapper(ScalarCollector<VALUE_T>::viewKeysStruct::scalarPath, &m_scalar_path, false)->
        setInputFlag(InputFlags::OPTIONAL)->
        setDescription("A scalar value to collect for time history output.");
    }

    virtual void AddToSpec( DataSpec & data_spec ) const override
    {
      data_spec.Append(1,1,sizeof(VALUE_T),typeid(VALUE_T),m_data_title);
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

  template < typename ARRAY_T, typename ENABLE = void >
  class ArrayTimeHistoryCollector;

  template <typename ARRAY_T >
  class ArrayTimeHistoryCollector< ARRAY_T, typename std::enable_if< is_array< ARRAY_T > && can_history_io< ARRAY_T::value_type> >::type > : public TimeHistoryCollector
  {
  public:
    ArrayTimeHistoryCollector ( string const & name, Group * parent )
    : TimeHistoryCollector( name, parent )
    { }

    virtual void AddToSpec ( DataSpec & data_spec ) const override
    {
      typename ARRAY_T::view_type & av = Group::getReference(m_array_path);
      AppendArraySpec(data_spec,av);
    }

    virtual void Collect ( real64 const GEOSX_UNUSED_PARAM(time_n), real64 const GEOSX_UNUSED_PARAM(dt), buffer_unit_type *& buffer ) override
    {
      typename ARRAY_T::view_type & av = Group::getReference(m_array_path);
      bufferOps::PackDevice<true>(buffer,av);
    }
  private:
    string m_array_path;
  };

  template < typename ARRAY_T, typename INDEX_ARRAY_T, typename ENABLE = void >
  class ArrayIndexedTimeHistoryCollector;

  template <typename ARRAY_T, typename INDEX_ARRAY_T >
  class ArrayIndexedTimeHistoryCollector< ARRAY_T, INDEX_ARRAY_T, typename std::enable_if< is_array< ARRAY_T > &&
                                                                                           can_history_io< ARRAY_T::value_type> &&
                                                                                           is_array< INDEX_ARRAY_T > &&
                                                                                           LvArray::is_integer< typename INDEX_ARRAY_T::value_type >::value >::type > : public TimeHistoryCollector
  {
  public:
    ArrayIndexedTimeHistoryCollector ( string const & name, Group * parent ) :
      TimeHistoryCollector( name, parent )
    { }

    virtual void AddToSpec ( DataSpec & data_spec ) const override
    {
      // Probe<INDEX_ARRAY_T> & probe = Group::getReference(m_probe_path);
      typename INDEX_ARRAY_T::view_type /*&*/ idxv;// = Group::GetGroupByPath(m_probe_path).getReference(Probe<INDEX_ARRAY_T>::viewKeysStruct::indicesKey);
      typename ARRAY_T::view_type & av = Group::getReference(m_array_path);
      AppendArrayIndicesSpec(data_spec,av,idxv.size( ));
    }

    virtual void Collect ( real64 const GEOSX_UNUSED_PARAM(time_n), real64 const GEOSX_UNUSED_PARAM(dt), buffer_unit_type *& buffer ) override
    {
      // Probe<INDEX_ARRAY_T> & probe = Group::getReference(m_probe_path);
      typename INDEX_ARRAY_T::view_type /*&*/ idxv;// = Group::GetGroupByPath(m_probe_path).getReference(Probe<INDEX_ARRAY_T>::viewKeysStruct::indicesKey);
      typename ARRAY_T::view_type & av = Group::getReference(m_array_path);
      bufferOps::PackByIndexDevice<true>(buffer,av,idxv.size( ));
    }
  private:
    string m_probe_path;
    string m_array_path;
  };
}

#endif