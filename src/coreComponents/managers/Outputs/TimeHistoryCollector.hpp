#ifndef GEOSX_HistoryCollection_HPP_
#define GEOSX_HistoryCollection_HPP_

#include "managers/Tasks/TaskBase.hpp"
#include "managers/TimeHistory/HistoryDataSpec.hpp"
#include "dataRepository/BufferOpsDevice.hpp"

#include <functional>

namespace geosx
{
  using namespace dataRepository;

  class HistoryCollection : public TaskBase
  {
  public:
    HistoryCollection( string const & name, Group * parent ) :
      TaskBase( name, parent ),
      m_buffer_call()
    {   }

    virtual HistoryMetadata GetMetadata( ) const = 0;
    virtual void Collect( real64 const time_n, real64 const dt, buffer_unit_type *& buffer ) = 0;

    virtual void Execute( real64 const time_n,
                          real64 const dt,
                          integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                          integer const GEOSX_UNUSED_PARAM( eventCounter ),
                          real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                          dataRepository::Group * GEOSX_UNUSED_PARAM( domain ) ) override
    {
      //assert that m_buffer_call is associated correctly
      buffer_unit_type * buffer = m_buffer_call();
      Collect( time_n, dt, buffer );
    }

    void RegisterBufferCall( std::function<buffer_unit_type*()> buffer_call )
    {
      m_buffer_call = buffer_call;
    }

  private:
    std::function<buffer_unit_type*()> m_buffer_call;
  };

  class TimeCollection : public HistoryCollection
  {
  public:
    TimeCollection( string const & name, Group * parent ) :
      HistoryCollection( name, parent )
    {  }
    virtual HistoryMetadata GetMetadata( ) const
    {
      return HistoryMetadata("Time",1,std::type_index(typeid(real64)));
    }
    virtual void Collect( real64 const time_n, real64 const GEOSX_UNUSED_PARAM( dt ), buffer_unit_type *& buffer )
    {
      bufferOps::Pack<true>( buffer, time_n );
    }
  };

  class PackCollection : public HistoryCollection
  {
    public:
    PackCollection ( string const & name, Group * parent )
    : HistoryCollection( name, parent )
    {
      registerWrapper(PackCollection::viewKeysStruct::arrayPath, &m_target_path, false)->
       setInputFlag(InputFlags::REQUIRED)->
       setDescription("The path of a packable data structure to collect for time history output.");
    }

    virtual HistoryMetadata GetMetadata( ) const override
    {
      WrapperBase const * target = Group::getWrapperBase( m_target_path );
      return target->getBufferedIOMetadata();
    }

    virtual void Collect( real64 const GEOSX_UNUSED_PARAM(time_n), real64 const GEOSX_UNUSED_PARAM(dt), buffer_unit_type*& buffer ) override
    {
      WrapperBase const * target = Group::getWrapperBase( m_target_path );
      // assert(target->isPackable());
      target->Pack( buffer );
    }

    static string CatalogName() { return "PackCollection"; }

    struct viewKeysStruct
    {
      static constexpr auto arrayPath = "target";
    } keys;

  protected:
    string m_target_path;
  };

  class PackByIndexCollection : public PackCollection
  {
  public:
    PackByIndexCollection ( string const & name, Group * parent ) :
      PackCollection( name, parent )
    {
      registerWrapper(PackByIndexCollection::viewKeysStruct::probePath, &m_indexer_path, false)->
       setInputFlag(InputFlags::REQUIRED)->
       setDescription("The path of a probe to supply indices for time history collection.");
    }

    virtual HistoryMetadata GetMetadata( ) const override
    {
      WrapperBase const * target = Group::getWrapperBase( m_target_path );
      return target->getBufferedIOMetadata();
      // todo : reduce first dimension extent of metadata object based on number of indices
      // typename INDICES_T::view_type & iv = Group::getReference(m_indexer_path);
      // return ArrayByIndexMetadata(m_data_title,av,num_indices);
    }

    virtual void Collect( real64 const GEOSX_UNUSED_PARAM(time_n), real64 const GEOSX_UNUSED_PARAM(dt), buffer_unit_type*& buffer ) override
    {
      WrapperBase const * target = Group::getWrapperBase( m_target_path );
      // assert(target->isPackable());
      //ArrayIndexer const & probe = Group::getReference<Probe>
      //arrayView1d< localIndex const > const & indices = probe.getIndices( );
      //target->PackByIndex( buffer, indices );
      target->Pack( buffer );
    }

    static string CatalogName() { return "PackByIndexCollection"; }

    struct viewKeysStruct
    {
      static constexpr auto probePath = "probe";
    } keys;

  protected:
    string m_indexer_path;
  };
}

#endif
