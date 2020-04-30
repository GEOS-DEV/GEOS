#ifndef GEOSX_HistoryCollection_HPP_
#define GEOSX_HistoryCollection_HPP_

#include "managers/Tasks/TaskBase.hpp"
#include "managers/TimeHistory/HistoryDataSpec.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/ProblemManager.hpp"
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

    virtual HistoryMetadata GetMetadata( Group * domain ) const = 0;
    virtual void Collect( Group * domain, real64 const time_n, real64 const dt, buffer_unit_type *& buffer ) = 0;

    virtual void Execute( real64 const time_n,
                          real64 const dt,
                          integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                          integer const GEOSX_UNUSED_PARAM( eventCounter ),
                          real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                          Group * domain ) override
    {
      // std::function defines the == and =! comparable against nullptr_t to check the
      //  function pointer is actually assigned (an error would be thrown on the call attempt even so)
      GEOSX_ERROR_IF( m_buffer_call == nullptr,
                      "History collection buffer retrieval function is unassigned, did you declare a related TimeHistoryOutput event?" );
      // using GEOSX_ERROR_IF_EQ causes type issues since the values are used in iostreams
      buffer_unit_type * buffer = m_buffer_call();
      Collect( domain, time_n, dt, buffer );

      int rank = MpiWrapper::Comm_rank();
      if ( rank == 0 )
      {
        buffer_unit_type * time_buffer = m_time_buffer_call();
        memcpy( time_buffer, &time_n, sizeof(decltype(time_n)) );
      }
    }

    void RegisterBufferCall( std::function<buffer_unit_type*()> buffer_call )
    {
      m_buffer_call = buffer_call;
    }

    HistoryMetadata GetTimeMetadata( ) const
    {
      return HistoryMetadata("Time",1,std::type_index(typeid(real64)));
    }

    void RegisterTimeBufferCall( std::function<buffer_unit_type*()> time_buffer_call )
    {
      m_time_buffer_call = time_buffer_call;
    }

  private:
    std::function<buffer_unit_type*()> m_time_buffer_call;
    std::function<buffer_unit_type*()> m_buffer_call;
  };

  class PackCollection : public HistoryCollection
  {
    public:
    PackCollection ( string const & name, Group * parent )
    : HistoryCollection( name, parent )
    {
      registerWrapper(PackCollection::viewKeysStruct::objectPath, &m_object_path, false)->
       setInputFlag(InputFlags::REQUIRED)->
       setDescription("The name of the object from which to retrieve field values.");

      registerWrapper(PackCollection::viewKeysStruct::fieldName, &m_field_name, false)->
        setInputFlag(InputFlags::REQUIRED)->
        setDescription("The name of the (packable) field associated with the specified object to retrieve data from");

      registerWrapper(PackCollection::viewKeysStruct::setNames, &m_set_names, false)->
        setInputFlag(InputFlags::OPTIONAL)->
        setDescription("The set(s) for which to retrieve data.");
    }

    virtual HistoryMetadata GetMetadata( Group * problem_group ) const override
    {
      ProblemManager const * pm = dynamicCast< ProblemManager const * >( problem_group );
      DomainPartition const & domain = *pm->getDomainPartition( );
      MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

      Group const * target_object = meshLevel.GetGroupByPath( m_object_path );
      WrapperBase const * m_target = target_object->getWrapperBase( m_field_name );
      GEOSX_ERROR_IF( ! m_target->isPackable(), "The object targeted for collection must be packable!" );
      return m_target->getBufferedIOMetadata(); // probably a better name for this function
    }

    virtual void Collect( Group * domain_group,
                          real64 const GEOSX_UNUSED_PARAM(time_n),
                          real64 const GEOSX_UNUSED_PARAM(dt),
                          buffer_unit_type*& buffer ) override
    {
      DomainPartition const & domain = dynamicCast< DomainPartition const &>( *domain_group );
      MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

      // todo: assert(m_object_path == "nodeManager" || "edgeManager" || "FaceManager" || "elementSubRegions" ...
      Group const * target_object = meshLevel.GetGroupByPath( m_object_path );
      if ( m_set_names.size( ) > 0 )
      {
        Group const * set_group = target_object->GetGroup( ObjectManagerBase::groupKeyStruct::setsString );
        for( auto & set_name : m_set_names )
        {
          WrapperBase const * set_target = set_group->getWrapperBase( set_name );
          if( set_target != nullptr )
          {
            set_target->Pack( buffer, false, true );
          }
        }
      }
      else
      {
        WrapperBase const * m_target = target_object->getWrapperBase( m_field_name );
        m_target->Pack( buffer, false, true );
      }
    }

    static string CatalogName() { return "PackCollection"; }

    struct viewKeysStruct
    {
      static constexpr auto objectPath = "objectPath";
      static constexpr auto fieldName = "fieldName";
      static constexpr auto setNames = "setNames";
    } keys;

  protected:
    string m_object_path;
    string m_field_name;
    string_array m_set_names;
  };

}

#endif
