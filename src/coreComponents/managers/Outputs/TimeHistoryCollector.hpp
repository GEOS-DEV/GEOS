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

  // todo : there is some refactoring that can be done to these classes to simplify the usage
  //        re: writing/collecting metadata once and repeated history collection
  class HistoryCollection : public TaskBase
  {
  public:
    HistoryCollection( string const & name, Group * parent ) :
      TaskBase( name, parent ),
      m_buffer_call()
    {   }

    virtual HistoryMetadata GetMetadata( Group * domain ) const
    {
      GEOSX_UNUSED_VAR( domain );
      return HistoryMetadata( "null", 0, std::type_index(typeid(nullptr)) );
    }

    virtual localIndex GetNumMetaCollectors( ) const { return 0; }
    virtual std::unique_ptr<HistoryCollection> GetMetaCollector( Group * problem_group, localIndex meta_idx ) const
    {
      GEOSX_UNUSED_VAR( problem_group );
      GEOSX_UNUSED_VAR( meta_idx );
      return std::unique_ptr<HistoryCollection>(nullptr);
    }

    virtual void Collect( Group * domain, real64 const time_n, real64 const dt, buffer_unit_type *& buffer )
    {
      GEOSX_UNUSED_VAR( domain );
      GEOSX_UNUSED_VAR( time_n );
      GEOSX_UNUSED_VAR( dt );
      GEOSX_UNUSED_VAR( buffer );
    }

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

  // todo : this makes some small assumptions about the output file indices (dense, contiguous, ordered by rank, etc)
  class SetMetaCollection : public HistoryCollection
  {
  public:
    SetMetaCollection( string const & object_path, string const & set_name, localIndex set_index_offset ) :
      HistoryCollection( "SetMetaCollection", nullptr ),
      m_object_path( object_path ),
      m_set_name( set_name ),
      m_local_index_offset( set_index_offset ),
      m_global_index_offset( 0 )
    { }

    void SetGlobalOffset( localIndex offset ) { m_global_index_offset = offset; }

    virtual HistoryMetadata GetMetadata( Group * problem_group ) const override
    {
      ProblemManager const * pm = dynamicCast< ProblemManager const * >( problem_group );
      DomainPartition const & domain = *pm->getDomainPartition( );
      MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
      Group const * target_object = meshLevel.GetGroupByPath( m_object_path );
      Group const * set_group = target_object->GetGroup( ObjectManagerBase::groupKeyStruct::setsString );
      dataRepository::Wrapper< SortedArray< localIndex > > const * const set_wrapper = set_group->getWrapper< SortedArray< localIndex > >( m_set_name );
      return set_wrapper->getBufferedIOMetadata( );
    }

    virtual void Collect( Group * domain_group,
                          real64 const GEOSX_UNUSED_PARAM(time_n),
                          real64 const GEOSX_UNUSED_PARAM(dt),
                          buffer_unit_type*& buffer ) override
    {
      DomainPartition const & domain = dynamicCast< DomainPartition const &>( *domain_group );
      MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
      Group const * target_object = meshLevel.GetGroupByPath( m_object_path );
      Group const * set_group = target_object->GetGroup( ObjectManagerBase::groupKeyStruct::setsString );
      dataRepository::Wrapper< SortedArray< localIndex > > const * const set_wrapper = set_group->getWrapper< SortedArray< localIndex > >( m_set_name );
      if( set_wrapper != nullptr )
      {
        SortedArrayView< localIndex const > const & set = set_wrapper->reference();
        localIndex sz = set.size( );
        array1d< localIndex > meta_idx( sz );
        for( localIndex ii = 0; ii < sz; ++ii )
        {
          meta_idx[ii] = m_global_index_offset + m_local_index_offset + ii;
        }
        bufferOps::PackDevice< true >( buffer, meta_idx.toView() );
      }
    }
  protected:
    string m_object_path;
    string m_set_name;

    localIndex m_local_index_offset;
    localIndex m_global_index_offset;
  };

  class PackCollection : public HistoryCollection
  {
    public:
    PackCollection ( string const & name, Group * parent )
    : HistoryCollection( name, parent )
    {
      registerWrapper(PackCollection::viewKeysStruct::objectPath, &m_object_path)->
       setInputFlag(InputFlags::REQUIRED)->
       setDescription("The name of the object from which to retrieve field values.");

      registerWrapper(PackCollection::viewKeysStruct::fieldName, &m_field_name)->
        setInputFlag(InputFlags::REQUIRED)->
        setDescription("The name of the (packable) field associated with the specified object to retrieve data from");

      registerWrapper(PackCollection::viewKeysStruct::setNames, &m_set_names)->
        setInputFlag(InputFlags::OPTIONAL)->
        setDescription("The set(s) for which to retrieve data.");
    }

    virtual HistoryMetadata GetMetadata( Group * problem_group ) const override
    {
      localIndex num_indices = CountLocalSetIndices( problem_group );
      ProblemManager const * pm = dynamicCast< ProblemManager const * >( problem_group );
      DomainPartition const & domain = *pm->getDomainPartition( );
      MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
      Group const * target_object = meshLevel.GetGroupByPath( m_object_path );
      WrapperBase const * target = target_object->getWrapperBase( m_field_name );
      return target->getBufferedIOMetadata( num_indices );
    }

    inline localIndex CountLocalSetIndices( Group * problem_group ) const
    {
      return CountLocalSetIndicesExclusive( problem_group, m_set_names.size( ) );
    }

    localIndex CountLocalSetIndicesExclusive( Group * problem_group, localIndex set_idx = 0 ) const
    {
      GEOSX_ERROR_IF( m_set_names.size() == 0 && set_idx > 0, "No set names specified in input, but trying to sum local set indices." );
      ProblemManager const * pm = dynamicCast< ProblemManager const * >( problem_group );
      DomainPartition const & domain = *pm->getDomainPartition( );
      MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
      Group const * target_object = meshLevel.GetGroupByPath( m_object_path );
      WrapperBase const * target = target_object->getWrapperBase( m_field_name );
      GEOSX_ERROR_IF( ! target->isPackable(), "The object targeted for collection must be packable!" );
      localIndex num_indices = -1;
      if ( m_set_names.size( ) > 0 )
      {
        num_indices = 0;
        Group const * set_group = target_object->GetGroup( ObjectManagerBase::groupKeyStruct::setsString );
        localIndex idx = 0;
        for( auto & set_name : m_set_names )
        {
          if ( idx >= set_idx ) break;
          dataRepository::Wrapper< SortedArray< localIndex > > const * const set_wrapper = set_group->getWrapper< SortedArray< localIndex > >( set_name );
          if( set_wrapper != nullptr )
          {
            SortedArrayView< localIndex const > const & set = set_wrapper->reference();
            num_indices += set.size( );
          }
          ++idx;
        }
      }
      return num_indices;
    }

    virtual localIndex GetNumMetaCollectors( ) const override
    {
      return m_set_names.size( );
    }

    virtual std::unique_ptr<HistoryCollection> GetMetaCollector( Group * problem_group, localIndex meta_idx ) const override
    {
      localIndex set_idx_off = CountLocalSetIndicesExclusive( problem_group, meta_idx );
      return std::unique_ptr<HistoryCollection>(new SetMetaCollection(m_object_path,m_set_names[meta_idx],set_idx_off));
    }

    virtual void Collect( Group * domain_group,
                          real64 const GEOSX_UNUSED_PARAM(time_n),
                          real64 const GEOSX_UNUSED_PARAM(dt),
                          buffer_unit_type*& buffer ) override
    {
      DomainPartition const & domain = dynamicCast< DomainPartition const &>( *domain_group );
      MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
      Group const * target_object = meshLevel.GetGroupByPath( m_object_path );
      WrapperBase const * target = target_object->getWrapperBase( m_field_name );
      if ( m_set_names.size( ) > 0 )
      {
        Group const * set_group = target_object->GetGroup( ObjectManagerBase::groupKeyStruct::setsString );
        for( auto & set_name : m_set_names )
        {
          dataRepository::Wrapper< SortedArray< localIndex > > const * const set_wrapper = set_group->getWrapper< SortedArray< localIndex > >( set_name );
          if( set_wrapper != nullptr )
          {
            SortedArrayView< localIndex const > const & set = set_wrapper->reference();
            localIndex sz = set.size( );
            if ( sz > 0 )
            {
              array1d< localIndex > set_idxs;
              set_idxs.insert( 0, set.begin( ), set.end( ) );
              target->PackByIndex( buffer, set_idxs, false, true );
            }
          }
        }
      }
      else
      {
        target->Pack( buffer, false, true );
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
