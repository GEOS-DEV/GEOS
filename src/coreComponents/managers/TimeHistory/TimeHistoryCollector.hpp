/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TimeHistoryCollector.hpp
 */

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

  /**
   * @class HistoryCollection
   *
   * A task class for serializing time history data into a buffer for later I/O.
   */
  class HistoryCollection : public TaskBase
  {
  public:
    /// @copydoc geosx::dataRepository::Group::Group(std::string const & name, Group * const parent)
    HistoryCollection( string const & name, Group * parent ) :
      TaskBase( name, parent ),
      m_buffer_call()
    {   }

    /**
     * @brief Get the metadata for what this collector collects.
     * @param domain The ProblemDomain cast to a group.
     * @retunr A HistoryMetadata object describing  the history data being collected by this collector.
     */
    virtual HistoryMetadata GetMetadata( Group * domain )
    {
      GEOSX_UNUSED_VAR( domain );
      return HistoryMetadata( "null", 0, std::type_index(typeid(nullptr)) );
    }

    /**
     * @brief Get the number of collectors of meta-information (set indices, etc) writing time-independent information during initialization.
     * @return The number of collectors of meta-information for this collector.
     */
    virtual localIndex GetNumMetaCollectors( ) const { return 0; }

    /**
     * @brief Get a pointer to a collector of meta-information for this collector.
     * @param problem_group The ProblemManager cast to a group.
     * @param meta_idx Which of the meta-info collectors to return. (see HistoryCollection::GetNumMetaCollectors( ) ).
     * @param meta_rank_offset The offset for this rank for the meta-info collector, used to number index metadata consistently across the simulation.
     * @return A unique pointer to the HistoryCollection object used for meta-info collection. Intented to fall out of scope and desctruct immediately
     *         after being used to perform output during simulation initialization.
     */
    virtual std::unique_ptr<HistoryCollection> GetMetaCollector( Group * problem_group, localIndex meta_idx, globalIndex meta_rank_offset )
    {
      GEOSX_UNUSED_VAR( problem_group );
      GEOSX_UNUSED_VAR( meta_idx );
      GEOSX_UNUSED_VAR( meta_rank_offset );
      return std::unique_ptr<HistoryCollection>(nullptr);
    }

    /**
     * @brief Collect history information into the provided buffer. Typically called from HistoryCollection::Execute .
     * @param domain The ProblemDomain cast to a group.
     * @param time_n The current simulation time.
     * @param dt The current simulation time delta.
     * @param buffer A properly-sized buffer to serialize history data into.
     */
    virtual void Collect( Group * domain, real64 const time_n, real64 const dt, buffer_unit_type *& buffer )
    {
      GEOSX_UNUSED_VAR( domain );
      GEOSX_UNUSED_VAR( time_n );
      GEOSX_UNUSED_VAR( dt );
      GEOSX_UNUSED_VAR( buffer );
    }

    /**
     * @brief Collects history data.
     * @copydoc EventBase::Execute()
     */
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
      if ( rank == 0 && m_time_buffer_call )
      {
        buffer_unit_type * time_buffer = m_time_buffer_call();
        memcpy( time_buffer, &time_n, sizeof(decltype(time_n)) );
      }
    }

    /**
     * @brief Register a callback that gives the current head of the time history data buffer.
     * @param buffer_call A functional that when invoked returns a pointer to the head of a buffer at least large enough to
     *                    serialize one timestep of history data into.
     * @note This is typically meant to callback to BufferedHistoryIO::GetBufferHead( )
     */
    void RegisterBufferCall( std::function<buffer_unit_type*()> buffer_call )
    {
      m_buffer_call = buffer_call;
    }

    /**
     * @brief Get a metadata object relating the the Time variable itself.
     * @param A HistroyMetadata object describing the Time variable.
     */
    HistoryMetadata GetTimeMetadata( ) const
    {
      return HistoryMetadata("Time",1,std::type_index(typeid(real64)));
    }

    /**
     * @brief Register a callback that gives the current head of the time data buffer.
     * @param buffer_call A functional that when invoked returns a pointer to the head of a buffer at least large enough to
     *                    serialize one instance of the Time variable into.
     * @note This is typically meant to callback to BufferedHistoryIO::GetBufferHead( )
     */
    void RegisterTimeBufferCall( std::function<buffer_unit_type*()> time_buffer_call )
    {
      m_time_buffer_call = time_buffer_call;
    }

  private:
    std::function<buffer_unit_type*()> m_time_buffer_call;
    std::function<buffer_unit_type*()> m_buffer_call;
  };

  /**
   * @class SetMetaCollection
   *
   * A task class for serializing set collection-index information (indices into HistoryOutput related to a Set,
   *   not the simulation indices of the Set).
   */
  class SetMetaCollection : public HistoryCollection
  {
  public:
    /**
     * @brief Constructor
     * @param set_index_offset An offset constituting all 'prior' rank-sets + local-set index counts to correctly provide history output set indices.
     * @copydoc geosx::dataRepository::Group::Group(std::string const & name, Group * const parent)
     */
    SetMetaCollection( string const & object_path, string const & set_name, globalIndex set_index_offset ) :
      HistoryCollection( "SetMetaCollection", nullptr ),
      m_object_path( object_path ),
      m_set_name( set_name ),
      m_set_index_offset( set_index_offset )
    { }

    /// @copydoc HistoryCollection::GetMetadata
    virtual HistoryMetadata GetMetadata( Group * problem_group ) override
    {
      ProblemManager const * pm = dynamicCast< ProblemManager const * >( problem_group );
      DomainPartition const & domain = *pm->getDomainPartition( );
      MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
      Group const * target_object = meshLevel.GetGroupByPath( m_object_path );
      Group const * set_group = target_object->GetGroup( ObjectManagerBase::groupKeyStruct::setsString );
      dataRepository::Wrapper< SortedArray< localIndex > > const * const set_wrapper = set_group->getWrapper< SortedArray< localIndex > >( m_set_name );
      HistoryMetadata meta = set_wrapper->getBufferedIOMetadata( );
      string outname = meta.getName( );
      meta.setName( outname + " Indices" );
      meta.setType( std::type_index(typeid(globalIndex)) );
      return meta;
    }

    /// @copydoc HistoryCollection::Collect
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
        globalIndex sz = set.size( );
        if ( sz != 0 )
        {
	  // this might wind up having to move the meta_idx array back and forth to device a lot
          array1d< globalIndex > meta_idx( sz );
          for( localIndex ii = 0; ii < sz; ++ii )
          {
            meta_idx[ii] = m_set_index_offset + LvArray::integerConversion< globalIndex >( ii );
          }
          bufferOps::PackDataDevice< true >( buffer, meta_idx.toViewConst() );
        }
      }
    }
  protected:
    string m_object_path;
    string m_set_name;
    globalIndex m_set_index_offset;
  };


  /**
   * @class PackCollection
   *
   * A task class for serializing history information using the MPI communication packing routines.
   */
  class PackCollection : public HistoryCollection
  {
    public:
    PackCollection ( string const & name, Group * parent )
      : HistoryCollection( name, parent )
      , m_sets_indices( )
      , m_object_path( )
      , m_field_name( )
      , m_set_names( )
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

    /**
     * @brief Catalog name interface
     * @return This type's catalog name
     */
    static string CatalogName() { return "PackCollection"; }

    void InitializePostSubGroups( Group * const group ) override
    {
      UpdateSetsIndices( group );
    }

    /// @copydoc HistoryCollection::GetMetadata
    virtual HistoryMetadata GetMetadata( Group * problem_group ) override
    {
      localIndex num_indices = CountLocalSetIndices( problem_group );
      ProblemManager const * pm = dynamicCast< ProblemManager const * >( problem_group );
      DomainPartition const & domain = *pm->getDomainPartition( );
      MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
      Group const * target_object = meshLevel.GetGroupByPath( m_object_path );
      WrapperBase const * target = target_object->getWrapperBase( m_field_name );
      return target->getBufferedIOMetadata( num_indices );
    }

    void UpdateSetsIndices( Group * problem_group )
    {
      ProblemManager const * pm = dynamicCast< ProblemManager const * >( problem_group );
      DomainPartition const & domain = *pm->getDomainPartition( );
      MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
      Group const * target_object = meshLevel.GetGroupByPath( m_object_path );
      WrapperBase const * target = target_object->getWrapperBase( m_field_name );
      GEOSX_ERROR_IF( ! target->isPackable( false ), "The object targeted for collection must be packable!" );
      localIndex num_sets = m_set_names.size( );
      if ( num_sets > 0 )
      {
        Group const * set_group = target_object->GetGroup( ObjectManagerBase::groupKeyStruct::setsString );
	m_sets_indices.resize( num_sets );
	localIndex set_idx = 0;
        for( auto & set_name : m_set_names )
        {
          dataRepository::Wrapper< SortedArray< localIndex > > const * const set_wrapper = set_group->getWrapper< SortedArray< localIndex > >( set_name );
          if( set_wrapper != nullptr )
          {
            SortedArrayView< localIndex const > const & set = set_wrapper->reference();
	    m_sets_indices[ set_idx ].insert( 0, set.begin(), set.end() );
          }
	  set_idx++;
	}
      }
    }

    /**
     * @brief Count the total number of indices being collected by this process with this collector.
     * @param problem_group The ProblemManager cast to a Group.
     * @return The number of local indices being collected.
     */
    inline localIndex CountLocalSetIndices( Group * problem_group )
    {
      return CountLocalSetIndicesExclusive( problem_group, m_set_names.size( ) );
    }

    /**
     * @brief Count the number of indices being collected by this process for all sets up to and
     *        excluding the specified set index (see HistoryCollection::GetNumMetaCollectors).
     * @param problem_group The ProblemManager cast to a Group.
     * @param set_idx The index of the Set to count all other local indices prior to.
     * @return The number of indices associate with all Sets locally prior to the specified Set.
     */
    localIndex CountLocalSetIndicesExclusive( Group * problem_group, localIndex last_set_idx = 0 )
    {
      GEOSX_ERROR_IF( m_set_names.size() == 0 && last_set_idx > 0, "No set names specified in input, but trying to sum local set indices." );
      int num_sets = m_sets_indices.size( );
      if ( num_sets == 0 )
      {
	UpdateSetsIndices( problem_group );
      }
      localIndex num_indices = 0;
      for ( localIndex set_idx = 0; set_idx < last_set_idx; ++set_idx )
      {
	num_indices += m_sets_indices[ set_idx ].size( );
      } 
      return num_indices;
    }

    /// @copydoc HistoryCollection::GetNumMetaCollectors
    virtual localIndex GetNumMetaCollectors( ) const override
    {
      return m_set_names.size( );
    }

    /// @copydoc HistoryCollection::GetMetaCollector
    virtual std::unique_ptr<HistoryCollection> GetMetaCollector( Group * problem_group, localIndex meta_idx, globalIndex meta_rank_offset ) override
    {
      globalIndex local_set_offset = LvArray::integerConversion<globalIndex>( CountLocalSetIndicesExclusive( problem_group, meta_idx ) );
      return std::unique_ptr<HistoryCollection>(new SetMetaCollection(m_object_path,m_set_names[meta_idx],meta_rank_offset + local_set_offset));
    }

    /// @copydoc HistoryCollection::Collect
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
	localIndex set_idx = 0;
        for( auto & set_name : m_set_names )
        {
          dataRepository::Wrapper< SortedArray< localIndex > > const * const set_wrapper = set_group->getWrapper< SortedArray< localIndex > >( set_name );
          if( set_wrapper != nullptr )
          {
            SortedArrayView< localIndex const > const & set = set_wrapper->reference();
            localIndex sz = set.size( );
            if ( sz > 0 )
            {
	      // if we could directly transfer a sorted array to an array1d including on device this wouldn't require data movement
              target->PackByIndex( buffer, m_sets_indices[ set_idx ], false, true );
            }
          }
	  set_idx++;
        }
      }
      else
      {
        target->Pack( buffer, false, true );
      }
    }

    /// @cond DO_NOT_DOCUMENT
    struct viewKeysStruct
    {
      static constexpr auto objectPath = "objectPath";
      static constexpr auto fieldName = "fieldName";
      static constexpr auto setNames = "setNames";
    } keys;
    /// @endcond

  protected:
    // todo : replace this with a vector of references to the actual set sortedarrays (after packing rework to allow sorted arrays to be used for indexing)
    std::vector< array1d< localIndex > > m_sets_indices;
    string m_object_path;
    string m_field_name;
    string_array m_set_names;
  };
}

#endif
