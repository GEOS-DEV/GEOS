#include "PackCollection.hpp"
#include "managers/GeosxState.hpp"

namespace geosx
{
PackCollection::PackCollection ( string const & name, Group * parent )
  : HistoryCollection( name, parent )
  , m_setsIndices( )
  , m_objectPath( )
  , m_fieldName( )
  , m_setNames( )
  , m_minimumSetSize( )
{
  registerWrapper( PackCollection::viewKeysStruct::objectPath, &m_objectPath ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The name of the object from which to retrieve field values." );

  registerWrapper( PackCollection::viewKeysStruct::fieldName, &m_fieldName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The name of the (packable) field associated with the specified object to retrieve data from" );

  registerWrapper( PackCollection::viewKeysStruct::setNames, &m_setNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The set(s) for which to retrieve data." );

  registerWrapper( PackCollection::viewKeysStruct::minSetSize, &m_minimumSetSize ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( -1 ).
    setDescription( "The minimum size of the set(s) to be collected (use for sets that expand during the simulation)." );
}

void PackCollection::initializePostSubGroups()
{
  localIndex numSets = m_setNames.size( );
  m_collectionCount = numSets == 0 ? 1 : numSets;
  DomainPartition & domain = getGlobalState().getProblemManager().getDomainPartition();
  updateSetsIndices( domain );
  HistoryCollection::initializePostSubGroups();
}

HistoryMetadata PackCollection::getMetadata( ProblemManager & pm, localIndex collectionIdx )
{
  DomainPartition & domain = pm.getDomainPartition();
  Group const * target_object = this->getTargetObject( domain );
  WrapperBase const & target = target_object->getWrapperBase( m_fieldName );
  if( m_setNames.size() != 0 )
  {
    GEOSX_ERROR_IF( collectionIdx >= m_setNames.size(), "Invalid collection index specified." );
    localIndex num_indices = m_setsIndices[ collectionIdx ].size( );
    num_indices = m_minimumSetSize > num_indices ? m_minimumSetSize : num_indices;
    HistoryMetadata metadata = target.getHistoryMetadata( num_indices );
    metadata.setName( metadata.getName() + " " + m_setNames[ collectionIdx ] );
    return metadata;
  }
  else
  {
    localIndex num_indices = m_setsIndices[ 0 ].size( );
    num_indices = m_minimumSetSize > num_indices ? m_minimumSetSize : num_indices;
    return target.getHistoryMetadata( num_indices );
  }
}

void PackCollection::updateSetsIndices( DomainPartition & domain )
{
  ObjectManagerBase const * target_object = this->getTargetObject( domain );
  WrapperBase const & target = target_object->getWrapperBase( m_fieldName );
  GEOSX_ERROR_IF( !target.isPackable( false ), "The object targeted for collection must be packable!" );
  localIndex num_sets = m_setNames.size( );

  if( num_sets > 0 )
  {
    // if sets are specified we retrieve the field only from those sets

    Group const & setGroup = target_object->getGroup( ObjectManagerBase::groupKeyStruct::setsString() );
    m_setsIndices.resize( num_sets );
    localIndex set_idx = 0;
    for( auto & set_name : m_setNames )
    {
      if( setGroup.hasWrapper( set_name ) )
      {
        SortedArrayView< localIndex const > const & set = setGroup.getReference< SortedArray< localIndex > >( set_name );
        m_setsIndices[ set_idx ].resize( 0 );
        if( set.size() > 0 )
        {
          m_setsIndices[ set_idx ].insert( 0, set.begin(), set.end() );
        }
      }
      set_idx++;
    }
  }
  else
  {
    // if no set is specified we retrieve the entire field
    m_setsIndices.resize( 1 );
    m_setsIndices[0].resize( target_object->size());
    for( localIndex k=0; k <  target_object->size(); k++ )
    {
      m_setsIndices[0][k] = k;
    }
  }
}

void PackCollection::filterGhostIndices( localIndex const setIndex,
                                         array1d< localIndex > & set,
                                         arrayView1d< integer const > const & ghostRank )
{
  // Resize the indices array

  // 3. fill in the non ghost indices
  localIndex idx = 0;
  for( localIndex k=0; k < m_setsIndices[setIndex].size(); k++ )
  {
    if( ghostRank[m_setsIndices[setIndex][k]] < 0 )
    {
      set[idx] = m_setsIndices[setIndex][k];
      idx++;
    }
  }
}

// TODO : once we add additional history collectors, this should likely be pulled into a super-class
ObjectManagerBase const * PackCollection::getTargetObject( DomainPartition & domain )
{
  dataRepository::Group * targetGroup = &domain.getMeshBody( 0 ).getMeshLevel( 0 );
  string_array const targetTokens = stringutilities::Tokenize( m_objectPath, "/" );
  localIndex const targetTokenLength = LvArray::integerConversion< localIndex >( targetTokens.size() );

  for( localIndex pathLevel = 0; pathLevel < targetTokenLength; ++pathLevel )
  {
    dataRepository::Group * const elemRegionSubGroup = targetGroup->getGroupPointer( ElementRegionManager::groupKeyStruct::elementRegionsGroup() );
    if( elemRegionSubGroup != nullptr )
    {
      targetGroup = elemRegionSubGroup;
    }
    dataRepository::Group * const elemSubRegionSubGroup = targetGroup->getGroupPointer( ElementRegionBase::viewKeyStruct::elementSubRegions() );
    if( elemSubRegionSubGroup != nullptr )
    {
      targetGroup = elemSubRegionSubGroup;
    }
    if( targetTokens[pathLevel] == ElementRegionManager::groupKeyStruct::elementRegionsGroup() ||
        targetTokens[pathLevel] == ElementRegionBase::viewKeyStruct::elementSubRegions() )
    {
      continue;
    }

    targetGroup = &targetGroup->getGroup( targetTokens[pathLevel] );
  }

  return dynamicCast< ObjectManagerBase const * >( targetGroup );
}

void PackCollection::collect( DomainPartition & domain,
                              real64 const GEOSX_UNUSED_PARAM( time_n ),
                              real64 const GEOSX_UNUSED_PARAM( dt ),
                              localIndex collectionIdx,
                              buffer_unit_type * & buffer )
{
  GEOSX_MARK_FUNCTION;
  GEOSX_ERROR_IF( collectionIdx >= getCollectionCount( ), "Attempting to collection from an invalid collection index!" );
  ObjectManagerBase const * target_object = this->getTargetObject( domain );
  WrapperBase const & target = target_object->getWrapperBase( m_fieldName );

  arrayView1d< integer const > const ghostRank = target_object->ghostRank();
  localIndex numIndices = 0;

  //  count non ghost indices
  for( localIndex k=0; k <  m_setsIndices[collectionIdx].size(); k++ )
  {
    if( ghostRank[m_setsIndices[collectionIdx][k]] < 0 )
    {
      numIndices++;
    }
  }

  if( numIndices > 0 )
  {
    array1d< localIndex > setIndices( numIndices );
    filterGhostIndices( collectionIdx, setIndices, ghostRank );
    // if we could directly transfer a sorted array to an array1d including on device this wouldn't require storing a copy of the indices
    target.packByIndex( buffer, setIndices, false, true );
  }

}

REGISTER_CATALOG_ENTRY( TaskBase, PackCollection, string const &, Group * const )
}
