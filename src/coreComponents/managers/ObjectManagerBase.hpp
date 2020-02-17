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
 * @file ObjectManagerBase.hpp
 */

#ifndef GEOSX_MANAGERS_OBJECTMANAGERBASE_HPP_
#define GEOSX_MANAGERS_OBJECTMANAGERBASE_HPP_

#include "dataRepository/Group.hpp"
#include "common/TimingMacros.hpp"

namespace geosx
{
class SiloFile;

/**
 * @class ObjectManagerBase
 * @brief The ObjectManagerBase is the base object of all object managers in the mesh data hierachy.
 *
 */
class ObjectManagerBase : public dataRepository::Group
{
public:
  ObjectManagerBase() = delete;

  explicit ObjectManagerBase( std::string const & name,
                              dataRepository::Group * const parent );

  ~ObjectManagerBase() override;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  using CatalogInterface = dataRepository::CatalogInterface< ObjectManagerBase, std::string const &, dataRepository::Group * const >;
  static CatalogInterface::CatalogType& GetCatalog();

  virtual const string getCatalogName() const = 0;
  ///@}

  using dataRepository::Group::PackSize;
  using dataRepository::Group::Pack;

  virtual localIndex PackSize( string_array const & wrapperNames,
                               arrayView1d<localIndex const> const & packList,
                               integer const recursive,
                               bool on_device = false ) const override;


  virtual localIndex Pack( buffer_unit_type * & buffer,
                           string_array const & wrapperNames,
                           arrayView1d<localIndex const> const & packList,
                           integer const recursive,
                           bool on_device = false)  const override;

  virtual localIndex Unpack( buffer_unit_type const *& buffer,
                             arrayView1d<localIndex> & packList,
                             integer const recursive,
                             bool on_device = false ) override;

  template< bool DOPACK >
  localIndex PackSets( buffer_unit_type * & buffer,
                       arrayView1d<localIndex const> const & packList ) const;

  localIndex UnpackSets( buffer_unit_type const *& buffer );

  virtual void ViewPackingExclusionList( SortedArray<localIndex> & exclusionList ) const;


  virtual localIndex PackGlobalMapsSize( arrayView1d<localIndex> const & packList,
                                         integer const recursive ) const;

  virtual localIndex PackGlobalMaps( buffer_unit_type * & buffer,
                                     arrayView1d<localIndex> const & packList,
                                     integer const recursive ) const;

  void SetReceiveLists(  );



  virtual localIndex PackUpDownMapsSize( arrayView1d<localIndex const> const & GEOSX_UNUSED_PARAM( packList ) ) const
  { return 0; }

  virtual localIndex PackUpDownMaps( buffer_unit_type * & GEOSX_UNUSED_PARAM( buffer ),
                                     arrayView1d<localIndex const> const & GEOSX_UNUSED_PARAM( packList ) ) const
  { return 0; }

  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & GEOSX_UNUSED_PARAM( buffer ),
                                       array1d<localIndex> & GEOSX_UNUSED_PARAM( packList ),
                                       bool const GEOSX_UNUSED_PARAM( overwriteUpMaps ),
                                       bool const GEOSX_UNUSED_PARAM( overwriteDownMaps ) )
  { return 0; }


  virtual localIndex UnpackGlobalMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       integer const recursive );


  localIndex PackParentChildMapsSize( arrayView1d<localIndex const> const & packList ) const
  {
    buffer_unit_type * buffer = nullptr;
    return PackParentChildMapsPrivate<false>( buffer, packList );
  }

  localIndex PackParentChildMaps( buffer_unit_type * & buffer,
                                  arrayView1d<localIndex const> const & packList ) const
  {
    return PackParentChildMapsPrivate<true>( buffer, packList );
  }

  localIndex UnpackParentChildMaps( buffer_unit_type const * & buffer,
                                    localIndex_array & packList );


private:
  template< bool DOPACK >
  localIndex PackPrivate( buffer_unit_type * & buffer,
                          string_array const & wrapperNames,
                          arrayView1d<localIndex const> const & packList,
                          integer const recursive,
                          bool on_device) const;

  template< bool DOPACK >
  localIndex PackGlobalMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d<localIndex const> const & packList,
                                    integer const recursive ) const;

  template<bool DOPACK>
  localIndex PackParentChildMapsPrivate( buffer_unit_type * & buffer,
                                         arrayView1d<localIndex const> const & packList ) const;


  //**********************************************************************************************************************
  // functions for compatibility with old data structure
  // TODO Deprecate or modernize all these suckers

public:

  using ObjectType = string;
  localIndex resize( localIndex const newSize,
                     const bool /*assignGlobals*/ )
  {
    dataRepository::Group::resize(newSize);
    return 0;
  }

  using dataRepository::Group::resize;

  void WriteSilo( SiloFile& siloFile,
                  const std::string& meshname,
                  const int centering,
                  const int cycleNum,
                  const realT problemTime,
                  const bool isRestart,
                  const std::string& multiRoot,
                  const std::string& regionName = "none",
                  const localIndex_array& mask = localIndex_array() ) const;


  void ReadSilo( const SiloFile& siloFile,
                 const std::string& meshname,
                 const int centering,
                 const int cycleNum,
                 const realT problemTime,
                 const bool isRestart,
                 const std::string& regionName = "none",
                 const localIndex_array& mask = localIndex_array() );

  void CreateSet( const std::string& newSetName );

  /// builds a new set on this object given another objects set and the map
  // between them
  void ConstructSetFromSetAndMap( SortedArrayView<localIndex const> const & inputSet,
                                  const array2d<localIndex>& map,
                                  const std::string& newSetName );

  /// builds a new set on this object given another objects set and the map
  // between them
  void ConstructSetFromSetAndMap( SortedArrayView<localIndex const> const & inputSet,
                                  const array1d<localIndex_array>& map,
                                  const std::string& newSetName );

  void ConstructSetFromSetAndMap( SortedArrayView<localIndex const> const & inputSet,
                                  ArrayOfArraysView< localIndex const > const & map,
                                  const std::string& setName );

  void ConstructGlobalToLocalMap();

  void ConstructLocalListOfBoundaryObjects( localIndex_array & objectList ) const;
  void ConstructGlobalListOfBoundaryObjects( globalIndex_array & objectList ) const;

  virtual void ExtractMapFromObjectForAssignGlobalIndexNumbers( ObjectManagerBase const * const,
                                                                std::vector< std::vector< globalIndex > > & )
  {}

  void SetGhostRankForSenders( arrayView1d<localIndex> const & indicesToSend )
  {
    for( auto index : indicesToSend )
    {
      GEOSX_ERROR_IF( m_ghostRank[index] >= 0,
                     "trying to set ghostRank of non-locally owned index: "
                     "m_ghostRank[" << index << "]=" << m_ghostRank[index] );
      m_ghostRank[index] = -1;
    }
  }

  localIndex GetNumberOfGhosts() const;

  localIndex GetNumberOfLocalIndices() const;

  integer SplitObject( localIndex const indexToSplit,
                       int const rank,
                       localIndex & newIndex );

  /**
   * @brief sets the value of m_ghostRank to the value of the objects parent.
   * @param indices the list of indices for which to set the ghost rank
   *
   * This function takes a list of indices, and then sets the value of m_ghostRank for those indices to be equal to the
   * value of the "parent" index. This assumes that "parentIndex" is allocated and filled correctly.
   */
  void inheritGhostRankFromParent( std::set<localIndex> const & indices );

  void CopyObject( localIndex const source, localIndex const destination );

  void SetMaxGlobalIndex();

  template< typename TYPE_RELATION >
  static void FixUpDownMaps( TYPE_RELATION & relation,
                             map< localIndex, array1d<globalIndex> > & unmappedIndices,
                             bool const clearIfUnmapped );

  template< typename TYPE_RELATION >
  static void FixUpDownMaps( TYPE_RELATION & relation,
                             map< localIndex, SortedArray<globalIndex> > & unmappedIndices,
                             bool const clearIfUnmapped );

  static void FixUpDownMaps( ArrayOfSets< localIndex > & relation,
                             unordered_map<globalIndex,localIndex> const & globalToLocal,
                             map< localIndex, SortedArray<globalIndex> > & unmappedIndices,
                             bool const clearIfUnmapped );

  static void CleanUpMap( std::set<localIndex> const & targetIndices,
                          array1d<SortedArray<localIndex> > & upmap,
                          arrayView2d<localIndex const> const & downmap );

  static void CleanUpMap( std::set<localIndex> const & targetIndices,
                          ArrayOfSetsView< localIndex > const & upmap,
                          arrayView2d< localIndex const > const & downmap );

  static void CleanUpMap( std::set<localIndex> const & targetIndices,
                          array1d<SortedArray<localIndex> > & upmap,
                          arrayView1d< arrayView1d<localIndex const > const > const & downmap );

  static void CleanUpMap( std::set<localIndex> const & targetIndices,
                          ArrayOfSetsView< localIndex > const & upmap,
                          arrayView1d< arrayView1d< localIndex const > const > const & downmap );

  static void CleanUpMap( std::set<localIndex> const & targetIndices,
                          ArrayOfSetsView< localIndex > const & upmap,
                          ArrayOfArraysView< localIndex const > const & downmap );

  virtual void enforceStateFieldConsistencyPostTopologyChange( std::set<localIndex> const & targetIndices );

  static localIndex GetParentRecusive( arraySlice1d<localIndex const> const & parentIndices,
                                       localIndex const lookup )
  {
    localIndex rval = lookup;

    while( parentIndices[rval] != -1 )
    {
      rval = parentIndices[rval];
    }

    return rval;
  }


  //**********************************************************************************************************************

  /**
   * @brief struct to serve as a container for variable strings and keys
   * @struct viewKeyStruct
   *
   */
  struct viewKeyStruct
  {

    static constexpr auto adjacencyListString = "adjacencyList";
    static constexpr auto childIndexString = "childIndex";
    static constexpr auto domainBoundaryIndicatorString = "domainBoundaryIndicator";
    static constexpr auto externalSetString = "externalSet";
    static constexpr auto ghostRankString = "ghostRank";
    static constexpr auto ghostsToSendString = "ghostsToSend";
    static constexpr auto ghostsToReceiveString = "ghostsToReceive";
    static constexpr auto globalToLocalMapString = "globalToLocalMap";
    static constexpr auto isExternalString = "isExternal";
    static constexpr auto localToGlobalMapString = "localToGlobalMap";
    static constexpr auto matchedPartitionBoundaryObjectsString = "matchedPartitionBoundaryObjects";
    static constexpr auto parentIndexString = "parentIndex";

    dataRepository::ViewKey adjacencyList = { adjacencyListString };
    dataRepository::ViewKey childIndex = { childIndexString };
    dataRepository::ViewKey domainBoundaryIndicator = { domainBoundaryIndicatorString };
    dataRepository::ViewKey externalSet = { externalSetString };
    dataRepository::ViewKey ghostRank = { ghostRankString };
    dataRepository::ViewKey ghostsToSend = { ghostsToSendString };
    dataRepository::ViewKey ghostsToReceive = { ghostsToReceiveString };
    dataRepository::ViewKey globalToLocalMap = { globalToLocalMapString };
    dataRepository::ViewKey isExternal = { isExternalString };
    dataRepository::ViewKey localToGlobalMap = { localToGlobalMapString };
    dataRepository::ViewKey matchedPartitionBoundaryObjects = { matchedPartitionBoundaryObjectsString };
    dataRepository::ViewKey parentIndex = { parentIndexString };
  } m_ObjectManagerBaseViewKeys;


  /**
   * @brief struct to serve as a container for group strings and keys
   * @struct viewKeyStruct
   *
   */
  struct groupKeyStruct
  {
    static constexpr auto setsString = "sets";
    static constexpr auto neighborDataString = "neighborData";
    dataRepository::GroupKey sets = { setsString };
    dataRepository::GroupKey neighborData = { neighborDataString };
  } m_ObjectManagerBaseGroupKeys;


  virtual viewKeyStruct & viewKeys() { return m_ObjectManagerBaseViewKeys; }
  virtual viewKeyStruct const & viewKeys() const { return m_ObjectManagerBaseViewKeys; }

  virtual groupKeyStruct & groupKeys() { return m_ObjectManagerBaseGroupKeys; }
  virtual groupKeyStruct const & groupKeys() const { return m_ObjectManagerBaseGroupKeys; }



  Group * sets()             {return &m_sets;}
  Group const * sets() const {return &m_sets;}

  SortedArray<localIndex> & externalSet()
  {return m_sets.getReference<SortedArray<localIndex>>(m_ObjectManagerBaseViewKeys.externalSet);}

  SortedArray<localIndex> const & externalSet() const
  {return m_sets.getReference<SortedArray<localIndex>>(m_ObjectManagerBaseViewKeys.externalSet);}

  integer_array & isExternal()
  { return this->m_isExternal; }

  integer_array const & isExternal() const
  { return this->m_isExternal; }

  integer_array & GhostRank()
  { return this->m_ghostRank; }

  integer_array const & GhostRank() const
  { return this->m_ghostRank; }

  Group m_sets;

  globalIndex_array  m_localToGlobalMap;
  unordered_map<globalIndex,localIndex>  m_globalToLocalMap;
  integer_array m_isExternal;
  integer_array m_ghostRank;

  real64 m_overAllocationFactor = 1.1;

  globalIndex m_maxGlobalIndex = -1;

//  localIndex_array m_ghostToSend;
 // localIndex_array m_ghostToReceive;

};


template< typename TYPE_RELATION >
void ObjectManagerBase::FixUpDownMaps( TYPE_RELATION & relation,
                                       map< localIndex, array1d<globalIndex> > & unmappedIndices,
                                       bool const  )
{
  GEOSX_MARK_FUNCTION;

  bool allValuesMapped = true;
  unordered_map<globalIndex,localIndex> const & globalToLocal = relation.RelatedObjectGlobalToLocal();
  for( map< localIndex, array1d<globalIndex> >::iterator iter = unmappedIndices.begin() ;
       iter != unmappedIndices.end() ;
       ++iter )
  {
    localIndex const li = iter->first;
    array1d<globalIndex> const & globalIndices = iter->second;
    for( localIndex a=0 ; a<globalIndices.size() ; ++a )
    {
      if( globalIndices[a] != unmappedLocalIndexValue )
      {
        if( relation[li][a] == unmappedLocalIndexValue  )
        {
          relation[li][a] = globalToLocal.at(globalIndices[a]);
        }
        else
        {
          allValuesMapped = false;
        }
      }
      GEOSX_ERROR_IF( relation[li][a]==unmappedLocalIndexValue, "Index not set");
    }
  }
  GEOSX_ERROR_IF( !allValuesMapped, "some values of unmappedIndices were not used");
  unmappedIndices.clear();
}


template< typename TYPE_RELATION >
void ObjectManagerBase::FixUpDownMaps( TYPE_RELATION & relation,
                                       map< localIndex, SortedArray<globalIndex> > & unmappedIndices,
                                       bool const clearIfUnmapped )
{
  GEOSX_MARK_FUNCTION;

  unordered_map<globalIndex,localIndex> const & globalToLocal = relation.RelatedObjectGlobalToLocal();
  for( map< localIndex, SortedArray<globalIndex> >::iterator iter = unmappedIndices.begin() ;
       iter != unmappedIndices.end() ;
       ++iter )
  {
    localIndex const li = iter->first;
    if( clearIfUnmapped )
    {
      relation[li].clear();
    }
    else
    {
      SortedArray<globalIndex> const & globalIndices = iter->second;
      for( auto const newGlobalIndex : globalIndices )
      {
        // NOTE: This simply ignores if newGlobalIndex is not found. This is OK if this function is
        // used for an upmap and the object shouldn't exist on this rank. There should be a better
        // way to check this.
        auto iterG2L = globalToLocal.find(newGlobalIndex);
        if( iterG2L != globalToLocal.end() )
        {
          relation[li].insert( iterG2L->second );
        }
      }
    }
  }
  unmappedIndices.clear();
}

inline
void ObjectManagerBase::FixUpDownMaps( ArrayOfSets< localIndex > & relation,
                                       unordered_map<globalIndex,localIndex> const & globalToLocal,
                                       map< localIndex, SortedArray<globalIndex> > & unmappedIndices,
                                       bool const clearIfUnmapped )
{
  GEOSX_MARK_FUNCTION;

  for( map< localIndex, SortedArray<globalIndex> >::iterator iter = unmappedIndices.begin() ;
       iter != unmappedIndices.end() ;
       ++iter )
  {
    localIndex const li = iter->first;
    if( clearIfUnmapped )
    {
      relation.clearSet( li );
    }
    else
    {
      SortedArray<globalIndex> const & globalIndices = iter->second;
      for( globalIndex const newGlobalIndex : globalIndices )
      {
        // NOTE: This simply ignores if newGlobalIndex is not found. This is OK if this function is
        // used for an upmap and the object shouldn't exist on this rank. There should be a better
        // way to check this.
        auto iterG2L = globalToLocal.find(newGlobalIndex);
        if( iterG2L != globalToLocal.end() )
        {
          relation.insertIntoSet( li, iterG2L->second );
        }
      }
    }
  }
  unmappedIndices.clear();
}

} /* namespace geosx */



typedef geosx::ObjectManagerBase ObjectDataStructureBaseT;

#endif /* GEOSX_MANAGERS_OBJECTMANAGERBASE_HPP_ */
