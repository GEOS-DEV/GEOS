/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file ObjectManagerBase.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_OBJECTMANAGERBASE_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_OBJECTMANAGERBASE_HPP_

#include "dataRepository/ManagedGroup.hpp"

namespace geosx
{
class SiloFile;

/**
 * @class ObjectManagerBase
 * @brief The ObjectManagerBase is the base object of all object managers in the mesh data hierachy.
 *
 */
class ObjectManagerBase : public dataRepository::ManagedGroup
{
public:
  ObjectManagerBase() = delete;

  explicit ObjectManagerBase( std::string const & name,
                              dataRepository::ManagedGroup * const parent );

  ~ObjectManagerBase() override;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  using CatalogInterface = cxx_utilities::CatalogInterface< ObjectManagerBase, std::string const &, dataRepository::ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

  virtual const string getCatalogName() const = 0;
  ///@}

  using dataRepository::ManagedGroup::PackSize;
  using dataRepository::ManagedGroup::Pack;

  virtual localIndex PackSize( string_array const & wrapperNames,
                               arrayView1d<localIndex const> const & packList,
                               integer const recursive ) const override;


  virtual localIndex Pack( buffer_unit_type * & buffer,
                           string_array const & wrapperNames,
                           arrayView1d<localIndex const> const & packList,
                           integer const recursive )  const override;

  virtual localIndex Unpack( buffer_unit_type const *& buffer,
                             arrayView1d<localIndex> & packList,
                             integer const recursive )  override;

  template< bool DOPACK >
  localIndex PackSets( buffer_unit_type * & buffer,
                       arrayView1d<localIndex const> const & packList ) const;

  localIndex UnpackSets( buffer_unit_type const *& buffer );

  virtual void ViewPackingExclusionList( set<localIndex> & exclusionList ) const;


  virtual localIndex PackGlobalMapsSize( arrayView1d<localIndex> const & packList,
                                         integer const recursive ) const;

  virtual localIndex PackGlobalMaps( buffer_unit_type * & buffer,
                                     arrayView1d<localIndex> const & packList,
                                     integer const recursive ) const;

  void SetReceiveLists(  );



  virtual localIndex PackUpDownMapsSize( arrayView1d<localIndex const> const & packList ) const
  { return 0; }

  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d<localIndex const> const & packList ) const
  { return 0;}


  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                       array1d<localIndex> & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps )
  { return 0;}



  virtual localIndex UnpackGlobalMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       integer const recursive );

private:
  template< bool DOPACK >
  localIndex PackPrivate( buffer_unit_type * & buffer,
                          string_array const & wrapperNames,
                          arrayView1d<localIndex const> const & packList,
                          integer const recursive ) const;

  template< bool DOPACK >
  localIndex PackGlobalMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d<localIndex const> const & packList,
                                    integer const recursive ) const;

  //**********************************************************************************************************************
  // functions for compatibility with old data structure
  // TODO Deprecate or modernize all these suckers

public:

  using ObjectType = string;
  localIndex resize( localIndex const newSize,
                     const bool /*assignGlobals*/ )
  {
    dataRepository::ManagedGroup::resize(newSize);
    return 0;
  }

  using dataRepository::ManagedGroup::resize;

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
  void ConstructSetFromSetAndMap( const set<localIndex>& inputSet,
                                  const array2d<localIndex>& map,
                                  const std::string& newSetName );

  /// builds a new set on this object given another objects set and the map
  // between them
  void ConstructSetFromSetAndMap( const set<localIndex>& inputSet,
                                  const array1d<localIndex_array>& map,
                                  const std::string& newSetName );

  void ConstructGlobalToLocalMap();

  void ConstructLocalListOfBoundaryObjects( localIndex_array & objectList ) const;
  void ConstructGlobalListOfBoundaryObjects( globalIndex_array & objectList ) const;

  virtual void ExtractMapFromObjectForAssignGlobalIndexNumbers( ObjectManagerBase const * const ,
                                                                array1d<globalIndex_array>&  )
  {}

  void SetGhostRankForSenders( arrayView1d<localIndex> const & indicesToSend )
  {
    for( auto index : indicesToSend )
    {
//      GEOS_ERROR_IF( m_ghostRank[index] >= 0, "trying to set ghostRank of non-locally owned index: m_ghostRank["<<index<<"]="<<m_ghostRank[index] );
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
                             map< localIndex, set<globalIndex> > & unmappedIndices,
                             bool const clearIfUnmapped  );

  static void CleanUpMap( std::set<localIndex> const & targetIndices,
                          array1d<set<localIndex> > & upmap,
                          array2d<localIndex> const & downmap );

  static void CleanUpMap( std::set<localIndex> const & targetIndices,
                          array1d<set<localIndex> > & upmap,
                          array1d< array1d<localIndex > > const & downmap );

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



  ManagedGroup * sets()             {return &m_sets;}
  ManagedGroup const * sets() const {return &m_sets;}

  set<localIndex> & externalSet()
  {return m_sets.getReference<set<localIndex>>(m_ObjectManagerBaseViewKeys.externalSet);}

  set<localIndex> const & externalSet() const
  {return m_sets.getReference<set<localIndex>>(m_ObjectManagerBaseViewKeys.externalSet);}

  integer_array & isExternal()
  { return this->m_isExternal; }

  integer_array const & isExternal() const
  { return this->m_isExternal; }

  integer_array & GhostRank()
  { return this->m_ghostRank; }

  integer_array const & GhostRank() const
  { return this->m_ghostRank; }

  ManagedGroup m_sets;

  globalIndex_array  m_localToGlobalMap;
  map<globalIndex,localIndex>  m_globalToLocalMap;
  integer_array m_isExternal;
  integer_array m_ghostRank;

  real64 m_overAllocationFactor = 1.1;

  globalIndex m_maxGlobalIndex = -1;

//  localIndex_array m_ghostToSend;
 // localIndex_array m_ghostToReceive;

};


//template< typename T >
//void ObjectManagerBase::FixUpDownMaps()
//{
//
//}


template< typename TYPE_RELATION >
void ObjectManagerBase::FixUpDownMaps( TYPE_RELATION & relation,
                                       map< localIndex, array1d<globalIndex> > & unmappedIndices,
                                       bool const  )
{
  bool allValuesMapped = true;
  map<globalIndex,localIndex> const & globalToLocal = relation.RelatedObjectGlobalToLocal();
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
      GEOS_ERROR_IF( relation[li][a]==unmappedLocalIndexValue, "Index not set");
    }
  }
  GEOS_ERROR_IF( !allValuesMapped, "some values of unmappedIndices were not used");
  unmappedIndices.clear();
}


template< typename TYPE_RELATION >
void ObjectManagerBase::FixUpDownMaps( TYPE_RELATION & relation,
                                       map< localIndex, set<globalIndex> > & unmappedIndices,
                                       bool const clearIfUnmapped )
{
  map<globalIndex,localIndex> const & globalToLocal = relation.RelatedObjectGlobalToLocal();
  for( map< localIndex, set<globalIndex> >::iterator iter = unmappedIndices.begin() ;
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
      set<globalIndex> const & globalIndices = iter->second;
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



} /* namespace geosx */



typedef geosx::ObjectManagerBase ObjectDataStructureBaseT;

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_OBJECTMANAGERBASE_HPP_ */
