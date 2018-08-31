/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
#include "DocumentationNode.hpp"

namespace geosx
{
class SiloFile;
namespace dataRepository
{
namespace keys
{
string const sets("Sets");
}
}


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

//  explicit ObjectManagerBase( std::string const & name,
//                              dataRepository::ManagedGroup * const parent,
//                              cxx_utilities::DocumentationNode * docNode );

  ~ObjectManagerBase() override;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  using CatalogInterface = cxx_utilities::CatalogInterface< ObjectManagerBase, std::string const &, dataRepository::ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

  virtual const string getCatalogName() const = 0;
  ///@}

  virtual void FillDocumentationNode() override;

  virtual void InitializePostSubGroups( ManagedGroup * const ) override;

  using dataRepository::ManagedGroup::PackSize;
  using dataRepository::ManagedGroup::Pack;

  virtual localIndex PackSize( string_array const & wrapperNames,
                        localIndex_array const & packList,
                        integer const recursive ) const override;


  virtual localIndex Pack( buffer_unit_type * & buffer,
                    string_array const & wrapperNames,
                    localIndex_array const & packList,
                    integer const recursive )  const override;

//  virtual int Unpack( buffer_unit_type const *& buffer,
//                      integer const recursive )  override;

  virtual localIndex Unpack( buffer_unit_type const *& buffer,
                      localIndex_array & packList,
                      integer const recursive )  override;

  virtual void ViewPackingExclusionList( set<localIndex> & exclusionList ) const;


  virtual localIndex PackGlobalMapsSize( localIndex_array const & packList,
                                  integer const recursive ) const;

  virtual localIndex PackGlobalMaps( buffer_unit_type * & buffer,
                              localIndex_array const & packList,
                              integer const recursive ) const;

  void SetReceiveLists(  );



  virtual localIndex PackUpDownMapsSize( localIndex_array const & packList ) const
  { return 0; }

  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                              localIndex_array const & packList ) const
  { return 0;}


  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                localIndex_array const & packList )
  { return 0;}



  virtual localIndex UnpackGlobalMaps( buffer_unit_type const * & buffer,
                                localIndex_array & packList,
                                integer const recursive );

private:
  template< bool DOPACK >
  localIndex PackPrivate( buffer_unit_type * & buffer,
                   string_array const & wrapperNames,
                   localIndex_array const & packList,
                   integer const recursive ) const;

  template< bool DOPACK >
  localIndex PackGlobalMapsPrivate( buffer_unit_type * & buffer,
                             localIndex_array const & packList,
                             integer const recursive ) const;

//  template< bool DOPACK >
//  int UnpackPrivate( buffer_unit_type const *& buffer,
//                     localIndex_array const & packList,
//                     integer const recursive );


//  cxx_utilities::DocumentationNode * m_docNode;



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

//    localIndex m_DataLengths;
//
//    localIndex DataLengths() const { return size(); }




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



  /// returns reference to specified field
  template< FieldKey FIELDKEY>
  typename dataRepository::ViewWrapper< array1d< typename Field<FIELDKEY>::Type > >::rtype GetFieldData( )
  {
    return const_cast<typename dataRepository::ViewWrapper< array1d< typename Field<FIELDKEY>::Type > >::rtype>( static_cast<const ObjectManagerBase&>(*this).
                                                                                                               GetFieldData<FIELDKEY>());
  }


  /// returns const reference to specified field
  template< FieldKey FIELDKEY>
  typename dataRepository::ViewWrapper< array1d< typename Field<FIELDKEY>::Type > >::rtype_const GetFieldData( ) const
  {
    return this->getData< array1d< typename Field<FIELDKEY>::Type > >( string(Field<FIELDKEY>::Name()) );
  }


  /// returns reference to specified field
  template< typename TYPE >
  typename dataRepository::ViewWrapper< array1d< TYPE > >::rtype GetFieldData( const std::string& fieldName )
  {
    return const_cast<typename dataRepository::ViewWrapper<array1d<TYPE> >::rtype>( static_cast<const ObjectManagerBase&>(*this).GetFieldData<TYPE>(fieldName));
  }

  /// returns const reference to specified field
  template< typename TYPE >
  typename dataRepository::ViewWrapper< array1d< TYPE > >::rtype_const GetFieldData( const std::string& name ) const
  {
    return this->getData< array1d<TYPE> >( name );
  }



  /// returns reference to specified field
  template< FieldKey FIELDKEY>
  typename Field<FIELDKEY>::Type * GetFieldDataPointer( )
  {
    return &this->getReference< typename Field<FIELDKEY>::Type >( Field<FIELDKEY>::Name() );
  }


  /// returns const reference to specified field
  template< FieldKey FIELDKEY>
  typename Field<FIELDKEY>::Type const * GetFieldDataPointer( ) const
  {
    return &this->getReference< typename Field<FIELDKEY>::Type >( Field<FIELDKEY>::Name() );
  }

  /// returns reference to specified field
  template< typename TYPE >
  TYPE * GetFieldDataPointer( const std::string& fieldName )
  {
    return &this->getReference< TYPE >( fieldName );
  }

  /// returns const reference to specified field
  template< typename TYPE >
  TYPE const * GetFieldDataPointer( const std::string& fieldName ) const
  {
    return &this->getReference< TYPE >( fieldName );
  }



  /// add a data field to a member
  template< typename T >
  int AddKeylessDataField( const std::string& name, const bool restart = false, const bool plot = false )
  {
    this->RegisterViewWrapper<T>(name);
    (void)restart;
    (void)plot;
    return 0;
  }


  /// add a data field to a member
  template< FieldKey FIELDKEY >
  int AddKeyedDataField()
  {
    this->RegisterViewWrapper<typename Field<FIELDKEY>::Type>( Field<FIELDKEY>::Name() );
    return 0;
  }





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

  virtual void ExtractMapFromObjectForAssignGlobalIndexNumbers( ObjectManagerBase const & ,
                                                                array1d<globalIndex_array>&  )
  {

  }

  localIndex GetNumberOfGhosts() const;

  localIndex GetNumberOfLocalIndices() const;

  integer SplitObject( localIndex const indexToSplit,
                       int const rank,
                       localIndex & newIndex );

  void CopyObject( localIndex const source, localIndex const destination );



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

//  localIndex_array m_ghostToSend;
 // localIndex_array m_ghostToReceive;


};



} /* namespace geosx */

typedef geosx::ObjectManagerBase ObjectDataStructureBaseT;

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_OBJECTMANAGERBASE_HPP_ */
