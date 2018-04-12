/*
 * ObjectManagerBase.hpp
 *
 *  Created on: Sep 15, 2016
 *      Author: settgast1
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



class ObjectManagerBase : public dataRepository::ManagedGroup
{
public:
  ObjectManagerBase() = delete;

  explicit ObjectManagerBase( std::string const & name,
                              dataRepository::ManagedGroup * const parent );

//  explicit ObjectManagerBase( std::string const & name,
//                              dataRepository::ManagedGroup * const parent,
//                              cxx_utilities::DocumentationNode * docNode );

  ~ObjectManagerBase();

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  using CatalogInterface = cxx_utilities::CatalogInterface< ObjectManagerBase, std::string const &, dataRepository::ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

  virtual const string getCatalogName() const = 0;
  ///@}

  virtual void FillDocumentationNode() override;


  using dataRepository::ManagedGroup::PackSize;
  using dataRepository::ManagedGroup::Pack;

  virtual int PackSize( array<string> const & wrapperNames,
                        localIndex_array const & packList,
                        integer const includeGlobalIndices,
                        integer const recursive ) const;


  virtual int Pack( buffer_unit_type * & buffer,
                    array<string> const & wrapperNames,
                    localIndex_array const & packList,
                    integer const includeGlobalIndices,
                    integer const recursive ) const;

//  virtual int Unpack( buffer_unit_type const *& buffer,
//                      integer const recursive )  override;

  virtual int Unpack( buffer_unit_type const *& buffer,
                      localIndex_array const & packList,
                      integer const recursive )  override;

  virtual void ViewPackingExclusionList( set<localIndex> & exclusionList ) const;


private:
  template< bool DOPACK >
  int PackPrivate( buffer_unit_type * & buffer,
                   array<string> const & wrapperNames,
                   localIndex_array const & packList,
                   integer const includeGlobalIndices,
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
  typename dataRepository::ViewWrapper< array< typename Field<FIELDKEY>::Type > >::rtype GetFieldData( )
  {
    return const_cast<typename dataRepository::ViewWrapper< array< typename Field<FIELDKEY>::Type > >::rtype>( static_cast<const ObjectManagerBase&>(*this).
                                                                                                               GetFieldData<FIELDKEY>());
  }


  /// returns const reference to specified field
  template< FieldKey FIELDKEY>
  typename dataRepository::ViewWrapper< array< typename Field<FIELDKEY>::Type > >::rtype_const GetFieldData( ) const
  {
    return this->getData< array< typename Field<FIELDKEY>::Type > >( string(Field<FIELDKEY>::Name()) );
  }


  /// returns reference to specified field
  template< typename TYPE >
  typename dataRepository::ViewWrapper< array< TYPE > >::rtype GetFieldData( const std::string& fieldName )
  {
    return const_cast<typename dataRepository::ViewWrapper<array<TYPE> >::rtype>( static_cast<const ObjectManagerBase&>(*this).GetFieldData<TYPE>(fieldName));
  }

  /// returns const reference to specified field
  template< typename TYPE >
  typename dataRepository::ViewWrapper< array< TYPE > >::rtype_const GetFieldData( const std::string& name ) const
  {
    return this->getData< array<TYPE> >( name );
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
  void ConstructSetFromSetAndMap( const lSet& inputSet,
                                  const lArray2d& map,
                                  const std::string& newSetName );

  /// builds a new set on this object given another objects set and the map
  // between them
  void ConstructSetFromSetAndMap( const lSet& inputSet,
                                  const array<localIndex_array>& map,
                                  const std::string& newSetName );

  void ConstructGlobalToLocalMap();

  void ConstructLocalListOfBoundaryObjects( localIndex_array & objectList ) const;
  void ConstructGlobalListOfBoundaryObjects( globalIndex_array & objectList ) const;

  virtual void ExtractMapFromObjectForAssignGlobalIndexNumbers( ObjectManagerBase const & ,
                                                                array<globalIndex_array>&  )
  {

  }

  //**********************************************************************************************************************

  struct viewKeyStruct
  {

    static constexpr auto domainBoundaryIndicatorString = "domainBoundaryIndicator";
    static constexpr auto ghostRankString = "ghostRank";
    static constexpr auto ghostsToSendString = "ghostsToSend";
    static constexpr auto globalToLocalMapString = "globalToLocalMap";
    static constexpr auto isExternalString = "isExternal";
    static constexpr auto localToGlobalMapString = "localToGlobalMap";
    static constexpr auto matchedPartitionBoundaryObjectsString = "matchedPartitionBoundaryObjects";

    dataRepository::ViewKey domainBoundaryIndicator = { domainBoundaryIndicatorString };
    dataRepository::ViewKey ghostRank = { ghostRankString };
    dataRepository::ViewKey ghostsToSend = { ghostsToSendString };
    dataRepository::ViewKey globalToLocalMap = { globalToLocalMapString };
    dataRepository::ViewKey isExternal = { isExternalString };
    dataRepository::ViewKey localToGlobalMap = { localToGlobalMapString };
    dataRepository::ViewKey matchedPartitionBoundaryObjects = { matchedPartitionBoundaryObjectsString };



  } viewKeys;

  struct groupKeyStruct
  {
    static constexpr auto setsString = "sets";
    static constexpr auto neighborDataString = "neighborData";
    dataRepository::GroupKey sets = { setsString };
    dataRepository::GroupKey neighborData = { neighborDataString };
  } groupKeys;


  dataRepository::view_rtype<integer_array> GhostRank()
  { return this->getData<integer_array>(viewKeys.ghostRank); }

  dataRepository::view_rtype_const<integer_array> GhostRank() const
  { return this->getData<integer_array>(viewKeys.ghostRank); }



  ManagedGroup m_sets;

  globalIndex_array  m_localToGlobalMap;
  map<globalIndex,localIndex>  m_globalToLocalMap;
  integer_array m_isExternal;
  integer_array m_ghostRank;

};


template< bool DOPACK >
int ObjectManagerBase::PackPrivate( buffer_unit_type * & buffer,
                                    array<string> const & wrapperNames,
                                    localIndex_array const & packList,
                                    integer const includeGlobalIndices,
                                    integer const recursive ) const
{
  int packedSize = 0;
  packedSize += CommBufferOps::Pack<DOPACK>( buffer, this->getName() );

  int rank=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );
  packedSize += CommBufferOps::Pack<DOPACK>( buffer, rank );


  int const numPackedIndices = packList.size()==0 ? this->size() : packList.size();
  packedSize += CommBufferOps::Pack<DOPACK>( buffer, numPackedIndices );

  packedSize += CommBufferOps::Pack<DOPACK>( buffer, includeGlobalIndices );

  packedSize += CommBufferOps::Pack<DOPACK>( buffer, string("Wrappers") );


  if( includeGlobalIndices )
  {
    globalIndex_array globalIndices;
    globalIndices.resize(numPackedIndices);
    if( packList.size()==0 )
    {
      for( localIndex a=0 ; a<numPackedIndices ; ++a )
      {
        globalIndices[a] = this->m_localToGlobalMap[a];
      }
    }
    else
    {
      for( localIndex a=0 ; a<numPackedIndices ; ++a )
      {
        globalIndices[a] = this->m_localToGlobalMap[packList[a]];
      }
    }
    packedSize += CommBufferOps::Pack<DOPACK>( buffer, globalIndices );
  }

  array<string> wrapperNamesForPacking;
  if( wrapperNames.size()==0 )
  {
    set<localIndex> exclusionList;
    ViewPackingExclusionList(exclusionList);
    wrapperNamesForPacking.resize( this->wrappers().size() );
    localIndex count = 0;
    for( localIndex k=0 ; k<this->wrappers().size() ; ++k )
    {
      if( exclusionList.count(k) == 0)
      {
        string junk  = wrappers().values()[k].first;
        wrapperNamesForPacking[count++] = wrappers().values()[k].first;
      }
    }
    wrapperNamesForPacking.resize(count);
  }
  else
  {
    wrapperNamesForPacking.resize( wrapperNames.size() );
    for( localIndex count=0 ; count<wrapperNames.size() ; ++count )
    {
      wrapperNamesForPacking[count+1] = wrapperNames[count];
    }
  }

  packedSize += CommBufferOps::Pack<DOPACK,int>( buffer, wrapperNamesForPacking.size() );
  for( auto const & wrapperName : wrapperNamesForPacking )
  {
    dataRepository::ViewWrapperBase const * const wrapper = this->getWrapperBase(wrapperName);
    packedSize += CommBufferOps::Pack<DOPACK>( buffer, wrapperName );
    if( packList.empty() )
    {
      if(DOPACK)
      {
        packedSize += wrapper->Pack( buffer );
      }
      else
      {
        packedSize += wrapper->PackSize();
      }
    }
    else
    {
      if(DOPACK)
      {
        packedSize += wrapper->Pack( buffer, packList );
      }
      else
      {
        packedSize += wrapper->PackSize( packList );
      }
    }
  }


  if( recursive > 0 )
  {
    packedSize += CommBufferOps::Pack<DOPACK>( buffer, string("SubGroups") );
    packedSize += CommBufferOps::Pack<DOPACK>( buffer, this->GetSubGroups().size() );
    for( auto const & keyGroupPair : this->GetSubGroups() )
    {
      packedSize += CommBufferOps::Pack<DOPACK>( buffer, keyGroupPair.first );
      packedSize += keyGroupPair.second->Pack( buffer, wrapperNames, packList, recursive );
    }
  }

  return packedSize;
}


//template< bool DOPACK >
//int ObjectManagerBase::UnpackPrivate( buffer_unit_type const *& buffer,
//                                      localIndex_array const & packList,
//                                      integer const recursive )
//{
//  int unpackedSize=0;
//
//  string groupName;
//  unpackedSize += CommBufferOps::Unpack( buffer, groupName );
//  GEOS_ASSERT( groupName==this->getName(), "ManagedGroup::Unpack(): group names do not match")
//
//  string wrappersLabel;
//  unpackedSize += CommBufferOps::Unpack( buffer, wrappersLabel);
//  GEOS_ASSERT( wrappersLabel=="Wrappers", "ManagedGroup::Unpack(): wrapper label incorrect")
//
//
//  globalIndex_array globalIndices;
//
//  return unpackedSize;
//}

} /* namespace geosx */

typedef geosx::ObjectManagerBase ObjectDataStructureBaseT;

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_OBJECTMANAGERBASE_HPP_ */
