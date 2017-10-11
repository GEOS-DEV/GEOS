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


private:
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


  globalIndex_array & m_localToGlobalMap;
  map<globalIndex,localIndex> & m_globalToLocalMap;


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



  void ConstructLocalListOfBoundaryObjects( localIndex_array & objectList ) const;
  void ConstructGlobalListOfBoundaryObjects( globalIndex_array & objectList ) const;

  //**********************************************************************************************************************


};

} /* namespace geosx */

typedef geosx::ObjectManagerBase ObjectDataStructureBaseT;

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_OBJECTMANAGERBASE_HPP_ */
