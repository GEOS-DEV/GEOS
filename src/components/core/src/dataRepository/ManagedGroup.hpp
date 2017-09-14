/**
 * @file DataObjectManager.h
 * @date created on Nov 21, 2014
 * @author Randolph R. Settgast
 */


#ifndef MANAGEDGROUP_H_
#define MANAGEDGROUP_H_

#include <iostream>
#include <mpi.h>

#include "ObjectCatalog.hpp"
#include "ViewWrapper.hpp"

#include "depricated/Common.h"
#include "DocumentationNode.hpp"

#include "MappedVector.hpp"

//#include "CodingUtilities/ANSTexception.hpp"

#ifndef USE_DYNAMIC_CASTING
#define USE_DYNAMIC_CASTING 1;
#endif

#ifndef NOCHARTOSTRING_KEYLOOKUP
#define NOCHARTOSTRING_KEYLOOKUP 1
#endif

/**
 * namespace to encapsulate functions in simulation tools
 */
namespace geosx
{
namespace dataRepository
{

using keyType = string;
using indexType = int;
//using DataKey = DataKeyT<keyType,indexType>;

/**
 * @author Randolph R. Settgast
 *
 * class that encapsulates and manages a collection of DataObjects. Can be considered a "node" in a
 * hierarchy of managers that represent physical groupings of data.
 *
 */
class ManagedGroup
{
public:
//  using subGroupMap = map< string, std::unique_ptr<ManagedGroup> >;
  using subGroupMap = MappedVector< ManagedGroup, std::unique_ptr<ManagedGroup>, keyType, indexType  >;
  using viewWrapperMap = MappedVector< ViewWrapperBase, std::unique_ptr<ViewWrapperBase>, keyType, indexType  >;
  /**
   * @name constructors, destructor, copy, move, assignments
   */
  ///@{

  /**
   * @author Randolph R. Settgast
   * @param name the name of this object manager
   */
  explicit ManagedGroup( std::string const & name,
                         ManagedGroup * const parent );

//  explicit ManagedGroup( std::string const & name,
//                         ManagedGroup * const parent,
//                         cxx_utilities::DocumentationNode * docNode );

  /**
   *
   */
  virtual ~ManagedGroup();

  /**
   *
   * @param source source WrapperCollection
   */
  ManagedGroup( ManagedGroup&& source );


  ManagedGroup() = delete;
  ManagedGroup( ManagedGroup const & source ) = delete;
  ManagedGroup& operator=( ManagedGroup const & ) = delete;
  ManagedGroup& operator=(ManagedGroup&&) = delete;

  ///@}


  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  using CatalogInterface = cxx_utilities::CatalogInterface< ManagedGroup, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();
  ///@}



  virtual const std::type_info& get_typeid() const
  {
    return typeid(*this);
  }


  template< typename T = ManagedGroup, typename TBASE = ManagedGroup >
  T * RegisterGroup( std::string const & name, std::unique_ptr<TBASE> newObject );

  template< typename T = ManagedGroup, typename TBASE = ManagedGroup >
  T * RegisterGroup( std::string const & name )
  {
    return RegisterGroup<T>( name, std::move(std::make_unique< T >( name, this )) );
  }

  template< typename T = ManagedGroup, typename TBASE = ManagedGroup >
  T * RegisterGroup( subGroupMap::KeyIndex & keyIndex )
  {
    T * rval = RegisterGroup<T>( keyIndex.Key(), std::move(std::make_unique< T >( keyIndex.Key(), this )) );
    keyIndex.setIndex( this->m_subGroups.getIndex(keyIndex.Key()) );
    return rval;
  }

  template< typename T = ManagedGroup, typename TBASE = ManagedGroup >
  T * RegisterGroup( std::string const & name, std::string const & catalogName )
  {
    std::unique_ptr<TBASE> newGroup = TBASE::CatalogInterface::Factory(catalogName, name, this );
    return RegisterGroup<T,TBASE>( name, std::move(newGroup) );
  }



  template< typename T = ManagedGroup >
  T * GetGroup( std::string const & name )
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T *>( m_subGroups[name] );
#else
    return static_cast<T *>( m_subGroups[name] );
#endif
  }

  template< typename T = ManagedGroup >
  T const * GetGroup( std::string const & name ) const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T const *>( m_subGroups[name] );
#else
    return static_cast<T const *>( m_subGroups[name] );
#endif
  }


  template< typename T = ManagedGroup >
  T * GetGroup( subGroupMap::KeyIndex & key )
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T *>( m_subGroups[key] );
#else
    return static_cast<T *>( m_subGroups[key] );
#endif
  }

  template< typename T = ManagedGroup >
  T const * GetGroup( subGroupMap::KeyIndex & key ) const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T const *>( m_subGroups[key] );
#else
    return static_cast<T const *>( m_subGroups[key] );
#endif
  }


  subGroupMap & GetSubGroups()
  {
    return m_subGroups;
  }

  subGroupMap const & GetSubGroups() const
  {
    return m_subGroups;
  }

  template< typename T = ManagedGroup, typename LAMBDA >
  void forSubGroups( LAMBDA lambda )
  {
    for( auto& subGroupIter : m_subGroups )
    {
#ifdef USE_DYNAMIC_CASTING
       T * subGroup = dynamic_cast<T *>( subGroupIter.second.get() );
#else
       T * subGroup = static_cast<T *>( subGroupIter.second.get() );
#endif
       lambda( subGroup );
    }
  }

  template< typename T = ManagedGroup, typename LAMBDA >
  void forSubGroups( LAMBDA lambda ) const
  {
    for( auto const & subGroupIter : m_subGroups )
    {
#ifdef USE_DYNAMIC_CASTING
       T const * subGroup = dynamic_cast<T const *>( subGroupIter.second );
#else
       T const * subGroup = static_cast<T const *>( subGroupIter.second );
#endif
       lambda( subGroup );
    }
  }

  virtual void Initialize( ManagedGroup * const group );

  virtual void InitializationOrder( string_array & order );

  virtual void InitializePreSubGroups( ManagedGroup * const group ) {}

  virtual void InitializePostSubGroups( ManagedGroup * const group ) {}


  template< typename T , typename TBASE=T >
  ViewWrapper<TBASE> * RegisterViewWrapper( std::string const & name, viewWrapperMap::KeyIndex::index_type * const rkey = nullptr );

  template< typename T , typename TBASE=T >
  ViewWrapper<TBASE> * RegisterViewWrapper( ManagedGroup::viewWrapperMap::KeyIndex & viewKey );


  ViewWrapperBase * RegisterViewWrapper( std::string const & name, rtTypes::TypeIDs const & type );

  template< typename T >
  ViewWrapper<T> * RegisterViewWrapper( std::string const & name, std::unique_ptr<T> newObject );


  ///@}


  /**
   * @name Self Documentation Functions
   */
  ///@{

  cxx_utilities::DocumentationNode * getDocumentationNode()
  {
    return m_docNode;
  }

  void RegisterDocumentationNodes();


  ///@}
  ///


  void PrintDataHierarchy();

  virtual void ReadXML( xmlWrapper::xmlNode const & targetNode );

  virtual void ReadXMLsub( xmlWrapper::xmlNode const & );

  virtual void ReadXML_PostProcess() {}

  virtual void BuildDataStructure( dataRepository::ManagedGroup * const rootGroup );

  virtual void FillDocumentationNode( dataRepository::ManagedGroup * const group );

  void SetDocumentationNodes( dataRepository::ManagedGroup * const group );



  //***********************************************************************************************

  // user defined conversion doesn't work. can't infer template argument
//  class GetDataClass
//  {
//  public:
//    GetDataClass( ManagedGroup & parent ): m_parent( parent ) {}
//
//    inline GetDataClass& operator() ( std::string const & name )
//    {
//      m_name = name;
//      return *this;
//    }
//
//    template< typename T>
//    operator typename ViewWrapper<T>::rtype ()
//    {
//      return m_parent.getData<T>( m_name );
//    }
//  private:
//    ManagedGroup & m_parent;
//    std::string m_name;
//  };
//  GetDataClass GetData = {*this};


  ViewWrapperBase const * getWrapperBase( size_t const index ) const
  { return m_wrappers[index]; }

  ViewWrapperBase * getWrapperBase( size_t const index )
  { return m_wrappers[index]; }

  ViewWrapperBase const * getWrapperBase( std::string const & name ) const
  { return m_wrappers[name]; }

  ViewWrapperBase * getWrapperBase( std::string const & name )
  { return m_wrappers[name]; }

  ViewWrapperBase const * getWrapperBase( viewWrapperMap::KeyIndex & keyIndex ) const
  { return m_wrappers[keyIndex]; }

  ViewWrapperBase * getWrapperBase( viewWrapperMap::KeyIndex & keyIndex )
  { return m_wrappers[keyIndex]; }


  template< typename T >
  ViewWrapper<T> const * getWrapper( std::size_t const index ) const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast< ViewWrapper<T> const * >( (m_wrappers[index]) );
#else
    return static_cast< ViewWrapper<T> const * >( (m_wrappers[index]) );
#endif
  }

  template< typename T >
  ViewWrapper<T> * getWrapper( std::size_t const index )
  { return const_cast<ViewWrapper<T> *>( const_cast< ManagedGroup const *>(this)->getWrapper<T>( index ) ); }

  template< typename T >
  ViewWrapper<T> const * getWrapper( std::string const & name ) const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast< ViewWrapper<T> const * >( (m_wrappers[name]) );
#else
    return static_cast< ViewWrapper<T> const * >( (m_wrappers[name]) );
#endif
  }

  template< typename T >
  ViewWrapper<T> * getWrapper( std::string const & name )
  { return const_cast<ViewWrapper<T> *>( const_cast<const ManagedGroup*>(this)->getWrapper<T>( name ) ); }


  template< typename T >
  ViewWrapper<T> const * getWrapper( viewWrapperMap::KeyIndex & keyIndex ) const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast< ViewWrapper<T> const * >( (m_wrappers[keyIndex]) );
#else
    return static_cast< ViewWrapper<T> const * >( (m_wrappers[keyIndex]) );
#endif
  }

  template< typename T >
  ViewWrapper<T> * getWrapper( viewWrapperMap::KeyIndex & keyIndex )
  { return const_cast<ViewWrapper<T> *>( const_cast<const ManagedGroup*>(this)->getWrapper<T>( keyIndex ) ); }


  template< typename T >
  ViewWrapper<T> const * getWrapper( viewWrapperMap::KeyIndex const & keyIndex ) const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast< ViewWrapper<T> const * >( (m_wrappers[keyIndex]) );
#else
    return static_cast< ViewWrapper<T> const * >( (m_wrappers[keyIndex]) );
#endif
  }

  template< typename T >
  ViewWrapper<T> * getWrapper( viewWrapperMap::KeyIndex const & keyIndex )
  { return const_cast<ViewWrapper<T> *>( const_cast<const ManagedGroup*>(this)->getWrapper<T>( keyIndex ) ); }






  template< typename T >
  view_rtype_const<T> getData( size_t const index ) const
  { return getWrapper<T>(index)->data(); }

  template< typename T >
  view_rtype<T> getData( size_t const index )
  { return getWrapper<T>(index)->data(); }

  template< typename T >
  view_rtype_const<T> getData( std::string const & name ) const
  { return getWrapper<T>( name )->data(); }

  template< typename T >
  view_rtype<T> getData( std::string const & name )
  { return getWrapper<T>( name )->data(); }

  template< typename T >
  view_rtype_const<T> getData( viewWrapperMap::KeyIndex & keyIndex ) const
  { return getWrapper<T>( keyIndex )->data(); }

  template< typename T >
  view_rtype<T> getData( viewWrapperMap::KeyIndex & keyIndex )
  { return getWrapper<T>( keyIndex )->data(); }

  template< typename T >
  view_rtype_const<T> getData( viewWrapperMap::KeyIndex const & keyIndex ) const
  { return getWrapper<T>( keyIndex )->data(); }

  template< typename T >
  view_rtype<T> getData( viewWrapperMap::KeyIndex const & keyIndex )
  { return getWrapper<T>( keyIndex )->data(); }

//  /**
//   *
//   * @param keyIndex
//   * @return
//   * @note BREAKS const correctness for keyIndex
//   */
//  template< typename T >
//  view_rtype_const<T> getData( viewWrapperMap::KeyIndex const & keyIndex ) const
//  { return getWrapper<T>( const_cast<viewWrapperMap::KeyIndex & >(keyIndex) )->data(); }
//
//  /**
//   *
//   * @param keyIndex
//   * @return
//   * @note BREAKS const correctness for keyIndex
//   */
//  template< typename T >
//  view_rtype<T> getData( viewWrapperMap::KeyIndex const & keyIndex )
//  { return getWrapper<T>( const_cast<viewWrapperMap::KeyIndex & >(keyIndex) )->data(); }


  template< typename T >
  T const & getReference( std::size_t const index ) const
  { return getWrapper<T>(index)->reference(); }

  template< typename T >
  T& getReference( std::size_t const index )
  { return const_cast<T&>( const_cast<const ManagedGroup*>(this)->getReference<T>( index ) ); }

  template< typename T >
  T const & getReference( std::string const & name ) const
  { return getWrapper<T>(name)->reference(); }

  template< typename T >
  T & getReference( std::string const & name )
  { return getWrapper<T>(name)->reference(); }

  template< typename T >
  T const & getReference( viewWrapperMap::KeyIndex & keyIndex ) const
  { return getWrapper<T>(keyIndex)->reference(); }

  template< typename T >
  T & getReference( viewWrapperMap::KeyIndex & keyIndex )
  { return getWrapper<T>(keyIndex)->reference(); }


  bool hasGroup( std::string const & name ) const
  {
    return (m_subGroups[name] != nullptr);
  }

  bool hasView( std::string const & name ) const
  {
    return (m_wrappers[name] != nullptr);
  }

  inline const string getName() const
  {
    return m_name;
  }

  virtual void resize( localIndex newsize );

  inline localIndex size() const
  {
    return m_size;
  }


#if ATK_FOUND
  axom::sidre::Group * getSidreGroup()              { return m_sidreGroup; }
  axom::sidre::Group const * getSidreGroup() const  { return m_sidreGroup; }

  static axom::sidre::Group * setSidreGroup( string const& name,
                                                       ManagedGroup * const parent );
#endif

  ManagedGroup * getParent()             { return m_parent; }
  ManagedGroup const * getParent() const { return m_parent; }

  ManagedGroup * setParent( ManagedGroup * const parent )
  {
    m_parent = parent;
#if ATK_FOUND
    m_sidreGroup = m_parent->getSidreGroup();
#endif

    return m_parent;
  }

  viewWrapperMap const & wrappers() const
  {
    return m_wrappers;
  }

  viewWrapperMap & wrappers()
  {
    return m_wrappers;
  }


  void writeRestart(int num_files, const string & path, const string & protocol, MPI_Comm comm);

#if ATK_FOUND
  void reconstructSidreTree(const string & root_path, const string & protocol, MPI_Comm comm);

  void loadSidreExternalData(const string & root_path, MPI_Comm comm);
#endif

protected:
  cxx_utilities::DocumentationNode * m_docNode = nullptr;

private:

#if ATK_FOUND
  void registerSubViews();

  void createSizeViews();

  void loadSizeViews();

  void unregisterSubViews();

  void resizeSubViews();

  void storeSizedFromParent();

  void loadSizedFromParent();
#endif
  

  ManagedGroup* m_parent = nullptr;
  viewWrapperMap m_wrappers;
  subGroupMap m_subGroups;

#if ATK_FOUND
  axom::sidre::Group* m_sidreGroup;
#endif

  int32 m_size;

  string m_name;


  /**
   * @name functions to disallow construction of strings from char const *
   */
  ///@{
#if NOCHARTOSTRING_KEYLOOKUP == 1

  template< typename T = ManagedGroup >
  T const * GetGroup( char const * ) const;

  template< typename T = ManagedGroup >
  T * GetGroup( char const * name );


  template< typename T >
  typename ViewWrapper<T>::rtype_const getData( char const * ) const;

  template< typename T >
  typename ViewWrapper<T>::rtype getData( char const * );

  template< typename T >
  T const & getReference( char const * ) const;

  template< typename T >
  T & getReference( char const * );


  template< typename T >
  ViewWrapper<T> const & getWrapper( char const * ) const;

  template< typename T >
  ViewWrapper<T>& getWrapper( char const * );

  ///@}


#endif
};

using GroupKey = ManagedGroup::subGroupMap::KeyIndex;
using ViewKey = ManagedGroup::viewWrapperMap::KeyIndex;



template < typename T, typename TBASE >
T * ManagedGroup::RegisterGroup( std::string const & name, std::unique_ptr<TBASE> newObject )
{
  #ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T*>( m_subGroups.insert( name, std::move(newObject) ) );
  #else
    return static_cast<T*>( m_subGroups.insert( name, std::move(newObject) ) );
  #endif
}



template< typename T , typename TBASE >
ViewWrapper<TBASE> * ManagedGroup::RegisterViewWrapper( std::string const & name, ViewKey::index_type * const rkey )
{
  m_wrappers.insert( name, std::move(ViewWrapper<TBASE>::template Factory<T>(name,this) ) );
  ViewWrapper<TBASE> * const rval = getWrapper<TBASE>(name);
  if( rval->sizedFromParent() == 1 )
  {
    rval->resize(this->size());
  }
  return rval;
}

template< typename T , typename TBASE >
ViewWrapper<TBASE> * ManagedGroup::RegisterViewWrapper( ViewKey & viewKey )
{
  ViewKey::index_type index;
  ViewWrapper<TBASE> * const rval = RegisterViewWrapper<T,TBASE>( viewKey.Key(), &index );
  viewKey.setIndex(index);

  return rval;
}


template < typename T >
ViewWrapper<T> * ManagedGroup::RegisterViewWrapper( std::string const & name, std::unique_ptr<T> newObject )
{
  m_wrappers.insert( name, std::make_unique< ViewWrapper<T> >( name, this, std::move(newObject) ) );

  ViewWrapper<T> * const rval = getWrapper<T>(name);
  if( rval->sizedFromParent() == 1 )
  {
    rval->resize(this->size());
  }
  return rval;
}


} /* end namespace dataRepository */
} /* end namespace geosx */


//typedef geosx::dataRepository::ManagedGroup ObjectDataStructureBaseT;

#endif /* MANAGEDGROUP_H_ */
