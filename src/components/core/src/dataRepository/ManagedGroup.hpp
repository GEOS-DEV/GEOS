/**
 * @file DataObjectManager.h
 * @date created on Nov 21, 2014
 * @author Randolph R. Settgast
 */


#ifndef MANAGEDGROUP_H_
#define MANAGEDGROUP_H_

#include <iostream>
#include <slic/slic.hpp>
#include <mpi.h>

#include "ObjectCatalog.hpp"
#include "ViewWrapper.hpp"

#include "depricated/Common.h"
#include "DocumentationNode.hpp"


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
  using subGroupMap = map< string, std::unique_ptr<ManagedGroup> >;

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
  T& RegisterGroup( std::string const & name, std::unique_ptr<TBASE> newObject );

  template< typename T = ManagedGroup, typename TBASE = ManagedGroup >
  T& RegisterGroup( std::string const & name )
  {
//    T* temp = dynamic_cast<T*>(this);
    return RegisterGroup<T>( name, std::move(std::make_unique< T >( name, this )) );
//    return RegisterGroup<T>( name, std::move(std::make_unique< T >( name, this )) );
  }

  template< typename T = ManagedGroup, typename TBASE = ManagedGroup >
  T& RegisterGroup( std::string const & name, std::string const & catalogName )
  {
//    T* temp = dynamic_cast<T*>(this);
    std::unique_ptr<TBASE> newGroup = TBASE::CatalogInterface::Factory(catalogName, name, this );
//    std::unique_ptr<T> newGroup = T::CatalogInterface::Factory(catalogName, name, (this) );
    return RegisterGroup<T,TBASE>( name, std::move(newGroup) );
  }



  template< typename T = ManagedGroup >
  T& GetGroup( std::string const & name )
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T&>( *(m_subGroups.at(name)) );
#else
    return static_cast<T&>( *(m_subGroups.at(getName)) );
#endif
  }



  template< typename T = ManagedGroup >
  T const & GetGroup( std::string const & name ) const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T const &>( *(m_subGroups.at(name)) );
#else
    return static_cast<T const &>( *(m_subGroups.at(getName)) );
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
       T & subGroup = dynamic_cast<T &>( *(subGroupIter.second) );
#else
       T & subGroup = static_cast<T &>( *(subGroupIter.second) );
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
       T const & subGroup = dynamic_cast<T const &>( *(subGroupIter.second) );
#else
       T const & subGroup = static_cast<T const &>( *(subGroupIter.second) );
#endif
       lambda( subGroup );
    }
  }

  virtual void Initialize( ManagedGroup * const group );

  virtual void InitializationOrder( string_array & order );

  virtual void InitializePreSubGroups( ManagedGroup * const group ) {}

  virtual void InitializePostSubGroups( ManagedGroup * const group ) {}


  template< typename T , typename TBASE=T >
  ViewWrapper<TBASE>& RegisterViewWrapper( std::string const & name, std::size_t * const rkey = nullptr );


  ViewWrapperBase& RegisterViewWrapper( std::string const & name, rtTypes::TypeIDs const & type );

  template< typename T >
  ViewWrapper<T>& RegisterViewWrapper( std::string const & name, std::unique_ptr<T> newObject );


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


  ViewWrapperBase const & getWrapperBase( size_t const index ) const
  {
    return *(m_wrappers[index]);
  }

  ViewWrapperBase & getWrapperBase( size_t const index )
  {
    return *(m_wrappers[index]);
  }

  ViewWrapperBase const & getWrapperBase( std::string const & name ) const
  {
    auto index = m_keyLookup.at(name);
    return getWrapperBase(index);
  }

  ViewWrapperBase & getWrapperBase( std::string const & name )
  {
    auto index = m_keyLookup.at(name);
    return getWrapperBase(index);
  }

  template< typename T >
  ViewWrapper<T> const * getWrapperPtr( std::size_t const index ) const
  {
#ifdef USE_DYNAMIC_CASTING
    ViewWrapperBase const * const temp = m_wrappers[index].get();
    ViewWrapper<T> const * temp2 = dynamic_cast<ViewWrapper<T> const * >(temp);
    return dynamic_cast< ViewWrapper<T> const * >( (m_wrappers[index]).get() );
#else
    return static_cast< ViewWrapper<T> const * >( (m_wrappers[index]).get() );
#endif
  }

  template< typename T >
  ViewWrapper<T> * getWrapperPtr( std::size_t const index )
  {
    return const_cast<ViewWrapper<T> *>( const_cast< ManagedGroup const *>(this)->getWrapperPtr<T>( index ) );
  }

  template< typename T >
  ViewWrapper<T> const * getWrapperPtr( std::string const & name ) const
  {
    auto index = m_keyLookup.at(name);
    return getWrapperPtr<T>(index);
  }

  template< typename T >
  ViewWrapper<T> * getWrapperPtr( std::string const & name )
  { return const_cast<ViewWrapper<T> *>( const_cast<const ManagedGroup*>(this)->getWrapperPtr<T>( name ) ); }



  template< typename T >
  ViewWrapper<T> const & getWrapper( std::size_t const index ) const
  {
    return *getWrapperPtr<T>(index);
  }

  template< typename T >
  ViewWrapper<T> & getWrapper( std::size_t const index )
  {
    return const_cast<ViewWrapper<T>&>( const_cast< ManagedGroup const *>(this)->getWrapper<T>( index ) );
  }

  template< typename T >
  ViewWrapper<T> const & getWrapper( std::string const & name ) const
  {
    auto index = m_keyLookup.at(name);
    return getWrapper<T>(index);
  }

  template< typename T >
  ViewWrapper<T>& getWrapper( std::string const & name )
  { return const_cast<ViewWrapper<T>&>( const_cast<const ManagedGroup*>(this)->getWrapper<T>( name ) ); }



  template< typename T >
  typename ViewWrapper<T>::rtype_const getData( size_t const index ) const
  {
    return getWrapper<T>(index).data();
  }

  template< typename T >
  typename ViewWrapper<T>::rtype getData( size_t const index )
  {
    return getWrapper<T>(index).data();
//    return const_cast<typename WrapperView<T>::rtype>( const_cast<const SynchronizedGroup*>(this)->getData<T>( index ) );
  }

  template< typename T >
  typename ViewWrapper<T>::rtype_const getData( std::string const & name ) const
  {
    auto index = m_keyLookup.at(name);
    return getData<T>( index );
  }

  template< typename T >
  typename ViewWrapper<T>::rtype getData( std::string const & name )
  {
    auto index = m_keyLookup.at(name);
    return getData<T>( index );
  }



  template< typename T >
  T const & getReference( std::size_t const index ) const
  {
    return getWrapper<T>(index).reference();
  }

  template< typename T >
  T& getReference( std::size_t const index )
  {
    return const_cast<T&>( const_cast<const ManagedGroup*>(this)->getReference<T>( index ) );
  }

  template< typename T >
  T const & getReference( std::string const & name ) const
  {
    auto index = m_keyLookup.at(name);
    return getReference<T>( index );
  }

  template< typename T >
  T & getReference( std::string const & name )
  {
    auto index = m_keyLookup.at(name);
    return getReference<T>( index );
  }


  ManagedGroup * hasGroup( std::string const & name )
  {
    ManagedGroup * rval = nullptr;
    auto iter = m_subGroups.find(name);
    if( iter!=m_subGroups.end() )
    {
      rval = iter->second.get();
    }
    return rval;
  }

  ManagedGroup const * hasGroup( std::string const & name ) const
  {
    ManagedGroup const * rval = nullptr;
    auto iter = m_subGroups.find(name);
    if( iter!=m_subGroups.end() )
    {
      rval = iter->second.get();
    }
    return rval;
  }

  bool hasView( std::string const & name )
  {
    return m_keyLookup.count(name);
  }
  bool hasView( std::string const & name ) const
  {
    return m_keyLookup.count(name);
  }

  inline string getName() const
  {
    return getData<string>(keys::Name);
  }

  virtual void resize( localIndex newsize );

  inline localIndex size() const
  {
    return *(getData<localIndex>(keys::Size));
  }



  axom::sidre::Group * getSidreGroup()              { return m_sidreGroup; }
  axom::sidre::Group const * getSidreGroup() const  { return m_sidreGroup; }

  static axom::sidre::Group * setSidreGroup( string const& name,
                                                       ManagedGroup * const parent );

  ManagedGroup * getParent()             { return m_parent; }
  ManagedGroup const * getParent() const { return m_parent; }

  ManagedGroup * setParent( ManagedGroup * const parent )
  {
    m_parent = parent;
    m_sidreGroup = m_parent->getSidreGroup();

    return m_parent;
  }

  std::vector< std::unique_ptr<ViewWrapperBase> > const & wrappers() const
  {
    return m_wrappers;
  }

  std::vector< std::unique_ptr<ViewWrapperBase> > & wrappers()
  {
    return m_wrappers;
  }


void registerSubViews();

void unregisterSubViews();

void writeRestart(int num_files, const string & path, const string & protocol, MPI_Comm comm);

void reconstructSidreTree(const string & root_path, const string & protocol, MPI_Comm comm);

void resizeSubViews();

void loadSidreExternalData(const string & root_path, MPI_Comm comm);



protected:
  cxx_utilities::DocumentationNode * m_docNode = nullptr;

private:
  unordered_map<string,size_t> m_keyLookup;
  std::vector< std::unique_ptr<ViewWrapperBase> > m_wrappers;

  ManagedGroup* m_parent = nullptr;
  subGroupMap m_subGroups;

  axom::sidre::Group* m_sidreGroup;

  int32 const & m_size;
  string const & m_name;
//  string const & m_path;


  /**
   * @name functions to disallow construction of strings from char const *
   */
  ///@{
#if NOCHARTOSTRING_KEYLOOKUP == 1

  template< typename T = ManagedGroup >
  T const & GetGroup( char const * ) const;

  template< typename T = ManagedGroup >
  T& GetGroup( char const * name );


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



template < typename T, typename TBASE >
T& ManagedGroup::RegisterGroup( std::string const & name, std::unique_ptr<TBASE> newObject )
{
  auto iterKeyLookup = m_subGroups.find(name);

  // if the key was not found, make DataObject<T> and insert
  if( iterKeyLookup == m_subGroups.end() )
  {
    auto insertResult = m_subGroups.insert( std::make_pair( name, std::move(newObject) ) );

    if( !insertResult.second )
    {
      std::cout<<LOCATION<<std::endl;
      throw std::exception();
    }
    iterKeyLookup = insertResult.first;
//    iterKeyLookup->second.get()->
  }
  // if key was found, make sure that they are the same type
  else
  {
    if( typeid(T) != iterKeyLookup->second->get_typeid() )
    {
      std::cout<<LOCATION<<std::endl;
      throw std::exception();
    }
  }
#ifdef USE_DYNAMIC_CASTING
  return *(dynamic_cast<T*>( (iterKeyLookup->second).get() ) );
#else
  return *(static_cast<T*>( (iterKeyLookup->second).get() ) );
#endif
}



template< typename T , typename TBASE >
ViewWrapper<TBASE>& ManagedGroup::RegisterViewWrapper( std::string const & name, std::size_t * const rkey )
{
  std::size_t key = static_cast<std::size_t>(-1);

  auto iterKeyLookup = m_keyLookup.find(name);

  // if the key was not found, make DataObject<T> and insert
  if( iterKeyLookup == m_keyLookup.end() )
  {
    m_wrappers.push_back( std::move( ViewWrapper<TBASE>::template Factory<T>(name,this) ) );
    key = m_wrappers.size() - 1;
    m_keyLookup.insert( std::make_pair(name,key) );
    if( m_wrappers[key]->sizedFromParent() == 1 )
    {
      m_wrappers[key]->resize(this->size());
    }
  }
  // if key was found, make sure that they are the same type
  else
  {
    key = m_keyLookup.at(name);
    auto& basePtr = m_wrappers[key];
    if( typeid(T) != basePtr->get_typeid() )
    {
      std::string error = string("Call to Group::RegisterViewWrapper( ")
                          +name+string(", std::size_t * const ) attempts to re-register ViewWrapper<")
                          +basePtr->get_typeid().name()+string(", but with different type") ;
      SLIC_ERROR(error);
//      throw std::exception();
    }
  }

  if( rkey != nullptr )
  {
    *rkey = key;
  }
  return getWrapper<TBASE>(key);
}


template < typename T >
ViewWrapper<T>& ManagedGroup::RegisterViewWrapper( std::string const & name, std::unique_ptr<T> newObject )
{
  std::size_t key = static_cast<std::size_t>(-1);
  auto iterKeyLookup = m_keyLookup.find(name);
//  auto iterKeyLookup = m_wrappers.find(name);

  // if the key was not found, make DataObject<T> and insert
  if( iterKeyLookup == m_keyLookup.end() )
  {
    std::unique_ptr< ViewWrapper<T> > newWrapper = std::make_unique< ViewWrapper<T> >( name, this, std::move(newObject) );
    m_wrappers.push_back( std::move( newWrapper ) );
    key = m_wrappers.size() - 1;
    m_keyLookup.insert( std::make_pair(name,key) );
    m_wrappers[key]->resize(this->size());
  }
  // if key was found, make sure that they are the same type
  else
  {
    key = m_keyLookup.at(name);
    auto& basePtr = m_wrappers[key];
    if( typeid(T) != basePtr->get_typeid() )
    {
      std::string error = string("Call to Group::RegisterViewWrapper( ")
                          +name+string(", std::size_t * const ) attempts to re-register ViewWrapper<")
                          +basePtr->get_typeid().name()+string(", but with different type") ;
      SLIC_ERROR(error);
    }
  }
  return getWrapper<T>(key);
}

} // namespace dataRepository
} /* namespace geosx */


//typedef geosx::dataRepository::ManagedGroup ObjectDataStructureBaseT;

#endif /* MANAGEDGROUP_H_ */
