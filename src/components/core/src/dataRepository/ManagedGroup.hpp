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

#include "MapVectorContainer.hpp"

#include "DataKey.hpp"

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
using DataKey = DataKeyT<keyType,indexType>;

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
  using subGroupMap = MapVectorContainer< ManagedGroup, std::unique_ptr<ManagedGroup>, keyType, indexType  >;
  using viewWrapperMap = MapVectorContainer< ViewWrapperBase, std::unique_ptr<ViewWrapperBase>, keyType, indexType  >;
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
  T * GetGroupPtr( std::string const & name )
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T *>( m_subGroups[name] );
#else
    return static_cast<T *>( m_subGroups[name] );
#endif
  }

  template< typename T = ManagedGroup >
  T const * GetGroupPtr( std::string const & name ) const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T const *>( m_subGroups[name] );
#else
    return static_cast<T const *>( m_subGroups[name] );
#endif
  }




  template< typename T = ManagedGroup >
  T& GetGroup( std::string const & name )
  { return *(GetGroupPtr<T>(name)); }

  template< typename T = ManagedGroup >
  T const & GetGroup( std::string const & name ) const
  { return *(GetGroupPtr<T>(name)); }





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
    for( auto& subGroupIter : m_subGroups.objects() )
    {
#ifdef USE_DYNAMIC_CASTING
       T & subGroup = dynamic_cast<T &>( *(subGroupIter) );
#else
       T & subGroup = static_cast<T &>( *(subGroupIter) );
#endif
       lambda( subGroup );
    }
  }

  template< typename T = ManagedGroup, typename LAMBDA >
  void forSubGroups( LAMBDA lambda ) const
  {
    for( auto const & subGroupIter : m_subGroups.objects() )
    {
#ifdef USE_DYNAMIC_CASTING
       T const & subGroup = dynamic_cast<T const &>( *(subGroupIter) );
#else
       T const & subGroup = static_cast<T const &>( *(subGroupIter) );
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
  { return *(m_wrappers[index]); }

  ViewWrapperBase & getWrapperBase( size_t const index )
  { return *(m_wrappers[index]); }

  ViewWrapperBase const & getWrapperBase( std::string const & name ) const
  { return *(m_wrappers[name]); }

  ViewWrapperBase & getWrapperBase( std::string const & name )
  { return *(m_wrappers[name]); }

  ViewWrapperBase const & getWrapperBase( DataKey & dataKey ) const
  { return *(m_wrappers[dataKey]); }

  ViewWrapperBase & getWrapperBase( DataKey & dataKey )
  { return *(m_wrappers[dataKey]); }


  template< typename T >
  ViewWrapper<T> const * getWrapperPtr( std::size_t const index ) const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast< ViewWrapper<T> const * >( (m_wrappers[index]) );
#else
    return static_cast< ViewWrapper<T> const * >( (m_wrappers[index]) );
#endif
  }

  template< typename T >
  ViewWrapper<T> * getWrapperPtr( std::size_t const index )
  { return const_cast<ViewWrapper<T> *>( const_cast< ManagedGroup const *>(this)->getWrapperPtr<T>( index ) ); }

  template< typename T >
  ViewWrapper<T> const * getWrapperPtr( std::string const & name ) const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast< ViewWrapper<T> const * >( (m_wrappers[name]) );
#else
    return static_cast< ViewWrapper<T> const * >( (m_wrappers[name]) );
#endif
  }

  template< typename T >
  ViewWrapper<T> * getWrapperPtr( std::string const & name )
  { return const_cast<ViewWrapper<T> *>( const_cast<const ManagedGroup*>(this)->getWrapperPtr<T>( name ) ); }

  template< typename T >
  ViewWrapper<T> const * getWrapperPtr( DataKey & dataKey ) const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast< ViewWrapper<T> const * >( (m_wrappers[dataKey]) );
#else
    return static_cast< ViewWrapper<T> const * >( (m_wrappers[dataKey]) );
#endif
  }

  template< typename T >
  ViewWrapper<T> * getWrapperPtr( DataKey & dataKey )
  { return const_cast<ViewWrapper<T> *>( const_cast<const ManagedGroup*>(this)->getWrapperPtr<T>( dataKey ) ); }





  template< typename T >
  ViewWrapper<T> const & getWrapper( std::size_t const index ) const
  { return *getWrapperPtr<T>(index); }

  template< typename T >
  ViewWrapper<T> & getWrapper( std::size_t const index )
  { return *getWrapperPtr<T>(index); }

  template< typename T >
  ViewWrapper<T> const & getWrapper( std::string const & name ) const
  { return *getWrapperPtr<T>(name);  }

  template< typename T >
  ViewWrapper<T>& getWrapper( std::string const & name )
  { return *getWrapperPtr<T>(name);  }

  template< typename T >
  ViewWrapper<T> const & getWrapper( DataKey & dataKey ) const
  { return *getWrapperPtr<T>(dataKey);  }

  template< typename T >
  ViewWrapper<T>& getWrapper( DataKey & dataKey )
  { return *getWrapperPtr<T>(dataKey);  }



  template< typename T >
  typename ViewWrapper<T>::rtype_const getData( size_t const index ) const
  { return getWrapper<T>(index).data(); }

  template< typename T >
  typename ViewWrapper<T>::rtype getData( size_t const index )
  { return getWrapper<T>(index).data(); }

  template< typename T >
  typename ViewWrapper<T>::rtype_const getData( std::string const & name ) const
  { return getWrapper<T>( name ).data(); }

  template< typename T >
  typename ViewWrapper<T>::rtype getData( std::string const & name )
  { return getWrapper<T>( name ).data(); }

  template< typename T >
  typename ViewWrapper<T>::rtype_const getData( DataKey & dataKey ) const
  { return getWrapper<T>( dataKey ).data(); }

  template< typename T >
  typename ViewWrapper<T>::rtype getData( DataKey & dataKey )
  { return getWrapper<T>( dataKey ).data(); }



  template< typename T >
  T const & getReference( std::size_t const index ) const
  { return getWrapper<T>(index).reference(); }

  template< typename T >
  T& getReference( std::size_t const index )
  { return const_cast<T&>( const_cast<const ManagedGroup*>(this)->getReference<T>( index ) ); }

  template< typename T >
  T const & getReference( std::string const & name ) const
  { return getWrapper<T>(name).reference(); }

  template< typename T >
  T & getReference( std::string const & name )
  { return getWrapper<T>(name).reference(); }

  template< typename T >
  T const & getReference( DataKey & dataKey ) const
  { return getWrapper<T>(dataKey).reference(); }

  template< typename T >
  T & getReference( DataKey & dataKey )
  { return getWrapper<T>(dataKey).reference(); }


  ManagedGroup const * hasGroup( std::string const & name ) const
  {
    return (m_subGroups[name]);
  }

  bool hasView( std::string const & name ) const
  {
    return m_wrappers[name];
  }

  inline string getName() const
  {
    return m_sidreGroup->getName();
  }

  virtual void resize( localIndex newsize );

  inline localIndex size() const
  {
    return m_size;
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
    return m_wrappers.objects();
  }

  std::vector< std::unique_ptr<ViewWrapperBase> > & wrappers()
  {
    return m_wrappers.objects();
  }


  void writeRestart(int num_files, const string & path, const string & protocol, MPI_Comm comm);

  void reconstructSidreTree(const string & root_path, const string & protocol, MPI_Comm comm);

  void loadSidreExternalData(const string & root_path, MPI_Comm comm);


protected:
  cxx_utilities::DocumentationNode * m_docNode = nullptr;

private:


  void registerSubViews();

  void createSizeViews();

  void loadSizeViews();

  void unregisterSubViews();

  void resizeSubViews();

  void storeSizedFromParent();

  void loadSizedFromParent();
  

  ManagedGroup* m_parent = nullptr;
  viewWrapperMap m_wrappers;
  subGroupMap m_subGroups;

  axom::sidre::Group* m_sidreGroup;

  int32 m_size;


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
  #ifdef USE_DYNAMIC_CASTING
    return *(dynamic_cast<T*>( m_subGroups.insert( name, std::move(newObject) ) ) );
  #else
    return *(static_cast<T*>( m_subGroups.insert( name, std::move(newObject) ) ) );
  #endif
}



template< typename T , typename TBASE >
ViewWrapper<TBASE>& ManagedGroup::RegisterViewWrapper( std::string const & name, std::size_t * const rkey )
{
  m_wrappers.insert( name, std::move(ViewWrapper<TBASE>::template Factory<T>(name,this) ) );
  ViewWrapper<T> & rval = getWrapper<TBASE>(name);
  if( rval.sizedFromParent() == 1 )
  {
    rval.resize(this->size());
  }
  return rval;
}


template < typename T >
ViewWrapper<T>& ManagedGroup::RegisterViewWrapper( std::string const & name, std::unique_ptr<T> newObject )
{
  m_wrappers.insert( name, std::make_unique< ViewWrapper<T> >( name, this, std::move(newObject) ) );

  ViewWrapper<T> & rval = getWrapper<T>(name);
  if( rval.sizedFromParent() == 1 )
  {
    rval.resize(this->size());
  }
  return rval;
}


} // namespace dataRepository
} /* namespace geosx */


//typedef geosx::dataRepository::ManagedGroup ObjectDataStructureBaseT;

#endif /* MANAGEDGROUP_H_ */
