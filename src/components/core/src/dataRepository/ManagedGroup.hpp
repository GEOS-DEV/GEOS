/**
 * @file DataObjectManager.h
 * @date created on Nov 21, 2014
 * @author Randolph R. Settgast
 */


#ifndef MANAGEDGROUP_H_
#define MANAGEDGROUP_H_

#include <iostream>
#include <slic/slic.hpp>

#include "ObjectCatalog.hpp"
#include "ViewWrapper.hpp"

#include "depricated/Common.h"


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


  using CatalogInterface = cxx_utilities::CatalogInterface< ManagedGroup, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();


  virtual void Registration( dataRepository::ManagedGroup * const )
  {}

  virtual const std::type_info& get_typeid() const
  {
    return typeid(*this);
  }


  template< typename T = ManagedGroup >
  T& RegisterGroup( std::string const & name, std::unique_ptr<T> newObject );

  template< typename T = ManagedGroup >
  T& RegisterGroup( std::string const & name )
  {
    return RegisterGroup<T>( name, std::move(std::make_unique< T >( name, this )) );
  }

  template< typename T = ManagedGroup >
  T& RegisterGroup( std::string const & name, std::string const & catalogName )
  {

    std::unique_ptr<T> newGroup = T::CatalogInterface::Factory(catalogName, name, this );
    return RegisterGroup<T>( name, std::move(newGroup) );
  }



  template< typename T = ManagedGroup >
  T& GetGroup( std::string const & name )
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T&>( *(m_subObjectManagers.at(name)) );
#else
    return static_cast<T&>( *(m_subObjectManagers.at(name)) );
#endif
  }



  template< typename T = ManagedGroup >
  T const & GetGroup( std::string const & name ) const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T const &>( *(m_subObjectManagers.at(name)) );
#else
    return static_cast<T const &>( *(m_subObjectManagers.at(name)) );
#endif
  }


  template< typename T >
  ViewWrapper<T>& RegisterViewWrapper( std::string const & name, std::size_t * const rkey = nullptr );


  ViewWrapperBase& RegisterViewWrapper( std::string const & name, rtTypes::TypeIDs const & type );


  //***********************************************************************************************

  // user defined conversion doesn't work. can't infer template argument
  class GetDataClass
  {
  public:
    GetDataClass( ManagedGroup & parent ): m_parent( parent ) {}

    inline GetDataClass& operator() ( std::string const & name )
    {
      m_name = name;
      return *this;
    }

    template< typename T>
    operator typename ViewWrapper<T>::rtype ()
    {
      return m_parent.getData<T>( m_name );
    }
  private:
    ManagedGroup & m_parent;
    std::string m_name;
  };
  GetDataClass GetData = {*this};



  template< typename T >
  ViewWrapper<T> const & getWrapper( std::size_t const index ) const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast< ViewWrapper<T> const & >( *(m_wrappers[index]) );
#else
    return static_cast< ViewWrapper<T> const & >( *(m_wrappers[index]) );
#endif
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
  typename ViewWrapper<T>::rtype_const getData( int32 const index ) const
  {
    return getWrapper<T>(index).data();
  }

  template< typename T >
  typename ViewWrapper<T>::rtype getData( int32 const index )
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


  inline string name() const
  {
    return getData<string>(keys::Name);
  }

  void resize( std::size_t newsize );

  inline localIndex size() const
  {
    return *(getData<localIndex>(keys::Size));
  }



  asctoolkit::sidre::DataGroup * getSidreGroup()              { return m_sidreGroup; }
  asctoolkit::sidre::DataGroup const * getSidreGroup() const  { return m_sidreGroup; }

  asctoolkit::sidre::DataGroup * setSidreGroup()
  {
    return m_sidreGroup;
  }

  ManagedGroup * getParent()             { return m_parent; }
  ManagedGroup const * getParent() const { return m_parent; }

  ManagedGroup * setParent( ManagedGroup * const parent )
  {
    m_parent = parent;
    m_sidreGroup = m_parent->getSidreGroup();

    return m_parent;
  }

private:
  std::unordered_map<std::string,std::size_t> m_keyLookup;
  std::vector< std::unique_ptr<ViewWrapperBase> > m_wrappers;

  ManagedGroup* m_parent = nullptr;
  std::unordered_map< std::string, std::unique_ptr<ManagedGroup> > m_subObjectManagers;

  asctoolkit::sidre::DataGroup* m_sidreGroup;


//****************************************************
// functions for compatibility with old data structure
// TODO Deprecate or modernize all these suckers

public:

  using ObjectType = string;
  class SiloFile;
  localIndex resize( localIndex const newSize,
                     const bool assignGlobals );

  localIndex m_DataLengths;


  localIndex DataLengths() const { return size(); }

  void WriteSilo( SiloFile& siloFile,
                  const std::string& meshname,
                  const int centering,
                  const int cycleNum,
                  const realT problemTime,
                  const bool isRestart,
                  const std::string& multiRoot,
                  const std::string& regionName = "none",
                  const lArray1d& mask = lArray1d() ) const;


  void ReadSilo( const SiloFile& siloFile,
                 const std::string& meshname,
                 const int centering,
                 const int cycleNum,
                 const realT problemTime,
                 const bool isRestart,
                 const std::string& regionName = "none",
                 const lArray1d& mask = lArray1d() );



  /// returns reference to specified field
  template< FieldKey FIELDKEY>
  typename ViewWrapper< Array1dT< typename Field<FIELDKEY>::Type > >::rtype GetFieldData( )
  {
    return const_cast<typename ViewWrapper< Array1dT< typename Field<FIELDKEY>::Type > >::rtype>( static_cast<const ManagedGroup&>(*this).GetFieldData<FIELDKEY>());
  }


  /// returns const reference to specified field
  template< FieldKey FIELDKEY>
  typename ViewWrapper< Array1dT< typename Field<FIELDKEY>::Type > >::rtype_const GetFieldData( ) const
  {
    return this->getData< Array1dT< typename Field<FIELDKEY>::Type >::rtype_const >( Field<FIELDKEY>::Name() );
  }


  /// returns reference to specified field
  template< typename TYPE >
  typename ViewWrapper< TYPE >::rtype GetFieldData( const std::string& fieldName )
  {
    return const_cast<typename ViewWrapper<TYPE>::rtype>( static_cast<const ManagedGroup&>(*this).GetFieldData<TYPE>(fieldName));
  }

  /// returns const reference to specified field
  template< typename TYPE >
  const Array1dT<TYPE>& GetFieldData( const std::string& name ) const
  {
    return this->getData< TYPE >( name );
  }








  /// returns reference to specified field
  template< FieldKey FIELDKEY>
  typename ViewWrapper< typename Field<FIELDKEY>::Type >::rtype* GetFieldDataPointer( )
  {
    return const_cast<typename ViewWrapper<typename Field<FIELDKEY>::Type>::rtype*>( static_cast<const ManagedGroup&>(*this).GetFieldDataPointer<FIELDKEY>());
  }


  /// returns const reference to specified field
  template< FieldKey FIELDKEY>
  typename ViewWrapper< typename Field<FIELDKEY>::Type >::rtype_const* GetFieldDataPointer( ) const
  {
    return this->getData< typename Field<FIELDKEY>::Type >( Field<FIELDKEY>::Name() );
  }

  /// returns reference to specified field
  template< typename TYPE >
  typename ViewWrapper< TYPE >::rtype* GetFieldDataPointer( const std::string& fieldName )
  {
    return this->getData< TYPE >( fieldName );
  }

  /// returns const reference to specified field
  template< typename TYPE >
  typename ViewWrapper< TYPE >::rtype_const* GetFieldDataPointer( const std::string& name ) const;



//**********************************************************************************************************************


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



template< typename T >
ViewWrapper<T>& ManagedGroup::RegisterViewWrapper( std::string const & name, std::size_t * const rkey )
{
  std::size_t key = static_cast<std::size_t>(-1);

  auto iterKeyLookup = m_keyLookup.find(name);

  // if the key was not found, make DataObject<T> and insert
  if( iterKeyLookup == m_keyLookup.end() )
  {
    m_wrappers.push_back( std::move( ViewWrapper<T>::Factory(name,this) ) );
    key = m_wrappers.size() - 1;
    m_keyLookup.insert( std::make_pair(name,key) );
    m_wrappers.back()->resize(this->size());
  }
  // if key was found, make sure that they are the same type
  else
  {
    key = m_keyLookup.at(name);
    auto& basePtr = m_wrappers[key];
    if( typeid(T) != basePtr->get_typeid() )
    {
      std::string error = string("Call to Group::RegisterViewWrapper( ")
                          +name+string(", std::size_t * const ) attempts to re-register ViewWrapper, but with different type") ;
//      SLIC_ERROR(error);
      throw std::exception();
    }
  }

  if( rkey != nullptr )
  {
    *rkey = key;
  }
  return getWrapper<T>(key);
}

template< typename T >
T& ManagedGroup::RegisterGroup( std::string const & name,
                                     std::unique_ptr<T> newObject )
{
  auto iterKeyLookup = m_subObjectManagers.find(name);

  // if the key was not found, make DataObject<T> and insert
  if( iterKeyLookup == m_subObjectManagers.end() )
  {
    auto insertResult = m_subObjectManagers.insert( std::make_pair( name, std::move(newObject) ) );

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

} // namespace dataRepository
} /* namespace geosx */


typedef geosx::dataRepository::ManagedGroup ObjectDataStructureBaseT;

#endif /* MANAGEDGROUP_H_ */
