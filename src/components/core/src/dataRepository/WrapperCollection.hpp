/**
 * @file DataObjectManager.h
 * @date created on Nov 21, 2014
 * @author Randolph R. Settgast
 */


#ifndef DATAOBJECTMANAGER_H_
#define DATAOBJECTMANAGER_H_

#include <iostream>

#include "Wrapper.hpp"
#include "ObjectCatalog.hpp"

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
class WrapperCollection
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
  explicit WrapperCollection( std::string const & name,
                              WrapperCollection * const parent );

  /**
   *
   */
  virtual ~WrapperCollection();

  /**
   *
   * @param source source WrapperCollection
   */
  WrapperCollection( WrapperCollection&& source );


  WrapperCollection() = delete;
  WrapperCollection( WrapperCollection const & source ) = delete;
  WrapperCollection& operator=( WrapperCollection const & ) = delete;
  WrapperCollection& operator=(WrapperCollection&&) = delete;

  ///@}


  using CatalogInterface = cxx_utilities::CatalogInterface< WrapperCollection, std::string const &, WrapperCollection * const >;
  static CatalogInterface::CatalogType& GetCatalog();


  virtual void Registration( dataRepository::WrapperCollection * const )
  {}

  virtual const std::type_info& get_typeid() const
  {
    return typeid(*this);
  }


  template< typename T = WrapperCollection >
  T& RegisterChildWrapperCollection( std::string const & name, std::unique_ptr<T> newObject );

  template< typename T = WrapperCollection >
  T& RegisterChildWrapperCollection( std::string const & name )
  {
    return RegisterChildWrapperCollection<T>( name, std::move(std::make_unique< T >( name, this )) );
  }


  template< typename T = WrapperCollection >
  T& GetChildWrapperCollection( std::string const & name )
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T&>( *(m_subObjectManagers.at(name)) );
#else
    return static_cast<T&>( *(m_subObjectManagers.at(name)) );
#endif
  }



template< typename T = WrapperCollection >
T const & GetChildWrapperCollection( std::string const & name ) const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T const &>( *(m_subObjectManagers.at(name)) );
#else
    return static_cast<T const &>( *(m_subObjectManagers.at(name)) );
#endif
  }


  template< typename T >
  Wrapper<T>& RegisterWrapper( std::string const & name, std::size_t * const rkey = nullptr );


  WrapperBase& RegisterWrapper( std::string const & name, rtTypes::TypeIDs const & type );


  //***********************************************************************************************

  template< typename T >
  Wrapper<T> const & getWrapper( std::size_t const index ) const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast< Wrapper<T> const & >( *(m_wrappers[index]) );
#else
    return static_cast< Wrapper<T> const & >( *(m_wrappers[index]) );
#endif
  }

  template< typename T >
  Wrapper<T> & getWrapper( std::size_t const index )
  {
    return const_cast<Wrapper<T>&>( const_cast< WrapperCollection const *>(this)->getWrapper<T>( index ) );
  }

  template< typename T >
  Wrapper<T> const & getWrapper( std::string const & name ) const
  {
    auto index = m_keyLookup.at(name);
    return getWrapper<T>(index);
  }

  template< typename T >
  Wrapper<T>& getWrapper( std::string const & name )
  { return const_cast<Wrapper<T>&>( const_cast<const WrapperCollection*>(this)->getWrapper<T>( name ) ); }



  template< typename T >
  typename Wrapper<T>::rtype_const getData( std::size_t const index ) const
  {
    return getWrapper<T>(index).data();
  }

  template< typename T >
  typename Wrapper<T>::rtype getData( std::size_t const index )
  {
    return const_cast<typename Wrapper<T>::rtype>( const_cast<const WrapperCollection*>(this)->getData<T>( index ) );
  }

  template< typename T >
  typename Wrapper<T>::rtype_const getData( std::string const & name ) const
  {
    auto index = m_keyLookup.at(name);
    return getData<T>( index );
  }

  template< typename T >
  typename Wrapper<T>::rtype getData( std::string const & name )
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
    return const_cast<T&>( const_cast<const WrapperCollection*>(this)->getReference<T>( index ) );
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

  void resize( std::size_t newsize );

  inline std::size_t size() const
  {
    return *(getData<std_size_t>(keys::size));
  }

  // compatibility functions
  inline std::size_t DataLengths() {return size();}
  template< typename T >
  int AddKeylessDataField( const std::string& name, const bool restart, const bool plot )
  {
    RegisterWrapper<T>(name);
    return 0;
  }

  template< typename TYPE >
  typename Wrapper< array<TYPE> >::rtype GetFieldData( const std::string& fieldName )
  {
    return ( this->getData< array<TYPE> >(fieldName) );
  }

//#include "Common/Common.h"
//  template< FieldKey FIELDKEY>
//  typename Wrapper<TYPE>::rtype GetFieldData( )
//  {
//    return GetFieldData< array<typename Field<FIELDKEY>::Type> >(Field<FIELDKEY>::Name());
//  }




  asctoolkit::sidre::DataGroup * getSidreGroup()              { return m_sidreGroup; }
  asctoolkit::sidre::DataGroup const * getSidreGroup() const  { return m_sidreGroup; }

private:
  std::unordered_map<std::string,std::size_t> m_keyLookup;
  std::vector< std::unique_ptr<WrapperBase> > m_wrappers;

  WrapperCollection* m_parent = nullptr;
  std::unordered_map< std::string, std::unique_ptr<WrapperCollection> > m_subObjectManagers;

  asctoolkit::sidre::DataGroup* m_sidreGroup;

#if NOCHARTOSTRING_KEYLOOKUP == 1

  template< typename T = WrapperCollection >
  T const & GetChildWrapperCollection( char const * ) const;

  template< typename T = WrapperCollection >
  T& GetChildWrapperCollection( char const * name );


  template< typename T >
  typename Wrapper<T>::rtype_const getData( char const * ) const;

  template< typename T >
  typename Wrapper<T>::rtype getData( char const * );

  template< typename T >
  T const & getReference( char const * ) const;

  template< typename T >
  T & getReference( char const * );


  template< typename T >
  Wrapper<T> const & getWrapper( char const * ) const;
  template< typename T >

  Wrapper<T>& getWrapper( char const * );


#endif
};



template< typename T >
Wrapper<T>& WrapperCollection::RegisterWrapper( std::string const & name, std::size_t * const rkey )
{
  std::size_t key = static_cast<std::size_t>(-1);

  auto iterKeyLookup = m_keyLookup.find(name);

  // if the key was not found, make DataObject<T> and insert
  if( iterKeyLookup == m_keyLookup.end() )
  {
    m_wrappers.push_back( std::move( Wrapper<T>::Factory(name,this) ) );
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
      std::cout<<LOCATION<<std::endl;
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
T& WrapperCollection::RegisterChildWrapperCollection( std::string const & name,
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

#endif /* DATAOBJECTMANAGER_H_ */
