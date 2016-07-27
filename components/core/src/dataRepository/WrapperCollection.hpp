/**
 * @file DataObjectManager.h
 * @date created on Nov 21, 2014
 * @author Randolph R. Settgast
 */


#ifndef DATAOBJECTMANAGER_H_
#define DATAOBJECTMANAGER_H_

#include <iostream>

#include "Wrapper.hpp"
//#include "CodingUtilities/ANSTexception.hpp"

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



  virtual const std::type_info& get_typeid() const
  {
    return typeid(*this);
  }


  template< typename T >
  T& RegisterChildWrapperCollection( std::string const & name, std::unique_ptr<T> newObject );

  template< typename T >
  T& RegisterChildWrapperCollection( std::string const & name )
  {
    return RegisterChildWrapperCollection<T>( name, std::move(std::make_unique< T >( name, this )) );
  }


  template< typename T >
  T& GetChildDataObjectManager( std::string const & name )
  {
    return *(m_subObjectManagers.at(name));
  }

  template< typename T >
  Wrapper<T>& RegisterWrapper( std::string const & name, std::size_t * const rkey = nullptr );


  WrapperBase& RegisterWrapper( std::string const & name, rtTypes::TypeIDs const & type );


  //***********************************************************************************************

  template< typename T >
  Wrapper<T> const & getWrapper( std::size_t const index ) const
  {
#ifdef USE_DYNAMIC_CASTING
    return dyanmic_cast< Wrapper<T>* >(m_wrappers[index].get());
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
  typename Wrapper<T>::rtype_const getWrappedObjectData( std::size_t const index ) const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T*>(m_wrappers[index].get())->data();
#else
    return static_cast<Wrapper<T>*>(m_wrappers[index].get())->data();
#endif
  }

  template< typename T >
  typename Wrapper<T>::rtype getWrappedObjectData( std::size_t const index )
  {
    return const_cast<typename Wrapper<T>::rtype>( const_cast<const WrapperCollection*>(this)->getWrappedObjectData<T>( index ) );
  }



  template< typename T >
  typename Wrapper<T>::rtype_const getWrappedObjectData( std::string const & name ) const
  {
    auto index = m_keyLookup.at(name);
    return getWrappedObjectData<T>( index );
  }
  template< typename T >
  typename Wrapper<T>::rtype getWrappedObjectData( std::string const & name )
  {
    auto index = m_keyLookup.at(name);
    return getWrappedObjectData<T>( index );
  }

//  { return static_cast<typename DataObject<T>::rtype>( static_cast<const DataObjectManager *>(this)->GetDataObjectData<T>( name ) ); }


  template< typename T >
  T& GetDataObjectManager( std::string const & name )
  {
    return *(m_subObjectManagers.at(name));
  }

  void resize( std::size_t newsize );

  inline std::size_t size() const
  {
    return *(getWrappedObjectData<std_size_t>("size"));
  }


  asctoolkit::sidre::DataGroup * getSidreGroup()              { return m_sidreGroup; }
  asctoolkit::sidre::DataGroup const * getSidreGroup() const  { return m_sidreGroup; }

private:
  std::unordered_map<std::string,std::size_t> m_keyLookup;
  std::vector< std::unique_ptr<WrapperBase> > m_wrappers;

  WrapperCollection* m_parent = nullptr;
  std::unordered_map< std::string, std::unique_ptr<WrapperCollection> > m_subObjectManagers;

  asctoolkit::sidre::DataGroup* m_sidreGroup;


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
