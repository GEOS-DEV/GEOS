/**
 * @file DataObjectManager.h
 * @date created on Nov 21, 2014
 * @author Randolph R. Settgast
 */


#ifndef DATAOBJECTMANAGER_H_
#define DATAOBJECTMANAGER_H_

#include <iostream>

#include "../intrinsic/DataObject.hpp"
//#include "CodingUtilities/ANSTexception.hpp"

/**
 * namespace to encapsulate functions in simulation tools
 */
namespace geosx
{

/**
 * @author Randolph R. Settgast
 *
 * class that encapsulates and manages a collection of DataObjects. Can be considered a "node" in a
 * hierarchy of managers that represent physical groupings of data.
 *
 */
class DataObjectManager
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
  explicit DataObjectManager( std::string const & name );

  /**
   *
   */
  virtual ~DataObjectManager();

  DataObjectManager( DataObjectManager const & source ) = delete;

  DataObjectManager( DataObjectManager&& source );


  DataObjectManager() = delete;
  DataObjectManager& operator=( DataObjectManager const & ) = delete;
  DataObjectManager& operator=(DataObjectManager&&) = delete;

  ///@}



  virtual const std::type_info& get_typeid() const
  {
    return typeid(*this);
  }


  template< typename T >
  T& RegisterChildDataObjectManager( std::string const & name );

  template< typename T >
  T& GetChildDataObjectManager( std::string const & name )
  {
    return *(m_subObjectManagers.at(name));
  }

  template< typename T >
  DataObjectBase * RegisterDataObject( std::string const & name, std::size_t * const rkey = nullptr );


  DataObjectBase * RegisterDataObject( std::string const & name, rtTypes::TypeIDs const & type )
  {
    return rtTypes::ApplyTypeLambda( type,
                                     [this, &name]( auto a ) -> DataObjectBase *
                                     {
                                       return this->RegisterDataObject<decltype(a)>(name);
                                     } );
  }


  //***********************************************************************************************

  template< typename T >
  T const & GetDataObject( std::size_t const index ) const
  {
    return m_dataObjects[index]->getObject<T>();
  }
  template< typename T >
  T& GetDataObject( std::size_t const index )
  {
    return const_cast<T&>( const_cast<const DataObjectManager*>(this)->GetDataObject<T>( index ) );
  }


  template< typename T >
  const DataObject<T>& GetDataObject( std::string const & name ) const
  {
    auto index = m_keyLookup.at(name);
    return m_dataObjects[index]->getObject<T>();
  }

  template< typename T >
  DataObject<T>& GetDataObject( std::string const & name )
  { return const_cast<DataObject<T>&>( const_cast<const DataObjectManager*>(this)->GetDataObject<T>( name ) ); }








  template< typename T >
  typename DataObject<T>::rtype_const GetDataObjectData( std::size_t const index ) const
  { return m_dataObjects[index]->getObjectData<T>(); }

  template< typename T >
  typename DataObject<T>::rtype GetDataObjectData( std::size_t const index )
  { return const_cast<T&>( const_cast<const DataObjectManager*>(this)->GetDataObjectData<T>( index ) ); }



  template< typename T >
  typename DataObject<T>::rtype_const GetDataObjectData( std::string const & name ) const
  {
    auto index = m_keyLookup.at(name);
    return m_dataObjects[index]->getObjectData<T>();
  }
  template< typename T >
  typename DataObject<T>::rtype GetDataObjectData( std::string const & name )
  {
    auto index = m_keyLookup.at(name);
    return m_dataObjects[index]->getObjectData<T>();
  }

//  { return static_cast<typename DataObject<T>::rtype>( static_cast<const DataObjectManager *>(this)->GetDataObjectData<T>( name ) ); }


  template< typename T >
  T& GetDataObjectManager( std::string const & name )
  {
    return *(m_subObjectManagers.at(name));
  }

  void resize( std::size_t const newsize );

  std::size_t size() const { return m_size; }

private:
  std::size_t m_size;
  std::string m_name{"name not set"};
  std::string m_path{"path not set"};
  std::unordered_map<std::string,std::size_t> m_keyLookup;
  std::vector< std::unique_ptr<DataObjectBase> > m_dataObjects;

  DataObjectManager* m_parent = nullptr;
  std::unordered_map< std::string, std::unique_ptr<DataObjectManager> > m_subObjectManagers;


};




template< typename T >
DataObjectBase * DataObjectManager::RegisterDataObject( std::string const & name, std::size_t * const rkey )
{
  std::size_t key = static_cast<std::size_t>(-1);

  auto iterKeyLookup = m_keyLookup.find(name);

  // if the key was not found, make DataObject<T> and insert
  if( iterKeyLookup == m_keyLookup.end() )
  {
    m_dataObjects.push_back( std::move( DataObjectBase::Factory<T>(name) ) );
    key = m_dataObjects.size() - 1;
    m_keyLookup.insert( std::make_pair(name,key) );
  }
  // if key was found, make sure that they are the same type
  else
  {
    key = m_keyLookup.at(name);
    auto& basePtr = m_dataObjects[key];
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
  return m_dataObjects[key].get();
}

template< typename T >
T& DataObjectManager::RegisterChildDataObjectManager( std::string const & name )
{
  auto iterKeyLookup = m_subObjectManagers.find(name);

  // if the key was not found, make DataObject<T> and insert
  if( iterKeyLookup == m_subObjectManagers.end() )
  {
    auto insertResult = m_subObjectManagers.insert( std::move(std::make_pair( name, std::move( std::make_unique< T >( name ) ) ) ) );

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
  return *(iterKeyLookup->second);
}

} /* namespace ODS */

#endif /* DATAOBJECTMANAGER_H_ */
