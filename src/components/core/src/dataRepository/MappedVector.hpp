/*
 * MapVectorContainer.hpp
 *
 *  Created on: Aug 23, 2017
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_MAPVECTORCONTAINER_HPP_
#define SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_MAPVECTORCONTAINER_HPP_

#include "DataKey.hpp"
#include "SFINAE_Macros.hpp"

template< typename T, typename T2=T*, typename KEY_TYPE=std::string, typename INDEX_TYPE = int >
class MappedVector
{
public:
  using LookupMapType = std::unordered_map<KEY_TYPE, INDEX_TYPE >;
  using iterMap = typename std::map<KEY_TYPE, T * >;
  using const_iterMap = typename std::map<KEY_TYPE, T const * >;
  using iterator = typename iterMap::iterator;
  using const_iterator = typename const_iterMap::const_iterator;

  MappedVector() = default;
  ~MappedVector() = default;
  MappedVector( MappedVector const & source ) = default ;
  MappedVector & operator=( MappedVector const & source ) = default;

  MappedVector( MappedVector && source ) = default;
  MappedVector & operator=( MappedVector && source ) = default;


  T * insert( KEY_TYPE const & keyName, T2 source );


  /**
   * @name accessor functions
   */
  ///@{


  /**
   *
   * @param index
   * @return
   */
  T const * operator[]( INDEX_TYPE index ) const
  {
    return ( index>-1 && index<static_cast<INDEX_TYPE>( m_objects.size() ) ) ? const_cast<T const *>(&(*m_objects[index])) : nullptr ;
  }

  /**
   *
   * @param index
   * @return
   */
  inline T * operator[]( INDEX_TYPE index )
  { return const_cast<T*>( const_cast< MappedVector<T,T2,KEY_TYPE,INDEX_TYPE> const * >(this)->operator[](index) ); }

  /**
   *
   * @param keyName
   * @return
   */
  T const * operator[]( KEY_TYPE const & keyName ) const
  {
    typename LookupMapType::const_iterator iter = m_keyLookup.find(keyName);
    return ( iter!=m_keyLookup.end() ? this->operator[]( iter->second ) : nullptr );
  }

  /**
   *
   * @param keyName
   * @return
   */
  inline T * operator[]( KEY_TYPE const & keyName )
  { return const_cast<T*>( const_cast< MappedVector<T,T2,KEY_TYPE,INDEX_TYPE> const * >(this)->operator[](keyName) ); }

  /**
   *
   * @param dataKey
   * @return
   */
  inline T const * operator[]( DataKeyT<KEY_TYPE,INDEX_TYPE> & dataKey ) const
  {
    INDEX_TYPE index = dataKey.Index();

    if( (index==-1) || (m_indexToKeys[index]!=dataKey.Key()) )
    {
      index = getIndex( dataKey.Key() );
      dataKey.setIndex(index);
    }

    return this->operator[]( index );
  }

  /**
   *
   * @param dataKey
   * @return
   */
  inline T * operator[]( DataKeyT<KEY_TYPE,INDEX_TYPE> & dataKey )
  { return const_cast<T*>( const_cast< MappedVector<T,T2,KEY_TYPE,INDEX_TYPE> const * >(this)->operator[](dataKey) ); }

  ///@}

  inline INDEX_TYPE getIndex( KEY_TYPE const & keyName ) const
  {
    typename LookupMapType::const_iterator iter = m_keyLookup.find(keyName);
    return ( iter!=m_keyLookup.end() ? iter->second : -1 );
  }


  void erase( INDEX_TYPE index )
  {
    m_objects[index] = nullptr;
    return;
  }

  void erase( KEY_TYPE const & keyName )
  {
    typename LookupMapType::const_iterator iter = m_keyLookup.find(keyName);
    if( iter!=m_keyLookup.end() )
    {
      m_objects[iter->second] = nullptr;
    }
    return;
  }

  void clear()
  {
    m_objects.clear();
  }

  inline INDEX_TYPE size() const
  {
    return m_objects.size();
  }

  std::vector<T2> & values()
  {
    return m_objects;
  }

  std::vector<T const *> const values() const
  {
    std::vector<T const *> rval(m_objects.size());
    for( unsigned int a=0 ; a<m_objects.size() ; ++a )
    {
      rval[a] = &(*m_objects[a]);
    }
    return rval;
  }

  LookupMapType const & keys() const
  {
    return m_keyLookup;
  }

  iterator begin()
  { return m_mapForIteration.begin(); }

  const_iterator begin() const
  { return m_mapForConstIteration.begin(); }

  iterator end()
  {
    return m_mapForIteration.end();
  }

  const_iterator end() const
  {
    return m_mapForConstIteration.end();
  }

private:
  std::vector<T2> m_objects;

  LookupMapType m_keyLookup;

  std::vector<KEY_TYPE> m_indexToKeys;

  iterMap m_mapForIteration;
  const_iterMap m_mapForConstIteration;

};

template< typename T, typename T2, typename KEY_TYPE, typename INDEX_TYPE >
T * MappedVector<T,T2,KEY_TYPE,INDEX_TYPE>::insert( KEY_TYPE const & keyName , T2 source )
{
  typename LookupMapType::iterator iterKeyLookup = m_keyLookup.find(keyName);

  INDEX_TYPE key = -1;
  // if the key was not found, make DataObject<T> and insert
  if( iterKeyLookup == m_keyLookup.end() )
  {
    m_objects.push_back( std::move( source ) );
    m_indexToKeys.push_back( keyName );

    key = m_objects.size() - 1;
    m_keyLookup.insert( std::make_pair(keyName,key) );

//    std::cout<<&(*m_objects[key])<<std::endl;
    m_mapForIteration.insert({keyName,&(*(m_objects[key]))});
    m_mapForConstIteration.insert({keyName,&(*(m_objects[key]))});
  }
  // if key was found, make sure it is empty
  else
  {
    key = iterKeyLookup->second;
    if( m_objects[key]==nullptr )
    {
      m_objects[key] = std::move(source);

      m_mapForIteration.insert({keyName,&(*(m_objects[key]))});
      m_mapForConstIteration.insert({keyName,&(*(m_objects[key]))});

    }
    else
    {
      // error?
    }
  }

return m_objects[key].get();
}

#endif /* SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_MAPVECTORCONTAINER_HPP_ */
