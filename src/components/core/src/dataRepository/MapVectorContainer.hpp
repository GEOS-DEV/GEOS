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
class MapVectorContainer
{
public:
  MapVectorContainer() = default;
  ~MapVectorContainer() = default;
  MapVectorContainer( MapVectorContainer const & source ) = default ;
  MapVectorContainer & operator=( MapVectorContainer const & source ) = default;

  MapVectorContainer( MapVectorContainer && source ) = default;
  MapVectorContainer & operator=( MapVectorContainer && source ) = default;


  T * insert( KEY_TYPE const & keyName, T2 source );



  /**
   * @name accessor functions
   */
  ///@{


  T const * operator[]( KEY_TYPE const & keyName ) const
  {
    typename MapType::const_iterator iter = m_keyLookup.find(keyName);
    return ( iter!=m_keyLookup.end() ? this->operator[]( iter->second ) : nullptr );
  }


  inline T * operator[]( KEY_TYPE const & keyName )
  {
    return const_cast<T*>( const_cast< MapVectorContainer<T,T2,KEY_TYPE,INDEX_TYPE> const * >(this)->operator[](keyName) );
  }


  inline T * operator[]( INDEX_TYPE index )
  {
    return const_cast<T*>( const_cast< MapVectorContainer<T,T2,KEY_TYPE,INDEX_TYPE> const * >(this)->operator[](index) );
  }

  HAS_MEMBER_FUNCTION( get, T *,const,,)
  template<class U = T2>
  typename std::enable_if<has_memberfunction_get<U>::value, T const *>::type
  operator[]( INDEX_TYPE index ) const
  {
    return m_objects[ index ].get();
  }
  template<class U = T2>
  typename std::enable_if<!has_memberfunction_get<U>::value, T const *>::type
  operator[]( INDEX_TYPE index ) const
  {
    return m_objects[ index ];
  }

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

  inline T * operator[]( DataKeyT<KEY_TYPE,INDEX_TYPE> & dataKey )
  {
    return const_cast<T*>( const_cast< MapVectorContainer<T,T2,KEY_TYPE,INDEX_TYPE> const * >(this)->operator[](dataKey) );
  }

  ///@}

  inline INDEX_TYPE getIndex( KEY_TYPE const & keyName ) const
  {
    typename MapType::const_iterator iter = m_keyLookup.find(keyName);
    return ( iter!=m_keyLookup.end() ? iter->second : -1 );
  }


  void erase( INDEX_TYPE index )
  {
    m_objects[index] = nullptr;
    return;
  }

  void erase( KEY_TYPE const & keyName )
  {
    typename MapType::const_iterator iter = m_keyLookup.find(keyName);
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

  std::vector<T2> & objects()
  {
    return m_objects;
  }
  std::vector<T2> const & objects() const
  {
    return m_objects;
  }

private:
  std::vector<T2> m_objects;

  using MapType = std::unordered_map<KEY_TYPE, INDEX_TYPE >;
  MapType m_keyLookup;

  std::vector<KEY_TYPE> m_indexToKeys;


};

template< typename T, typename T2, typename KEY_TYPE, typename INDEX_TYPE >
T * MapVectorContainer<T,T2,KEY_TYPE,INDEX_TYPE>::insert( KEY_TYPE const & keyName , T2 source )
{
  typename MapType::iterator iterKeyLookup = m_keyLookup.find(keyName);

  INDEX_TYPE key = -1;
  // if the key was not found, make DataObject<T> and insert
  if( iterKeyLookup == m_keyLookup.end() )
  {
    m_objects.push_back( std::move( source ) );
    m_indexToKeys.push_back( keyName );

    key = m_objects.size() - 1;
    m_keyLookup.insert( std::make_pair(keyName,key) );
    if( m_objects[key]->sizedFromParent() == 1 )
    {
      m_objects[key]->resize(this->size());
    }
  }
  // if key was found, make sure it is empty
  else
  {
    key = iterKeyLookup->second;
    if( m_objects[key]==nullptr )
    {
      m_objects[key] = std::move(source);
    }
  }

return m_objects[key].get();
}

#endif /* SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_MAPVECTORCONTAINER_HPP_ */
