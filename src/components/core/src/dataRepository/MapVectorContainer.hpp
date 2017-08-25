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
  using MapType = std::unordered_map<KEY_TYPE, INDEX_TYPE >;

//  template< typename U=T2,
//            bool isUniquePtr = std::is_same< U, std::unique_ptr<T> >::value >
//  struct const_T2
//  {
//    typedef std::unique_ptr<T const> type;
//  };
//
//  template<typename U>
//  struct const_T2<U,false>
//  {
//    typedef T const * type;
//  };



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


//  HAS_MEMBER_FUNCTION( get, T *,const,,)
//  template<class U = T2>
//  typename std::enable_if<has_memberfunction_get<U>::value, T const *>::type
  T const *
  operator[]( INDEX_TYPE index ) const
  {
    return ( index>-1 && index<static_cast<INDEX_TYPE>( m_objects.size() ) ) ? &(*m_objects[ index ]) : nullptr ;
  }

//  template<class U = T2>
//  typename std::enable_if<!has_memberfunction_get<U>::value, T const *>::type
//  operator[]( INDEX_TYPE index ) const
//  {
//    return ( index>-1 && index<static_cast<INDEX_TYPE>( m_objects.size() ) ) ? m_objects[ index ] : nullptr ;
//  }

  inline T * operator[]( INDEX_TYPE index )
  {
    return const_cast<T*>( const_cast< MapVectorContainer<T,T2,KEY_TYPE,INDEX_TYPE> const * >(this)->operator[](index) );
  }


  T const * operator[]( KEY_TYPE const & keyName ) const
  {
    typename MapType::const_iterator iter = m_keyLookup.find(keyName);
    return ( iter!=m_keyLookup.end() ? this->operator[]( iter->second ) : nullptr );
  }


  inline T * operator[]( KEY_TYPE const & keyName )
  {
    return const_cast<T*>( const_cast< MapVectorContainer<T,T2,KEY_TYPE,INDEX_TYPE> const * >(this)->operator[](keyName) );
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

  std::vector<T *> & objects()
  {
    return m_objects;
  }

  std::vector<T const *> const objects() const
  {
    std::vector<T const *> rval(m_objects.size());
    for( unsigned int a=0 ; a<m_objects.size() ; ++a )
    {
      rval[a] = &(*m_objects[a]);
    }
    return rval;
  }

  MapType const & keys() const
  {
    return m_keyLookup;
  }

//
//  template<class U = T2>
//  typename std::enable_if<has_memberfunction_get<U>::value, T *>::type
//  iterDeref( typename std::vector<T2>::iterator objIter )
//  {
//    return objIter->get();
//  }
//
//  template<class U = T2>
//  typename std::enable_if<!has_memberfunction_get<U>::value, T *>::type
//  iterDeref( typename std::vector<T2>::iterator objIter )
//  {
//    return *objIter;
//  }
//

  std::map<KEY_TYPE,T*> map()
  {
    std::map<KEY_TYPE,T*> rval;
    typename MapType::iterator keyIter = m_keyLookup.begin();
    typename std::vector<T2>::iterator objIter = m_objects.begin();

    for(  ; keyIter!=m_keyLookup.end() ; ++keyIter, ++objIter )
    {
      rval[keyIter->first] = &(*objIter);
    }
    return rval;
  }

  std::map<KEY_TYPE,T const*> map() const
  {
    std::map<KEY_TYPE,T const*> rval;
    typename MapType::iterator keyIter = m_keyLookup.begin();
    typename std::vector<T2>::const_iterator objIter = m_objects.begin();

    for(  ; keyIter!=m_keyLookup.end() ; ++keyIter, ++objIter )
    {
      rval[keyIter->first] = &(*objIter);
    }
    return rval;
  }

private:
  std::vector<T2> m_objects;

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
