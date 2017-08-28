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
  using LookupMapType = std::unordered_map<KEY_TYPE, INDEX_TYPE >;

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
    typename LookupMapType::const_iterator iter = m_keyLookup.find(keyName);
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

  std::vector<T2> & objects()
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

  LookupMapType const & keys() const
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
    typename LookupMapType::iterator keyIter = m_keyLookup.begin();
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
    typename LookupMapType::iterator keyIter = m_keyLookup.begin();
    typename std::vector<T2>::const_iterator objIter = m_objects.begin();

    for(  ; keyIter!=m_keyLookup.end() ; ++keyIter, ++objIter )
    {
      rval[keyIter->first] = &(*objIter);
    }

    std::map<string,int>::iterator iter;
    iter->first;
    (*iter).first;
    return rval;
  }


  struct iterator
  {
    typedef std::pair<KEY_TYPE,T *>  value_type;
    typedef std::pair<KEY_TYPE,T *> & reference;
    typedef std::pair<KEY_TYPE,T *> * pointer;

    iterator( typename LookupMapType::iterator lookupIter,
              std::vector<T2> * objectPtr ):
      m_lookupIter(lookupIter),
      m_objectPtr(objectPtr),
      m_pair({m_lookupIter->first, &(*((*m_objectPtr)[m_lookupIter->second]))})
    {

    }

    iterator( iterator const & source ):
      m_lookupIter(source.m_lookupIter),
      m_objectPtr(source.m_objectPtr),
      m_pair(source.m_pair)
    {}

    iterator( iterator && source ):
      m_lookupIter(std::move(source.m_lookupIter)),
      m_objectPtr(std::move(source.m_objectPtr)),
      m_pair(std::move(source.m_pair))
    {}

    iterator & operator=( iterator const & source )
    {
      m_lookupIter = source.m_lookupIter;
      m_objectPtr = source.m_objectPtr;
      m_pair = source.m_pair;

      return *this;
    }

    iterator & operator=( iterator && source )
    {
      m_lookupIter = std::move(source.m_lookupIter);
      m_objectPtr = std::move(source.m_objectPtr);
      m_pair = std::move(source.m_pair);
      return *this;
    }

    reference operator*()
    {
      return m_pair;
    }

    pointer operator->()
    {
      return m_pair;
    }

    iterator & operator++()
    {
      do
      {
        ++m_lookupIter;
      } while( ((*m_objectPtr)[m_lookupIter->second])==nullptr );
      m_pair.first = m_lookupIter->first;
      m_pair.second = &(*((*m_objectPtr)[m_lookupIter->second]));

      return *this;
    }

    bool
    operator==(const iterator& rhs) const _GLIBCXX_NOEXCEPT
    { return m_lookupIter == rhs.m_lookupIter; }

    bool
    operator!=(const iterator& rhs) const _GLIBCXX_NOEXCEPT
    { return m_lookupIter != rhs.m_lookupIter; }

    typename LookupMapType::iterator m_lookupIter;
    std::vector<T2> const *  m_objectPtr;
    std::pair<KEY_TYPE,T *> m_pair;
  };

  struct const_iterator
  {
    typedef std::pair<KEY_TYPE,T const *>  value_type;
    typedef std::pair<KEY_TYPE,T const *> const & reference;
    typedef std::pair<KEY_TYPE,T const *> const * pointer;

    const_iterator( typename LookupMapType::const_iterator lookupIter,
                    std::vector<T2> const * objectPtr ):
      m_lookupIter(lookupIter),
      m_objectPtr(objectPtr),
      m_pair({m_lookupIter->first, &(*((*m_objectPtr)[m_lookupIter->second]))})
    {

    }

    const_iterator( iterator const & source ):
      m_lookupIter(source.m_lookupIter),
      m_objectPtr(source.m_objectPtr),
      m_pair(source.m_pair)
    {}

    const_iterator( iterator && source ):
      m_lookupIter(std::move(source.m_lookupIter)),
      m_objectPtr(std::move(source.m_objectPtr)),
      m_pair(std::move(source.m_pair))
    {}

    const_iterator & operator=( iterator const & source )
    {
      m_lookupIter = source.m_lookupIter;
      m_objectPtr = source.m_objectPtr;
      m_pair = source.m_pair;

      return *this;
    }

    const_iterator & operator=( iterator && source )
    {
      m_lookupIter = std::move(source.m_lookupIter);
      m_objectPtr = std::move(source.m_objectPtr);
      m_pair = std::move(source.m_pair);
      return *this;
    }

    reference operator*()
    {
      return m_pair;
    }

    pointer operator->()
    {
      return m_pair;
    }

    const_iterator & operator++()
    {
      do
      {
        ++m_lookupIter;
      } while( ((*m_objectPtr)[m_lookupIter->second])==nullptr );
      m_pair.first = m_lookupIter->first;
      m_pair.second = &(*((*m_objectPtr)[m_lookupIter->second]));

      return *this;
    }

    bool
    operator==(const const_iterator& rhs) const _GLIBCXX_NOEXCEPT
    { return m_lookupIter == rhs.m_lookupIter; }

    bool
    operator!=(const const_iterator& rhs) const _GLIBCXX_NOEXCEPT
    { return m_lookupIter != rhs.m_lookupIter; }

    typename LookupMapType::const_iterator m_lookupIter;
    std::vector<T2> const *  m_objectPtr;
    std::pair<KEY_TYPE,T *> m_pair;
  };

  iterator begin()
  { return iterator( m_keyLookup.begin(), &m_objects); }

  const_iterator begin() const
  { return const_iterator( m_keyLookup.begin(), &m_objects); }

  iterator end()
  {
    return iterator( m_keyLookup.end(), &m_objects);
  }

  const_iterator end() const
  {
    return const_iterator( m_keyLookup.end(), &m_objects);
  }

private:
  std::vector<T2> m_objects;

  LookupMapType m_keyLookup;

  std::vector<KEY_TYPE> m_indexToKeys;
  std::map<string,int>::iterator junk;
};

template< typename T, typename T2, typename KEY_TYPE, typename INDEX_TYPE >
T * MapVectorContainer<T,T2,KEY_TYPE,INDEX_TYPE>::insert( KEY_TYPE const & keyName , T2 source )
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
