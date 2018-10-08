/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

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

template< typename T, typename T_PTR=T*, typename KEY_TYPE=std::string, typename INDEX_TYPE = int >
class MappedVector
{
public:
  using key_type      = KEY_TYPE;
  using mapped_type   = T_PTR;

  using container              = std::vector<T_PTR>;
  using LookupMapType          = std::unordered_map<KEY_TYPE, INDEX_TYPE >;
  using iterMap                = std::map<KEY_TYPE, T * >;
  using const_iterMap          = std::map<KEY_TYPE, T const * >;
  using iterator               = typename iterMap::iterator;
  using reverse_iterator       = typename iterMap::reverse_iterator;
  using const_reverse_iterator = typename const_iterMap::const_reverse_iterator;
  using const_iterator         = typename const_iterMap::const_iterator;
  using value_type             = typename iterMap::value_type;
  using pointer                = typename iterMap::pointer;
  using const_pointer          = typename const_iterMap::const_pointer;
  using reference              = typename iterMap::reference;
  using const_reference        = typename const_iterMap::const_reference;
  using size_type              = typename iterMap::size_type;

  using KeyIndex = KeyIndexT<KEY_TYPE,INDEX_TYPE>;


  MappedVector() = default;
  ~MappedVector() = default;
  MappedVector( MappedVector const & source ) = default;
  MappedVector & operator=( MappedVector const & source ) = default;

  MappedVector( MappedVector && source ) = default;
  MappedVector & operator=( MappedVector && source ) = default;


  T * insert( KEY_TYPE const & keyName, T_PTR source );


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
    return ( index>-1 && index<static_cast<INDEX_TYPE>( m_values.size() ) ) ? const_cast<T const *>(&(*m_values[index])) : nullptr;
  }

  /**
   *
   * @param index
   * @return
   */
  inline T * operator[]( INDEX_TYPE index )
  { return const_cast<T*>( const_cast< MappedVector<T,T_PTR,KEY_TYPE,INDEX_TYPE> const * >(this)->operator[](index) ); }

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
  { return const_cast<T*>( const_cast< MappedVector<T,T_PTR,KEY_TYPE,INDEX_TYPE> const * >(this)->operator[](keyName) ); }

  /**
   *
   * @param dataKey
   * @return
   */
  inline T const * operator[]( KeyIndexT<KEY_TYPE,INDEX_TYPE> & dataKey ) const
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
  inline T * operator[]( KeyIndexT<KEY_TYPE,INDEX_TYPE> & dataKey )
  { return const_cast<T*>( const_cast< MappedVector<T,T_PTR,KEY_TYPE,INDEX_TYPE> const * >(this)->operator[](dataKey) ); }

  ///@}

  inline INDEX_TYPE getIndex( KEY_TYPE const & keyName ) const
  {
    typename LookupMapType::const_iterator iter = m_keyLookup.find(keyName);
    return ( iter!=m_keyLookup.end() ? iter->second : -1 );
  }


  void erase( INDEX_TYPE index )
  {
    m_values[index] = nullptr;
    return;
  }

  void erase( KEY_TYPE const & keyName )
  {
    typename LookupMapType::const_iterator iter = m_keyLookup.find(keyName);
    if( iter!=m_keyLookup.end() )
    {
      m_values[iter->second] = nullptr;
    }
    return;
  }

  void clear()
  {
    m_values.clear();
  }

  inline INDEX_TYPE size() const
  {
    return m_values.size();
  }

  std::vector<T_PTR> & values()
  {
    return m_values;
  }

  std::vector<T const *> const values() const
  {
    std::vector<T const *> rval(m_values.size());
    for( unsigned int a=0 ; a<m_values.size() ; ++a )
    {
      rval[a] = &(*m_values[a]);
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
  std::vector<T_PTR> m_values;

  LookupMapType m_keyLookup;

  std::vector< KEY_TYPE> m_indexToKeys;

  iterMap m_mapForIteration;
  const_iterMap m_mapForConstIteration;

};

template< typename T, typename T_PTR, typename KEY_TYPE, typename INDEX_TYPE >
T * MappedVector<T,T_PTR,KEY_TYPE,INDEX_TYPE>::insert( KEY_TYPE const & keyName, T_PTR source )
{
  typename LookupMapType::iterator iterKeyLookup = m_keyLookup.find(keyName);

  INDEX_TYPE key = -1;
  // if the key was not found, make DataObject<T> and insert
  if( iterKeyLookup == m_keyLookup.end() )
  {
    m_values.push_back( std::move( source ) );
    m_indexToKeys.push_back( keyName );

    key = m_values.size() - 1;
    m_keyLookup.insert( std::make_pair(keyName,key) );

    m_mapForIteration.insert({keyName,&(*(m_values[key]))});
    m_mapForConstIteration.insert({keyName,&(*(m_values[key]))});
  }
  // if key was found, make sure it is empty
  else
  {
    key = iterKeyLookup->second;
    if( m_values[key]==nullptr )
    {
      m_values[key] = std::move(source);

      m_mapForIteration.insert({keyName,&(*(m_values[key]))});
      m_mapForConstIteration.insert({keyName,&(*(m_values[key]))});

    }
    else
    {
      // error?
    }
  }

  return m_values[key].get();
}

#endif /* SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_MAPVECTORCONTAINER_HPP_ */
