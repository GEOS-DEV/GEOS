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

/**
 * @file MappedVector.hpp
 * @date Aug 23, 2017
 * @authors settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_MAPPEDVECTOR_HPP_
#define SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_MAPPEDVECTOR_HPP_

#include "common/Logger.hpp"
#include "KeyIndexT.hpp"
#include "SFINAE_Macros.hpp"

namespace geosx
{
/**
 * @class MappedVector
 * This class defines a stl-like container that stores values in an stl vector,
 * and has a map lookup table to
 * access the values by a key. It combines the random access performance of a
 * vector when the index is known,
 * the flexibility of a mapped key lookup O(n) if only the key is known. In
 * addition, a keyIndex can be used for lookup,
 * which will give similar performance to an index lookup after the first use of
 * a keyIndex.
 */
template< typename T,
          typename T_PTR=T*,
          typename KEY_TYPE=std::string,
          typename INDEX_TYPE = int >
class MappedVector
{
public:
  static_assert( std::is_same< T_PTR, T * >::value || std::is_same< T_PTR, std::unique_ptr<T> >::value,
                 "invalid second template argument for MappedVector<T,T_PTR,KEY_TYPE,INDEX_TYPE>. Allowable types are T * and std::unique_ptr<T>." );


  using key_type      = KEY_TYPE;
  using mapped_type   = T_PTR;

  using LookupMapType          = std::unordered_map<KEY_TYPE, INDEX_TYPE >;
  using value_type             = typename std::pair< KEY_TYPE const, T_PTR >;
  using const_value_type       = typename std::pair< KEY_TYPE const, T const * >;
  using valueContainer         = std::vector<value_type>;
  using const_valueContainer   = std::vector<const_value_type>;
  using pointer                = typename valueContainer::pointer;
  using const_pointer          = typename valueContainer::const_pointer;
  using reference              = typename valueContainer::reference;
  using const_reference        = typename valueContainer::const_reference;
  using size_type              = typename valueContainer::size_type;


  using iterator               = typename valueContainer::iterator;
  using const_iterator         = typename const_valueContainer::const_iterator;
  using reverse_iterator       = typename valueContainer::reverse_iterator;
  using const_reverse_iterator = typename const_valueContainer::const_reverse_iterator;

  using KeyIndex = KeyIndexT<KEY_TYPE const,INDEX_TYPE>;


  MappedVector() = default;


  ~MappedVector()
  {
    clear();
  }



  MappedVector( MappedVector const & source ) = default;
  MappedVector & operator=( MappedVector const & source ) = default;
  MappedVector( MappedVector && source ) = default;
  MappedVector & operator=( MappedVector && source ) = default;



  /**
   * @name element access functions
   */
  ///@{


  /**
   * @param index
   * @return pointer to const T
   */
  inline T const * operator[]( INDEX_TYPE index ) const
  {
    return ( index>KeyIndex::invalid_index &&
             index<static_cast<INDEX_TYPE>( m_values.size() ) ) ? const_cast<T const *>(&(*(m_values[index].second))) : nullptr;
  }

  /**
   *
   * @param index
   * @return pointer to T
   */
  inline T * operator[]( INDEX_TYPE index )
  { return const_cast<T*>( const_cast< MappedVector<T,T_PTR,KEY_TYPE,INDEX_TYPE> const * >(this)->operator[](index) ); }

  /**
   *
   * @param keyName
   * @return pointer to const T
   */
  inline T const * operator[]( KEY_TYPE const & keyName ) const
  {
    typename LookupMapType::const_iterator iter = m_keyLookup.find(keyName);
    return ( iter!=m_keyLookup.end() ? this->operator[]( iter->second ) : nullptr );
  }

  /**
   *
   * @param keyName
   * @return pointer to T
   */
  inline T * operator[]( KEY_TYPE const & keyName )
  { return const_cast<T*>( const_cast< MappedVector<T,T_PTR,KEY_TYPE,INDEX_TYPE> const * >(this)->operator[](keyName) ); }

  /**
   *
   * @param keyIndex
   * @return pointer to const T
   */
  inline T const * operator[]( KeyIndex & keyIndex ) const
  {
    INDEX_TYPE index = keyIndex.Index();

    if( index==KeyIndex::invalid_index )
    {
      index = getIndex( keyIndex.Key() );
      keyIndex.setIndex(index);
    }
#ifdef MAPPED_VECTOR_RANGE_CHECKING
    else if (m_values[index].first!=keyIndex.Key() )
    {
      index = getIndex( keyIndex.Key() );
      keyIndex.setIndex(index);
    }
#endif

    return this->operator[]( index );
  }

  /**
   *
   * @param keyIndex
   * @return pointer to T
   */
  inline T * operator[]( KeyIndex & keyIndex )
  { return const_cast<T*>( const_cast< MappedVector<T,T_PTR,KEY_TYPE,INDEX_TYPE> const * >(this)->operator[](keyIndex) ); }



  /**
   *
   * @param keyIndex
   * @return pointer to const T
   */
  inline T const * operator[]( KeyIndex const & keyIndex ) const
  {
//    INDEX_TYPE index = keyIndex.Index();
//
//    if( index==KeyIndex::invalid_index )
//    {
//      GEOS_ERROR("MappedVector::operator[]( KeyIndex const & keyIndex ):
// invalid key index passed as const into accessor function\n");
//    }
//#if RANGE_CHECKING==1
//    else if (m_values[index].first!=keyIndex.Key() )
//    {
//      GEOS_ERROR("MappedVector::operator[]( KeyIndex const & keyIndex ):
// inconsistent key passed as const into accessor function\n")
//    }
//#endif
//
//    return this->operator[]( index );
    return this->operator[]( const_cast<KeyIndex&>(keyIndex) );
  }

  /**
   *
   * @param keyIndex
   * @return pointer to T
   */
  inline T * operator[]( KeyIndex const & keyIndex )
  { return const_cast<T*>( const_cast< MappedVector<T,T_PTR,KEY_TYPE,INDEX_TYPE> const * >(this)->operator[](keyIndex) ); }


  ///@}


  /**
   * @name iterator functions
   */
  ///@{

  /**
   * @return  a read/write iterator that points to the first
   *  element in the in m_objects.
   */
  iterator begin()
  { return m_values.begin(); }

  /**
   * @return  a read-only iterator that points to the first
   *  element in the in m_objects.
   */
  const_iterator begin() const
  { return m_constValues.begin(); }

  /**
   * @return  a read-only iterator that points to the first
   *  element in the in m_objects.
   */
  const_iterator cbegin()
  { return m_constValues.begin(); }

  /**
   * @return  a read/write iterator that points to the last
   *  element in the in m_objects.
   */
  iterator end()
  { return m_values.end(); }

  /**
   * @return  a read-only iterator that points to the last
   *  element in the in m_objects.
   */
  const_iterator end() const
  { return m_constValues.end(); }

  /**
   * @return  a read-only iterator that points to the last
   *  element in the in m_objects.
   */
  const_iterator cend()
  { return m_constValues.end(); }
  ///@}


  /**
   *
   * @param key value of the key to use in the lookup
   * @return index associated with key
   */
  inline INDEX_TYPE getIndex( KEY_TYPE const & key ) const
  {
    typename LookupMapType::const_iterator iter = m_keyLookup.find(key);
    return ( iter!=m_keyLookup.end() ? iter->second : KeyIndex::invalid_index );
  }


  /**
   * @name modifier functions
   */
  ///@{

  /**
   *
   * @param keyName
   * @param source
   * @param overwrite
   * @return
   */
  T * insert( KEY_TYPE const & keyName, T_PTR source,
              bool takeOwnership,
              bool overwrite = false );



  /**
   * @brief  Remove element at given index
   * @tparam dummy parameter to allow use of enable_if for multiple definitions
   * based on type.
   * @param  index  index of element to remove.
   * @return  void
   *
   *  This function will call delete on the pointer at given index and set it to
   * nullptr.
   */
  template< typename U = T_PTR >
  typename std::enable_if< std::is_same< U, T * >::value, void >::type
  erase( INDEX_TYPE index )
  {
    if( m_ownsValues[index] )
    {
      delete m_values[index].second;
    }
    m_values[index].second = nullptr;
    m_constValues[index].second = nullptr;
    return;
  }

  /**
   * @brief  Remove element at given index
   * @tparam dummy parameter to allow use of enable_if for multiple definitions
   * based on type.
   * @param  index  index of element to remove.
   * @return  void
   *
   *  This function will the pointer at given index and set it to nullptr.
   */
  template< typename U = T_PTR >
  typename std::enable_if< std::is_same< U, std::unique_ptr<T> >::value, void >::type
  erase( INDEX_TYPE index )
  {
    m_values[index].second = nullptr;
    m_constValues[index].second = nullptr;
    return;
  }

  /**
   *  @brief  Remove element at given key
   *  @param  key  key of element to remove.
   *  @return  void
   *
   *  This function will set the element at the given key to nullptr.
   */
  void erase( KEY_TYPE const & key )
  {
    typename LookupMapType::const_iterator iter = m_keyLookup.find(key);
    if( iter!=m_keyLookup.end() )
    {
      erase(iter->second);
    }
    return;
  }

  /**
   *  @brief  Remove element at given key
   *  @param  key  key of element to remove.
   *  @return  void
   *
   *  This function will set the element at the given key to nullptr.
   */
  void erase( KeyIndex & keyIndex )
  {
    INDEX_TYPE index = keyIndex.Index();

    if( (index==KeyIndex::invalid_index) || (m_values[index].first!=keyIndex.Key()) )
    {
      index = getIndex( keyIndex.Key() );
      keyIndex.setIndex(index);
    }
    erase( index );
  }


  void clear()
  {
    for( typename valueContainer::size_type a=0 ; a<m_values.size() ; ++a )
    {
      //TODO this needs to be a safe conversion
      erase(static_cast<INDEX_TYPE>(a));
    }
    m_keyLookup.clear();
    m_values.clear();
  }


  ///@}


  inline INDEX_TYPE size() const
  {
    //TODO this needs to be a safe conversion
    return static_cast<INDEX_TYPE>(m_values.size());
  }

  inline valueContainer const & values()
  { return this->m_values; }

  inline const_valueContainer const & values() const
  { return this->m_constValues; }

  inline LookupMapType const & keys() const
  { return m_keyLookup; }


private:
  valueContainer m_values;
  const_valueContainer m_constValues;

  LookupMapType m_keyLookup;

  std::vector<int> m_ownsValues;
};

template< typename T, typename T_PTR, typename KEY_TYPE, typename INDEX_TYPE >
T * MappedVector<T,T_PTR,KEY_TYPE,INDEX_TYPE>::insert( KEY_TYPE const & keyName,
                                                                 T_PTR source,
                                                                 bool takeOwnership,
                                                                 bool overwrite )
{
  INDEX_TYPE index = KeyIndex::invalid_index;
  typename LookupMapType::iterator iterKeyLookup = m_keyLookup.find(keyName);


  // if the key was not found, make DataObject<T> and insert
  if( iterKeyLookup == m_keyLookup.end() )
  {
    value_type newEntry = std::make_pair( keyName, std::move( source ) );
    m_values.push_back( std::move( newEntry ) );
    //TODO this needs to be a safe conversion
    index = static_cast<INDEX_TYPE>(m_values.size()) - 1;
    m_ownsValues.resize( index + 1 );
    if( takeOwnership )
    {
      m_ownsValues[index] = true;
    }

    m_keyLookup.insert( std::make_pair(keyName,index) );
    m_constValues.push_back( std::make_pair( keyName, &(*(m_values[index].second)) ) );

  }
  // if key was found
  else
  {
    index = iterKeyLookup->second;

    if( takeOwnership )
    {
      m_ownsValues[index] = true;
    }

    // if value is empty, then move source into value slot
    if( m_values[index].second==nullptr )
    {
      m_values[index].second = std::move( source );
      m_constValues[index].second =  &(*(m_values[index].second));
    }
    else
    {
      if( overwrite )
      {
        erase(index);
        m_values[index].second = std::move( source );
        m_constValues[index].second =  &(*(m_values[index].second));
      }
      else if( source->get_typeid() != m_values[index].second->get_typeid() )
      {
        string const message = "MappedVector::insert(): Tried to insert existing key that was not empty without overwrite flag\n";
        GEOS_ERROR(message);
      }
      else
      {
        delete source;
      }
    }
  }

  return &(*(m_values[index].second));
}
}

#endif /* SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_MAPPEDVECTOR_HPP_ */
