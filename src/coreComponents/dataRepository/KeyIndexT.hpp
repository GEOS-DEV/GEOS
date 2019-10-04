/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file KeyIndexT.hpp
 */

#ifndef GEOSX_DATAREPOSITORY_KEYINDEXT_HPP_
#define GEOSX_DATAREPOSITORY_KEYINDEXT_HPP_

#include <string>
#include <ostream>

/**
 * @class KeyIndexT
 * @tparam KEY_TYPE the key type
 * @tparam INDEX_TYPE the index type
 * @tparam INVALID_INDEX the value of an unset/invalid index
 *
 * This class is used for fast index based lookups for a MappedVector type. It
 * is templated on a contains a KEY_TYPE, which is defaulted to a string, an
 * INDEX_TYPE that defaults to an int. The key is const, while the index is set
 * upon first use. The intent is to use the index for lookups, and check the
 * key to confirm the key is correct.
 */
template< typename KEY_TYPE = std::string,
          typename INDEX_TYPE = int,
          int INVALID_INDEX = -1 >
class KeyIndexT
{
public:
  /// the type used for the map key
  using key_type      = KEY_TYPE;

  /// the type used for the index
  using index_type    = INDEX_TYPE;

  /// the value of an invalid index
  constexpr static INDEX_TYPE invalid_index = INVALID_INDEX;

  /**
   * deleted default constructor
   */
  KeyIndexT() = delete;

  /**
   * constructor sets the value of m_name
   * @param[in] key the key that defines the KeyIndex
   */
  KeyIndexT( KEY_TYPE const & key ):
    m_key( key ),
    m_index( INVALID_INDEX )
  {}

  /// default copy constructor
  KeyIndexT( KeyIndexT const & ) = default;

  /// default copy assignment operator
  KeyIndexT & operator=( KeyIndexT const & ) = default;

  /// default move constructor
  KeyIndexT( KeyIndexT && ) = default;

  /// default move assignment operator
  KeyIndexT & operator=( KeyIndexT && ) = default;

  /**
   * @brief comparison equals operator
   * @param key the key to compare
   * @return true of m_key==key
   */
  bool operator==( KEY_TYPE const & key ) const
  { return m_key == key; }

  /// default destructor
  virtual ~KeyIndexT() {}

  /**
   * access for m_key
   * @return a const reference of the key
   */
  KEY_TYPE const & Key() const
  { return m_key; }

  /**
   * access for m_index
   * @return a const reference to the index
   */
  INDEX_TYPE const & Index() const
  { return m_index; }

  /**
   * check to see of the index has been set
   * @return true if the index has been set
   */
  bool isIndexSet() const
  {
    return m_index==INVALID_INDEX ? false : true;
  }

  /**
   * function to set the index
   * @param[in] index new value of the index
   */
  void setIndex( INDEX_TYPE const & index )
  {
    m_index = index;
  }

private:
  /// const key value
  KEY_TYPE const m_key;

  /// index value
  INDEX_TYPE m_index;
};

/**
 * @tparam KEY_TYPE the type of the key
 * @tparam INDEX_TYPE the type of the index
 * @tparam INVALID_INDEX the value of an invalid index
 * @param os the stream
 * @param key the key to output
 * @return a reference to the stream \p os
 */
template< typename KEY_TYPE, typename INDEX_TYPE, int INVALID_INDEX >
std::ostream & operator<<( std::ostream & os, const KeyIndexT< KEY_TYPE, INDEX_TYPE, INVALID_INDEX > key )
{
  os << key.Key();
  return os;
}

#endif /* GEOSX_DATAREPOSITORY_KEYINDEXT_HPP_ */
