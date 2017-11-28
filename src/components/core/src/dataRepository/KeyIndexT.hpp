/**
 * @file KeyIndexT.hpp
 * @date Aug 23, 2017
 * @authors settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_KEYINDEXT_HPP_
#define SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_KEYINDEXT_HPP_

#include <string>

/**
 * @class KeyIndexT
 * @tparam KEY_TYPE the key type
 * @tparam INDEX_TYPE the index type
 * @tparam INVALID_INDEX the value of an unset/invalid index
 * This class is used for fast index based lookups for a MappedVector type. It
 * is templated on a contains a KEY_TYPE, which is defaulted
 * to a string, an INDEX_TYPE that defaults to an int. The key is const, while
 * the index is set upon first use. The intent
 * is to use the index for lookups, and check the key to confirm the key is
 * correct.
 */
template< typename KEY_TYPE = std::string, typename INDEX_TYPE = int, int INVALID_INDEX = -1 >
class KeyIndexT
{
public:
  using key_type      = KEY_TYPE;
  using index_type    = INDEX_TYPE;


  constexpr static INDEX_TYPE invalid_index = INVALID_INDEX;

  /**
   * deleted default constructor
   */
  KeyIndexT() = delete;

  /**
   * constructor sets the value of m_name
   * @param name the key that defines the KeyIndex
   */
  KeyIndexT( KEY_TYPE const & key ):
    m_key(key),
    m_index(INVALID_INDEX)
  {}

  /// default copy constructor
  KeyIndexT( KeyIndexT const & ) = default;

  /// default copy assignment operator
  KeyIndexT & operator=( KeyIndexT const & ) = default;

  /// default move constructor
  KeyIndexT( KeyIndexT && ) = default;

  /// default move assignment operator
  KeyIndexT & operator=( KeyIndexT && ) = default;

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
   * @param index new value of the index
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

#endif /* SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_KEYINDEXT_HPP_ */
