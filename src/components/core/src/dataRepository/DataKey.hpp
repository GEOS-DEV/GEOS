/*
 * DataKey.hpp
 *
 *  Created on: Aug 23, 2017
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_DATAKEY_HPP_
#define SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_DATAKEY_HPP_

#include <string>

template< typename KEY_TYPE = std::string, typename INDEX_TYPE = int, int INVALID_INDEX = -1 >
class DataKeyT
{
public:
  constexpr static int invalid_index = INVALID_INDEX;

  DataKeyT() = delete;
  DataKeyT( KEY_TYPE const & name ):
    m_key(name),
    m_index(INVALID_INDEX)
  {}

  virtual ~DataKeyT() {}

  KEY_TYPE const & Key() const
  {
    return m_key;
  }

  INDEX_TYPE const & Index() const
  {
    return m_index;
  }

  bool isIndexSet() const
  {
    return m_index==INVALID_INDEX ? false : true;
  }

  void setIndex( INDEX_TYPE const & index )
  {
    m_index = index;
  }

private:
  KEY_TYPE const m_key;
  INDEX_TYPE m_index;

};

//using DataKey = DataKeyT<>;
#endif /* SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_DATAKEY_HPP_ */
