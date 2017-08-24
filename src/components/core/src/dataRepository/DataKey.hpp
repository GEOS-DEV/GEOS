/*
 * DataKey.hpp
 *
 *  Created on: Aug 23, 2017
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_DATAKEY_HPP_
#define SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_DATAKEY_HPP_

template< typename T_KEY = std::string, typename T_INDEX = int >
class DataKeyT
{
public:
  DataKeyT() = delete;
  DataKeyT( T_KEY const & name ):
    m_key(name),
    m_index(-1)
  {}

  virtual ~DataKeyT() {}

  T_KEY const & Key() const
  {
    return m_key;
  }

  T_INDEX const & Index() const
  {
    return m_index;
  }

  bool isIndexSet() const
  {
    return m_index==-1 ? false : true;
  }

  void setIndex( T_INDEX const & index )
  {
    m_index = index;
  }

private:
  T_KEY const m_key;
  T_INDEX m_index;

};

//using DataKey = DataKeyT<>;
#endif /* SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_DATAKEY_HPP_ */
