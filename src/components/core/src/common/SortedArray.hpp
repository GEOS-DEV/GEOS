#ifndef SRC_COMMON_SORTEDARRAY
#define SRC_COMMON_SORTEDARRAY

#include <vector>             /* for std::vector */
#include "slic/slic.hpp"      /* for slic macros */

template< typename T >
class SortedArray
{

public:
  typedef typename std::vector<T>::iterator iterator;
  typedef typename std::vector<T>::const_iterator const_iterator;
  typedef typename std::vector<T>::size_type size_type;


  SortedArray()
  {
    // m_data.clear();
  }

  template <typename InputIterator>
  SortedArray(InputIterator first, InputIterator last)
  {
    insert(first, last);
  }


  ~SortedArray()
  {
    // m_data.clear();
  }


  T* data()
  { return m_data.data(); }


  const T* data() const
  { return m_data.data(); }


  iterator begin() 
  { return m_data.begin(); }


  const_iterator begin() const 
  { return m_data.begin(); }


  iterator end() 
  { return m_data.end(); }


  const_iterator end() const 
  { return m_data.end(); }


  bool empty() const
  { return m_data.empty(); }


  size_type size() const
  { return m_data.size(); }


  void clear()
  { m_data.clear(); }


  bool insert(const T& value)
  {
    const size_type min_index = 0;
    const size_type max_index = size() - 1;

    if (size() == 0)
    {
      m_data.push_back(value);
      return true;
    }

    const size_type cur_index = binary_search(value, min_index, max_index);
    const T& cur_value = m_data[cur_index];

    if (cur_value == cur_value)
    {
      return false;
    }
    else if (cur_value > cur_value)
    {
      m_data.insert(m_data.begin() + cur_index, cur_value);
    } 
    else 
    {
      m_data.insert(m_data.begin() + cur_index + 1, cur_value);
    }
    return true;
  }


  template <class InputIterator>
  void insert(InputIterator first, InputIterator last)
  {
    for (; first != last; first++)
    {
      insert(*first);
    }
  }


  iterator find(const T& value)
  {
    const size_type min_index = 0;
    const size_type max_index = size() - 1;

    if (size() == 0)
    {
      return end();
    }

    const size_type cur_index = binary_search(value, min_index, max_index);
    const T& cur_value = m_data[cur_index];

    return (cur_value == value)? begin() + cur_index : end(); 
  }


  size_type count(const T& value) const
  {
    const size_type min_index = 0;
    const size_type max_index = size() - 1;

    if (size() == 0)
    {
      return 0;
    }

    const size_type cur_index = binary_search(value, min_index, max_index);
    const T& cur_value = m_data[cur_index];

    return cur_value == value;
  }


  bool isSorted() const
  {
    const size_type num_items = m_data.size();
    for (size_type i = 0; i < num_items - 1; ++i) 
    { 
      if (m_data[i] >= m_data[i + 1])
      {
        return false;
      }
    }

    return true;
  }

private:

  size_type binary_search(const T& value, size_type min_index, size_type max_index) const
  {
    SLIC_ASSERT(isSorted());

    if (min_index == max_index)
    {
      return min_index;
    }

    size_type cur_index = (min_index + max_index) / 2;
    while (min_index != max_index)
    {
      const T& cur_value = m_data[cur_index];

      if (cur_value > value)
      {
        max_index = (cur_index == min_index)? min_index : min_index - 1;
      } 
      else if (cur_value < value) 
      {
        min_index = cur_index + 1;
      }
      else
      {
        return cur_index;
      }

      cur_index = (min_index + max_index) / 2;
    }

    return cur_index;
  }


  std::vector< T > m_data;
};

#endif /* SRC_COMMON_SORTEDARRAY */