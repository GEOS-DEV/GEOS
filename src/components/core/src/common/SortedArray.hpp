#ifndef SRC_COMMON_SORTEDARRAY
#define SRC_COMMON_SORTEDARRAY

#include <vector>             /* for std::vector */
#include <algorithm>          /* for std::binary_search, std::lower_bound */

template< typename T >
class SortedArray
{

public:
  typedef typename std::vector<T>::iterator iterator;
  typedef typename std::vector<T>::const_iterator const_iterator;
  typedef typename std::vector<T>::size_type size_type;


  SortedArray():
    m_data()
  {}

  template <typename InputIterator>
  SortedArray(InputIterator first, InputIterator last):
    m_data()
  { insert(first, last); }


  ~SortedArray()
  {}


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
    iterator it = find(value);
    if (it != end() && *it == value)
    {
      return false;
    }

    m_data.insert(it, value);
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
  { return std::lower_bound(begin(), end(), value); }


  const_iterator find(const T& value) const
  { return std::lower_bound(begin(), end(), value); }


  size_type count(const T& value) const
  { return std::binary_search(begin(), end(), value); }

private:

  std::vector< T > m_data;
};

#endif /* SRC_COMMON_SORTEDARRAY */