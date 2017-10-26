#ifndef SRC_COMMON_SORTEDARRAY
#define SRC_COMMON_SORTEDARRAY

#include <vector>             /* for std::vector */
#include <algorithm>          /* for std::binary_search, std::lower_bound */

template< typename T >
class SortedArray
{

public:
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;
  using size_type = typename std::vector<T>::size_type;
  using value_type = typename std::vector<T>::value_type;
  using pointer = typename std::vector<T>::pointer;
  using const_pointer = typename std::vector<T>::const_pointer;


  SortedArray():
    m_data()
  {}

  template <typename InputIterator>
  SortedArray(InputIterator first, InputIterator last):
    m_data()
  { insert(first, last); }


  ~SortedArray()
  {}


  T * data()
  { return m_data.data(); }


  const T * data() const
  { return m_data.data(); }


  T & operator[](size_type i)
  { return m_data[i]; }

  const T & operator[](size_type i) const
  { return m_data[i]; }


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


  void resize(size_type new_size)
  { return m_data.resize(new_size); }


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


  bool isSorted() const
  { return std::is_sorted(begin(), end()); }

private:

  std::vector< T > m_data;
};

#endif /* SRC_COMMON_SORTEDARRAY */