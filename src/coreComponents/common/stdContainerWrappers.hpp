#ifndef GEOS_COMMON_STD_CONTAINER_WRAPPERS_HPP
#define GEOS_COMMON_STD_CONTAINER_WRAPPERS_HPP

#include <vector>
#include <map>
#include <unordered_map>
#include <memory>

namespace geos
{


template< typename T, 
          typename Allocator = std::allocator<T>, 
          bool USE_STD_CONTAINER_BOUNDS_CHECKING = false 
        >
class StdVectorWrapper : public std::vector<T, Allocator>
{
public:
  using Base = std::vector<T, Allocator>;
  
  using Base::Base; // inherit constructors
  using Base::push_back;
  using Base::size;
  using Base::empty;

  // Override operator[] to toggle between at() and []
  T & operator[](size_t const index) 
  {
    if constexpr (USE_STD_CONTAINER_BOUNDS_CHECKING) 
    {
      return Base::at(index);  // Throws std::out_of_range if out of bounds
    } 
    else 
    {
      return Base::operator[](index);  // No bounds checking
    }
  }

  T const & operator[](size_t const index) const 
  {
    if constexpr (USE_STD_CONTAINER_BOUNDS_CHECKING) 
    {
      return Base::at(index);
    } 
    else
    {
      return Base::operator[](index);
    }
  }
};




template< typename Key, 
          typename Value, 
          template <typename, typename, typename> typename MapType, 
          typename CompareOrHash = std::less<Key>, 
          typename Allocator = std::allocator<std::pair<const Key, Value>>,
          bool USE_STD_CONTAINER_BOUNDS_CHECKING = false
        >
class MapWrapper : public MapType<Key, Value, CompareOrHash, Allocator>
{
public:
  using Base = MapType<Key, Value, CompareOrHash, Allocator>;
  using Base::Base;  // Inherit constructors

  // Override operator[]
  Value & operator[]( Key const & key) 
  {
    if constexpr(USE_STD_CONTAINER_BOUNDS_CHECKING) 
    {
      return this->at(key);  // Throws std::out_of_range if key is missing
    }
    else 
    {
      return Base::operator[](key);  // Inserts default-constructed value if missing
    }
  }

  Value const & operator[]( Key const & key) const 
  {
    if constexpr(USE_STD_CONTAINER_BOUNDS_CHECKING) 
    {
      return this->at(key);
    }
    else 
    {
      return Base::operator[](key);
    }
  }
};


}

#endif /* GEOS_COMMON_STD_CONTAINER_WRAPPERS_HPP */