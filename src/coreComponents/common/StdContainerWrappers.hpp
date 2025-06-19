#ifndef GEOS_COMMON_STD_CONTAINER_WRAPPERS_HPP
#define GEOS_COMMON_STD_CONTAINER_WRAPPERS_HPP

#include <vector>
#include <map>
#include <unordered_map>
#include <memory>

namespace geos
{

namespace internal
{

/**
 * Wrapper for std::vector that allows toggling between bounds-checked access
 * (using at()) and unchecked access (using operator[]).
 * @tparam T Type of elements in the vector.
 * @tparam Allocator Allocator type for the vector.
 * @tparam USE_STD_CONTAINER_BOUNDS_CHECKING If true, uses at() for bounds-checked access.
 * If false, uses operator[] for unchecked access.
 */
template< typename T,
          typename Allocator = std::allocator< T >,
          bool USE_STD_CONTAINER_BOUNDS_CHECKING = false
          >
class StdVectorWrapper : public std::vector< T, Allocator >
{
public:
  /// Type alias for the base class (i.e., std::vector)
  using Base = std::vector< T, Allocator >;

  /**
   * Access element at index with bounds checking if USE_STD_CONTAINER_BOUNDS_CHECKING is true.
   * Otherwise, uses operator[] for unchecked access.
   * @param index Index of the element to access.
   * @return Const reference to the element at the specified index.
   */
  T const & operator[]( size_t const index ) const
  {
    if constexpr (USE_STD_CONTAINER_BOUNDS_CHECKING)
    {
      return Base::at( index );
    }
    else
    {
      return Base::operator[]( index );
    }
  }

  /**
   * Access element at index with bounds checking if USE_STD_CONTAINER_BOUNDS_CHECKING is true.
   * Otherwise, uses operator[] for unchecked access.
   * @param index Index of the element to access.
   * @return Reference to the element at the specified index.
   */
  T & operator[]( size_t const index )
  {
    if constexpr (USE_STD_CONTAINER_BOUNDS_CHECKING)
    {
      return Base::at( index );  // Throws std::out_of_range if out of bounds
    }
    else
    {
      return Base::operator[]( index );  // No bounds checking
    }
  }
};
} // namespace internal

#if defined( GEOS_USE_BOUNDS_CHECK )
/**
 * type alias for internal::StdVectorWrapper with bounds checking enabled.
 * @tparam T Type of elements in the vector.
 * @tparam Allocator Allocator type for the vector.
 */
template< typename T, typename Allocator = std::allocator< T > >
using stdVector = internal::StdVectorWrapper< T, Allocator, true >;
#else
/**
 * type alias for std::vector with no bounds checking.
 * @tparam T Type of elements in the vector.
 * @tparam Allocator Allocator type for the vector.
 */
template< typename T, typename Allocator = std::allocator< T > >
using stdVector = std::vector< T, Allocator >;
#endif

// template< typename MapType,
//           bool USE_STD_CONTAINER_BOUNDS_CHECKING >
// class StdMapWrapper : public MapType
// {
// public:
//   using Base = MapType;
//   using Base::Base;  // Inherit constructors
//   using KeyType = typename Base::key_type;
//   using MappedType = typename Base::mapped_type;
//   using ValueType = typename Base::value_type;

//   // Override operator[]
//   MappedType & operator[]( KeyType const & key)
//   {
//     if constexpr(USE_STD_CONTAINER_BOUNDS_CHECKING)
//     {
//       return this->at(key);  // Throws std::out_of_range if key is missing
//     }
//     else
//     {
//       return Base::operator[](key);  // Inserts default-constructed value if missing
//     }
//   }

//   MappedType const & operator[]( KeyType const & key) const
//   {
//     if constexpr(USE_STD_CONTAINER_BOUNDS_CHECKING)
//     {
//       return this->at(key);
//     }
//     else
//     {
//       return Base::operator[](key);
//     }
//   }
// };

//} //namespace internal

// template< typename Key,
//           typename T,
//           typename Compare = std::less<Key>,
//           typename Allocator = std::allocator<std::pair<const Key, T>>,
//           bool USE_STD_CONTAINER_BOUNDS_CHECKING = false>
// using map = internal::StdMapWrapper< std::map< Key, T, Compare, Allocator> ,
//                                      USE_STD_CONTAINER_BOUNDS_CHECKING >;


// template< typename Key,
//           typename T,
//           typename Hash = std::hash<Key>,
//           typename KeyEqual = std::equal_to<Key>,
//           typename Allocator = std::allocator<std::pair<const Key, T>>,
//           bool USE_STD_CONTAINER_BOUNDS_CHECKING = false>
// using unordered_map = internal::StdMapWrapper< std::unordered_map< Key, T, Hash, KeyEqual, Allocator >,
//                                                USE_STD_CONTAINER_BOUNDS_CHECKING >;


} // namespace geos

#endif /* GEOS_COMMON_STD_CONTAINER_WRAPPERS_HPP */
