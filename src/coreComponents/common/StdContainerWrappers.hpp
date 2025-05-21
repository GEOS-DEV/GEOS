#ifndef GEOS_COMMON_STD_CONTAINER_WRAPPERS_HPP
#define GEOS_COMMON_STD_CONTAINER_WRAPPERS_HPP

#include <vector>
#include <map>
#include <unordered_map>
#include <memory>

namespace geos
{

#ifdef GEOS_USE_BOUNDS_CHECK
#define USE_STD_CONTAINER_BOUNDS_CHECKING true
#else
#define USE_STD_CONTAINER_BOUNDS_CHECKING false
#endif


namespace internal
{

/**
 * Wrapper for std::vector that allows toggling between bounds-checked access
 * (using at()) and unchecked access (using operator[]).
 * @tparam T Type of elements in the vector.
 * @tparam USE_STD_CONTAINER_BOUNDS_CHECKING If true, uses at() for bounds-checked access.
 * If false, uses operator[] for unchecked access.
 */
template< typename T,
          bool USE_BOUNDS_CHECKING = false
          >
class StdVectorWrapper : public std::vector< T >
{
public:
  /// Type alias for the base class (i.e., std::vector)
  using Base = std::vector< T >;
  // Inherit constructors
  using Base::Base;

  /**
   * @brief Conversion constructor for StdVectorWrapper.
   * @tparam T Type of elements in the vector.
   * @param vector std::vector of elements to copy into the StdVectorWrapper.
   */
  StdVectorWrapper( std::vector< T > const & vec ): Base( std::move(vec) ){}

  /**
   * @brief Conversion constructor for StdVectorWrapper.
   * @tparam T Type of elements in the vector.
   * @param vector std::vector of elements to copy into the StdVectorWrapper.
   */
  StdVectorWrapper( std::vector< T > && vec ): Base( vec ){}


  /**
   * Access element at index with bounds checking if USE_STD_CONTAINER_BOUNDS_CHECKING is true.
   * Otherwise, uses operator[] for unchecked access.
   * @param index Index of the element to access.
   * @return Const reference to the element at the specified index.
   */
  T const & operator[]( size_t const index ) const
  {
    if constexpr (USE_BOUNDS_CHECKING)
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
    if constexpr (USE_BOUNDS_CHECKING)
    {
      return Base::at( index );  // Throws std::out_of_range if out of bounds
    }
    else
    {
      return Base::operator[]( index );  // No bounds checking
    }
  }
};
}

/**
 * type alias for std::vector
 * @tparam T Type of elements in the vector.
 */
template< typename T>
using stdVector = internal::StdVectorWrapper< T, USE_STD_CONTAINER_BOUNDS_CHECKING >;



namespace internal
{

/**
 * Wrapper for the underlying map that allows toggling between bounds-checked access
 * (using at()) and unchecked access (using operator[]).
 * @tparam MapType The type of the underlying map (e.g., std::map).
 * @tparam Allocator Allocator type for the vector.
 * @tparam USE_STD_CONTAINER_BOUNDS_CHECKING A boolean flag to enable or disable bounds checking.
 */
template< typename MapType,
          bool USE_BOUNDS_CHECKING = false
          >
class StdMapWrapper : public MapType
{
public:
  /// Type alias for the base class (i.e., std::map)
  using Base = MapType;
  /// Inherit constructors
  using Base::Base;
  using KeyType = typename Base::key_type;
  using MappedType = typename Base::mapped_type;
  using ValueType = typename Base::value_type;

  // Override operator[]
  MappedType & operator[]( KeyType const & key )
  {
    if constexpr (USE_BOUNDS_CHECKING)
    {
      return this->at( key );  // Throws std::out_of_range if key is missing
    }
    else
    {
      return Base::operator[]( key );  // Inserts default-constructed value if missing
    }
  }

  MappedType const & operator[]( KeyType const & key ) const
  {
    if constexpr (USE_BOUNDS_CHECKING)
    {
      return this->at( key );
    }
    else
    {
      return Base::operator[]( key );
    }
  }
};

} //namespace internal

/**
 * type alias for std::map
 * @tparam Key The unique std::map key.
 * @tparam T Type of elements in the std::map.
 * @tparam Compare The comparison function used to order the keys. Defaults to std::less<Key>.
 * @tparam Allocator Allocator type for the map. Defaults to std::allocator<std::pair<const Key, T>>
 */
template< typename Key,
          typename T,
          typename Compare = std::less< Key >,
          typename Allocator = std::allocator< std::pair< const Key, T > > >
using stdMap = internal::StdMapWrapper< std::map< Key, T, Compare, Allocator >,
                                        USE_STD_CONTAINER_BOUNDS_CHECKING >;

/**
 * type alias for std::unordered_map
 * @tparam Key The unique std::unordered_map key.
 * @tparam T Type of elements in the std::unordered_map.
 * @tparam Hash The hash function to be used for the keys. Defaults to std::hash<Key>
 * @tparam KeyEqual The function used to compare keys for equality. Defaults to std::equal_to<Key>.
 * @tparam Allocator Allocator type for the map.  Defaults to std::allocator<std::pair<const Key, T>>
 */
template< typename Key,
          typename T,
          typename Hash = std::hash< Key >,
          typename KeyEqual = std::equal_to< Key >,
          typename Allocator = std::allocator< std::pair< const Key, T > > >
using stdUnorderedMap = internal::StdMapWrapper< std::unordered_map< Key, T, Hash, KeyEqual, Allocator >,
                                                 USE_STD_CONTAINER_BOUNDS_CHECKING >;

/**
 * @name Ordered and unordered map types.
 */
///@{

/**
 * @brief Base template for ordered and unordered maps.
 * @tparam TKEY key type
 * @tparam TVAL value type
 * @tparam SORTED a bool indicating whether map is ordered
 */
template< typename TKEY, typename TVAL, typename SORTED >
class mapBase
{};

/// @cond DO_NOT_DOCUMENT
template< typename TKEY, typename TVAL >
class mapBase< TKEY, TVAL, std::integral_constant< bool, true > > : public stdMap< TKEY, TVAL >
{
public:
  using stdMap< TKEY, TVAL >::stdMap;
};

template< typename TKEY, typename TVAL >
class mapBase< TKEY, TVAL, std::integral_constant< bool, false > > : public stdUnorderedMap< TKEY, TVAL >
{
  using stdUnorderedMap< TKEY, TVAL >::stdUnorderedMap;
};
/// @endcond

} // namespace geos

#endif /* GEOS_COMMON_STD_CONTAINER_WRAPPERS_HPP */
