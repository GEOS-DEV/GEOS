#ifndef GEOS_COMMON_STD_CONTAINER_WRAPPERS_HPP
#define GEOS_COMMON_STD_CONTAINER_WRAPPERS_HPP

#include <array>
#include <cstddef>
#include <array>
#include <cstddef>
#include <vector>
#include <map>
#include <unordered_map>
#include <memory>

#include <iostream>

namespace geos
{

/**
 * @def USE_STD_CONTAINER_BOUNDS_CHECKING
 * @brief Compile-time flag that enables or disables runtime bounds checking in GEOS container wrappers.
 */
#ifdef GEOS_USE_BOUNDS_CHECK
#define USE_STD_CONTAINER_BOUNDS_CHECKING true
#else
#define USE_STD_CONTAINER_BOUNDS_CHECKING false
#endif


/**
 * @file StdContainerWrappers.hpp
 *
 * @warning:
 * It is prohibited to use `std::map` or `std::unordered_map` in the geos repos.
 *
 * @section Usage
 * There are two ways to declare a map or unordered_map:
 *
 * 1. **stdMap** / **stdUnorderedMap**:
 *    - These types replace `std::map` with an overridden `operator[]` for bounds checking.
 *    - We cannot use the operator[] for the insertion.
 *
 * 2. **geos::map** / **geos::unordered_map**:
 *    - Use these types to ensure the compatibility with the geos packing, as `stdMap` is not yet compatible.
 */
namespace internal
{



/**
 * Default allocator type for std::vector.
 * This can be specialized if a different allocator is needed.
 * Required to avoid recursive evaluation in StdVectorWrapper.
 * @tparam T Type of elements in the vector.
 */
template< typename T >
using DefaultAllocator = std::allocator< T >;

/**
 * Wrapper for std::vector that allows toggling between bounds-checked access
 * (using at()) and unchecked access (using operator[]).
 * @tparam T Type of elements in the vector.
 * @tparam Allocator Allocator type for the vector.
 * @tparam USE_STD_CONTAINER_BOUNDS_CHECKING If true, uses at() for bounds-checked access.
 * If false, uses operator[] for unchecked access.
 */
template< typename T,
          typename Allocator  = DefaultAllocator< T >,
          bool USE_BOUNDS_CHECKING = false >
class StdVectorWrapper : public std::vector< T, Allocator >
{
public:
  /// Type alias for the base class (i.e., std::vector)
  using Base = std::vector< T, Allocator >;

  /// We have to declare explicitly all constructor because `using Base::Base` causes bugs during compilation
  /// @cond DO_NOT_DOCUMENT
  StdVectorWrapper() noexcept(noexcept(Allocator())): Base( Allocator()) {}

  explicit StdVectorWrapper( const Allocator & alloc ) noexcept: Base( alloc ) {}

  explicit StdVectorWrapper( size_t count, const Allocator & alloc = Allocator())
    : Base( count, alloc ) {}

  explicit StdVectorWrapper( size_t count, const T & value,
                             const Allocator & alloc = Allocator())
    : Base( count, value, alloc ) {}

  template< class InputIt >
  StdVectorWrapper( InputIt first, InputIt last,
                    const Allocator & alloc = Allocator())
    : Base( first, last, alloc ) {}

  StdVectorWrapper( const Base & other ): Base( other ) {}

  StdVectorWrapper( Base && other ): Base( other ) {}

  StdVectorWrapper( const Base & other, const Allocator & alloc )
    : Base( other, alloc ) {}

  StdVectorWrapper( Base && other, const Allocator & alloc )
    : Base( other, alloc ) {}

  StdVectorWrapper( std::initializer_list< T > init,
                    const Allocator & alloc = Allocator())
    : Base( init, alloc ) {}

  /// @endcond


  /**
   * Access element at index with bounds checking if USE_STD_CONTAINER_BOUNDS_CHECKING is true.
   * Otherwise, uses operator[] for unchecked access.
   * @param index Index of the element to access.
   * @return Const reference to the element at the specified index.
   * @throws std::out_of_range if index is out of bounds.
   */
  T const & operator[]( size_t const index ) const
  {
    if constexpr (USE_BOUNDS_CHECKING)
    {
      if( index >= this->size() )
      {
        std::cout<< "Index out of bounds in StdVectorWrapper::operator[]: index = " + std::to_string( index ) + ", size = " + std::to_string( this->size());
      }
      return Base::at( index );  // Throws std::out_of_range if out of bounds.
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
      if( index >= this->size() )
      {
        std::cout<< "Index out of bounds in StdVectorWrapper::operator[]: index = " + std::to_string( index ) + ", size = " + std::to_string( this->size());
      }
      return Base::at( index );  // Throws std::out_of_range if out of bounds.
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
 * @tparam Allocator Allocator type for the vector.
 */
template< typename T, typename Allocator = std::allocator< T > >
using stdVector = internal::StdVectorWrapper< T, Allocator, USE_STD_CONTAINER_BOUNDS_CHECKING >;



namespace internal
{

/**
 * Wrapper for the underlying map that allows toggling between bounds-checked access
 * (using at()) and unchecked access (using operator[]).
 * @tparam MapType The type of the underlying map (e.g., stdMap).
 * @tparam Allocator Allocator type for the vector.
 * @tparam USE_STD_CONTAINER_BOUNDS_CHECKING A boolean flag to enable or disable bounds checking.
 */
template< typename MapType,
          bool USE_BOUNDS_CHECKING = false >
class StdMapWrapper : public MapType
{
public:
  /// Type alias for the base class (i.e., stdMap)
  using Base = MapType;
  using KeyType = typename Base::key_type;
  using MappedType = typename Base::mapped_type;
  using ValueType = typename Base::value_type;

  /// Inherit constructors
  using Base::Base;

  /**
   * Access element at index with bounds checking if USE_STD_CONTAINER_BOUNDS_CHECKING is true.
   * Otherwise, uses operator[] for unchecked access.
   * @param index Index of the element to access.
   * @return Const reference to the element at the specified index.
   * @throws std::out_of_range if index is out of bounds.
   */
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

  /**
   * Access element at index with bounds checking if USE_STD_CONTAINER_BOUNDS_CHECKING is true.
   * Otherwise, uses operator[] for unchecked access.
   * @param index Index of the element to access.
   * @return Const reference to the element at the specified index.
   * @throws std::out_of_range if index is out of bounds.
   */
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

  /**
   * @brief If a key equivalent to k already exists in the container
   * Otherwise, inserts a new element into the container with key k
   * @param k The key used both to look up and to insert if not found
   * @return A reference to the mapped value
   */
  auto & get_inserted( KeyType && k )
  {
    return this->try_emplace( k ).first->second;
  }

  /**
   * @brief If a key equivalent to k already exists in the container
   * Otherwise, inserts a new element into the container with key k
   * @param k The key used both to look up and to insert if not found
   * @return A reference to the mapped value
   */
  auto & get_inserted( KeyType const & k )
  {
    return this->try_emplace( k ).first->second;
  }

};

} //namespace internal

/**
 * type alias for stdMap
 * @tparam Key The unique stdMap key.
 * @tparam T Type of elements in the stdMap.
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

/// @endcond

/**
 * @name Ordered and unordered map types.
 */
///@{

/**
 * @brief Base template for ordered and unordered maps.
 * @tparam TKEY key type
 * @tparam TVAL value type
 * @tparam SORTED a @p std::integral_constant<bool> indicating whether map is ordered
 */
template< typename TKEY, typename TVAL, typename SORTED >
class mapType
{};

/// @cond DO_NOT_DOCUMENT
template< typename TKEY, typename TVAL >
class mapType< TKEY, TVAL, std::integral_constant< bool, true > > : public stdMap< TKEY, TVAL >
{
  using stdMap< TKEY, TVAL >::stdMap; // enable list initialization
};

template< typename TKEY, typename TVAL >
class mapType< TKEY, TVAL, std::integral_constant< bool, false > > : public stdUnorderedMap< TKEY, TVAL >
{
  using stdUnorderedMap< TKEY, TVAL >::stdUnorderedMap; // enable list initialization
};
/// @endcond


/**
 * @name Ordered and unordered map types.
 *  OLD VERSION, TO BE REMOVED
 */
///@{

/**
 * @brief Base template for ordered and unordered maps.
 * @tparam TKEY key type
 * @tparam TVAL value type
 * @tparam SORTED a @p std::integral_constant<bool> indicating whether map is ordered
 */
template< typename TKEY, typename TVAL, typename SORTED >
class mapBase
{};

/// @cond DO_NOT_DOCUMENT
template< typename TKEY, typename TVAL >
class mapBase< TKEY, TVAL, std::integral_constant< bool, true > > : public std::map< TKEY, TVAL >
{
  using std::map< TKEY, TVAL >::map; // enable list initialization
};

template< typename TKEY, typename TVAL >
class mapBase< TKEY, TVAL, std::integral_constant< bool, false > > : public std::unordered_map< TKEY, TVAL >
{
  using std::unordered_map< TKEY, TVAL >::unordered_map; // enable list initialization
};
/// @endcond

namespace internal
{

/**
 * Wrapper for the underlying std::aray that allows toggling between bounds-checked access
 * (using at()) and unchecked access (using operator[]).
 * @tparam T Type of elements in stdArray
 * @tparam N The number of fixed element in the array
 * @tparam USE_STD_CONTAINER_BOUNDS_CHECKING A boolean flag to enable or disable bounds checking.
 */
template< class T,
          std::size_t N,
          bool USE_BOUNDS_CHECKING = false >
struct StdArrayWrapper : public std::array< T, N >
{
public:
  /// Type alias for the base class (i.e., stdArray)
  using Base = std::array< T, N >;

  /**
   * Access element at index with bounds checking if USE_STD_CONTAINER_BOUNDS_CHECKING is true.
   * Otherwise, uses operator[] for unchecked access.
   * @param index Index of the element to access.
   * @return Const reference to the element at the specified index.
   * @throws std::out_of_range if index is out of bounds.
   */
  constexpr T & operator[]( size_t const index )
  {
    if constexpr (USE_BOUNDS_CHECKING)
    {
      return this->at( index );
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
   * @return Const reference to the element at the specified index.
   * @throws std::out_of_range if index is out of bounds.
   */
  constexpr T const & operator[]( size_t const index ) const
  {
    if constexpr (USE_BOUNDS_CHECKING)
    {
      return this->at( index );
    }
    else
    {
      return Base::operator[]( index );
    }
  }

};
} //namespace internal

/**
 * @tparam T Type of elements in the array.
 * @tparam N The number of fixed element in the array.
 * @note we use a struct rather than an alias for taking advantage of deduction guide
 */
template< class T, std::size_t N >
struct stdArray : public internal::StdArrayWrapper< T, N, true >
{};

/**
 * @brief Deduction guide for stdArray
 * Allows the element type and array size to be automatically deduced from the initialization list.
 * @code
 * stdArray a{1, 2, 3}; // deduces stdArray<int, 3>
 * @endcode
 * @tparam _Tp Type of the first element provided.
 * @tparam _Up Types of the other elements provided.
 */
template< typename _Tp, typename ... _Up >
stdArray( _Tp, _Up ... )
  ->stdArray< std::enable_if_t< (std::is_same_v< _Tp, _Up > && ...), _Tp >,
              1 + sizeof...(_Up) >;

/**
 * @brief Helper function that convert a std::array into a stdArray by expanding its elements
 * @tparam T The type of the elements
 * @tparam N The number of fixed element in the array
 * @tparam Is An integer parameter pack representing the indices
 * @param arr The input std::array to be converted.
 * @return A constexpr stdArray< T, N >
 */
template< typename T, std::size_t N, std::size_t... Is >
constexpr stdArray< T, N > to_stdArray_impl( std::array< T, N > const & arr, std::index_sequence< Is... > )
{
  return {{arr[Is] ...}};
}

/**
 * @brief Convert an std::array to an stdArray.
 * @tparam T Type of elements in stdArray
 * @tparam N The number of fixed element in the array
 * @param arr The std::array to convert
 * @return A stdArray< T, N >
   @note we don't want implicitly convert an std::array into a stdArray.
 */
template< typename T, std::size_t N >
constexpr stdArray< T, N > to_stdArray( std::array< T, N > const & arr )
{
  return to_stdArray_impl( arr, std::make_index_sequence< N >{} );
}


} // namespace geos

/**
 * @namespace std
 * @brief Partial specialization for stdArray.
 * Those specialization serve to use geos::stdArray as tuples.
 */
namespace std
{

/**
 * @brief Provides access to the number of elements in a tuple as a compile-time constant expression.
 * @tparam T Type of elements in stdArray
 * @tparam N The number of fixed element in the array
 */
template< typename T, size_t N >
struct tuple_size< geos::stdArray< T, N > >
  : public integral_constant< size_t, N > { };

/**
 * @brief Provides compile-time indexed access to the types of the elements of the tuple.
 * @tparam Ind The index element of the tuple to access
 * @tparam Tp Type of elements in stdArray
 * @tparam Nm The number of fixed element in the array
 */
template< size_t Ind, typename Tp, size_t Nm >
struct tuple_element< Ind, geos::stdArray< Tp, Nm > >
{
  static_assert( Ind < Nm, "array index is in range" );
  /// The type of the Ind-th element, which is Tp for a homogeneous array.
  using type = Tp;
};

#if __cplusplus >= 201703L
/**
 * @brief Helper variable template
 * @tparam Tp Type of elements in stdArray
 * @tparam Nm The number of fixed element in the array
 */
template< typename Tp, size_t Nm >
inline constexpr size_t tuple_size_v< geos::stdArray< Tp, Nm > > = Nm;

template< typename Tp, size_t Nm >
inline constexpr size_t tuple_size_v< const geos::stdArray< Tp, Nm > > = Nm;
#endif

template< typename Tp, size_t Nm >
struct __is_tuple_like_impl< geos::stdArray< Tp, Nm > > : true_type
{ };
}

#endif /* GEOS_COMMON_STD_CONTAINER_WRAPPERS_HPP */
