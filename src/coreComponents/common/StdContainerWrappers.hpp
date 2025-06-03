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
          bool USE_BOUNDS_CHECKING = false
          >
class StdVectorWrapper : public std::vector< T, Allocator >
{
public:
  typedef typename std::_Vector_base< T, Allocator >::_Tp_alloc_type _Tp_alloc_type;
  typedef __gnu_cxx::__alloc_traits< _Tp_alloc_type > _Alloc_traits;
  /// Type alias for the base class (i.e., std::vector)
  using Base = std::vector< T, Allocator >;
  // Inherit constructors
  StdVectorWrapper(): std::vector< T, Allocator >()
  {}

  StdVectorWrapper( const Allocator & __a ): std::vector< T, Allocator >( __a )
  {}

  StdVectorWrapper( size_t __n, const Allocator & __a = Allocator()): std::vector< T, Allocator >( __n, __a )
  {}

  StdVectorWrapper( size_t __n, const T & __value,
                    const Allocator & __a = Allocator())
    : std::vector< T, Allocator >( __n, __value, __a )
  {}

  // StdVectorWrapper
  StdVectorWrapper( const StdVectorWrapper & __x )
    : std::vector< T, Allocator >( __x.size(),
                                   _Alloc_traits::_S_select_on_copy( __x._M_get_Tp_allocator()))
  {}

  StdVectorWrapper( const StdVectorWrapper & __x, const Allocator & __a )
    : std::vector< T, Allocator >( __x.size(), __a )
  {}

  // INITALIZER
  StdVectorWrapper( std::initializer_list< T > __l, const Allocator & __a =  Allocator())
    : std::vector< T, Allocator >( __l, __a )
  {}

  // MOVE
  StdVectorWrapper( StdVectorWrapper && __x )
    : std::vector< T, Allocator >( std::move( __x ))
  {}

  StdVectorWrapper( const StdVectorWrapper && __rv, const Allocator & __m )
    : std::vector< T, Allocator >( std::move( __rv ), __m )
  {}

  //Iterator
  template< typename _InputIterator >
  StdVectorWrapper( _InputIterator __first, _InputIterator __last,
                    const Allocator & __a = Allocator())
    : std::vector< T, Allocator >( __first, __last, __a )
  {}

  // operator
  StdVectorWrapper & operator=( const StdVectorWrapper & __x )
  {
    if( this != &__x )
    {
      std::vector< T, Allocator >::operator=( __x );  // Use base class assignment
    }
    return *this;
  }

  StdVectorWrapper & operator=( StdVectorWrapper && __x ) noexcept
  {
    if( this != &__x )
    {
      std::vector< T, Allocator >::operator=( std::move(__x));   // Move assignment
    }
    return *this;
  }

  StdVectorWrapper & operator=( std::initializer_list< T > __l )
  {
    std::vector< T, Allocator >::operator=( __l );  // Use base class assignment
    return *this;
  }

  // Copy constructor
  StdVectorWrapper( const std::vector< T, Allocator > & vec )
    : std::vector< T, Allocator >( vec ) {}

  // Move constructor
  StdVectorWrapper( std::vector< T, Allocator > && vec )
    : std::vector< T, Allocator >( std::move( vec )) {}

  // Conversion constructor
  template< typename U, typename A >
  StdVectorWrapper( std::vector< U, Allocator > & vec )
    : std::vector< T, Allocator >( vec.begin(), vec.end()) {}


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
 * @tparam Allocator Allocator type for the vector.
 */
template< typename T, typename Allocator = std::allocator< T > >
using stdVector = internal::StdVectorWrapper< T, Allocator, USE_STD_CONTAINER_BOUNDS_CHECKING >;



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
