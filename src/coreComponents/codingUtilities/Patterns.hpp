/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file Patterns.hpp
 */

#ifndef GEOS_CODINGUTILIIES_PATTERNS_HPP_
#define GEOS_CODINGUTILIIES_PATTERNS_HPP_

#include "traits.hpp"
#include "STLUtils.hpp"
#include <functional>
#include <memory>
#include <type_traits>
#include <vector>

namespace geos
{

/**
 * @brief Interface for observer objects.
 *
 * This template class defines the interface for observers that need to be notified of changes
 * by an Observable object of type T.
 *
 * @tparam T The type of data observed.
 */
template < typename T >
class Observer
{
public:
    /**
     * @brief Virtual destructor for the observer interface.
     */
    virtual ~Observer() {}

    /**
     * @brief Update method called by Observable when changes occur.
     *
     * This method is called to update the observer with the new state of the observable object.
     *
     * @param data The updated data from the observable.
     */
    virtual void update( const T & data ) = 0;
};

/**
 * @brief Class for objects that can be observed.
 *
 * This template class manages a list of observers and notifies them of changes.
 * Observers must implement the Observer interface.
 *
 * @tparam T The type of data to be observed.
 */
template < typename T >
class Observable
{
public:
  /**
   * @brief Destructor which clears observers.
   */
  virtual ~Observable()
  {
    m_observers.clear();
  }

  /**
   * @brief Adds an observer to the list.
   *
   * @param observer Pointer to the observer to add.
   */
  void observe( std::shared_ptr< Observer< T > > observer )
  {
    m_observers.push_back(observer);
  }

  /**
   * @brief Removes an observer from the list.
   *
   * @param observer Pointer to the observer to remove.
   */
  void forget( std::shared_ptr< Observer< T > > observer )
  {
    m_observers.erase(std::remove(m_observers.begin(), m_observers.end(), observer), m_observers.end());
  }

  /**
   * @brief Notifies all observers about a change.
   *
   * @param data The data to be passed to the observers' update method.
   */
  void notify( const T & data )
  {
    for ( auto & observer : m_observers )
    {
      observer->update(data);
    }
  }
private:
  std::vector< std::shared_ptr< Observer< T > > > m_observers;  ///< List of observers.
};

/**
 * @brief Template class for observing changes in objects of type T.
 *
 * This class templates an observer mechanism that allows for monitoring changes
 * on an object of type T. Changes are detected based on a predicate function
 * which should return true when a change of interest occurs. Upon detection,
 * a callback function is executed.
 *
 * Note: This is not a traditional observer as it requires actively monitoring
 *       the underlying object state. It is useful when the observed type
 *       cannot be made to inherit from the tradition Observable interface.
 *       But in other circumstances is less efficient.
 *
 * @tparam T The type of the object to observe.
 */
template< typename T, typename R >
class ActiveObserver
{
public:
  /**
   * @brief Constructs a new ActiveObserver object.
   *
   * @param observed The object to observe.
   * @param predicate A function that checks the observed object for changes. Should return true if a change that requires a callback execution is detected.
   * @param notification A function to be called when the predicate returns true. It is used to handle or respond to changes.
   */
  ActiveObserver(const T & observed, std::function<bool(const T&)> predicate, std::function<R(T&)> notification)
    : m_observed(observed),
      m_predicate(predicate),
      m_notification(notification)
  {}

  /**
   * @brief Checks the observed object for changes.
   *
   * Evaluates the predicate on the observed object. If the predicate returns true,
   * the notification function is executed. This function should be called periodically
   * or in scenarios where checking for changes is required.
   */
  R check()
  {
    if( m_predicate( m_observed ) )
    {
      return m_notification( m_observed );
    }
    return {};
  }

private:
  T m_observed;  ///< The observed object.
  std::function< bool ( const T & ) > m_predicate;  ///< Predicate function to detect changes.
  std::function< R ( T & ) > m_notification;  ///< Callback function to execute on change.
};

// Base template: assume T is not a pointer by default
template< typename T >
struct is_any_pointer : std::false_type {};

// Specialization for raw pointers
template< typename T >
struct is_any_pointer< T * > : std::true_type {};

// Specialization for std::shared_ptr
template< typename T >
struct is_any_pointer< std::shared_ptr< T > > : std::true_type {};

// Specialization for std::unique_ptr
template< typename T >
struct is_any_pointer< std::unique_ptr< T > > : std::true_type {};

template < typename T >
constexpr bool is_any_pointer_v = is_any_pointer< T >::value;

// This trait ensures that we only change the type if the node can actually be dereferenced.
template< typename T, bool = is_any_pointer_v< std::remove_cv_t< std::remove_reference_t< T > > > >
struct nested_dereference
{
  using type = std::remove_cv_t< std::remove_reference_t< T > >;  // If T is not dereferenceable, return T itself
};

// Specialize for dereferenceable types
template< typename T >
struct nested_dereference< T, true >
{
  using type = typename nested_dereference< std::remove_cv_t< std::remove_reference_t< decltype( *std::declval< T & >( ) ) > > >::type;
};

template< typename T >
using nested_dereference_t = typename nested_dereference< T >::type;

template < typename NODE >
decltype(auto) inline dereference( NODE & node )
{
  if constexpr ( is_any_pointer_v< NODE > )
  {
    if constexpr ( is_any_pointer_v< decltype( *node ) > )
    {
      return dereference( *node );
    }
    else
    {
      return *node;
    }
  }
  else
  {
    return node;
  }
}

/**
 * @brief A functor to get children using a member function 'children()'
 *        Assumes children() returns a container or iterable.
 */
struct DefaultChildAccessor
{
  template < typename Node >
  auto operator()( Node & node ) const -> decltype(auto)
  {
    return node.children();
  }
};

/**
 * @brief StaticTreeIteration is a structure for traversal and processing of tree structures which
 *        WILL NOT BE MODIFIED DURING traversal and processing of the nodes.
 *
 * @tparam NodeType The type of nodes in the tree.
 */
template < typename ChildAccessor = DefaultChildAccessor >
struct StaticTreeIteration
{
  template < typename NODE, typename LAMBDA >
  static void processTree( NODE & node, LAMBDA && lambda )
  {
    static_checks< NODE, LAMBDA >( );
    lambda( dereference( node ) );
    for ( auto & child : getChildren( node ) )
    {
      processTree( dereference( child ), std::forward< LAMBDA >( lambda ) );
    }
  }
private:
  template < typename NODE >
  static auto getChildren( NODE & node ) -> decltype(auto)
  {
    static ChildAccessor accessor;
    return accessor( dereference( node ) );
  }
  template < typename NODE, typename LAMBDA >
  static constexpr void static_checks()
  {
    // TODO (wrt/c++20) replace with concepts
    using DecayedNode = nested_dereference_t< NODE >;
    // Ensure that the lambda can be invoked with a reference to the dereferenced node type
    static_assert( std::is_invocable_v< LAMBDA, DecayedNode & >,
                    "The lambda must be invokable with a reference to the dereferenced node type" );
    static_assert( std::is_invocable_v< ChildAccessor, DecayedNode & >,
                    "ChildAccessor must be invokable on the dereferenced node type");

    using ChildrenType = decltype( std::declval< ChildAccessor >()( std::declval< DecayedNode & >() ) );
    using IteratorType = decltype( std::begin( std::declval< ChildrenType & >() ) );
    using DereferencedIteratorType = decltype( *std::declval< IteratorType >() );
    using ElementOfChildrenType = nested_dereference_t< DereferencedIteratorType >;
    static_assert( std::is_same_v< std::decay_t< ElementOfChildrenType >, DecayedNode >,
                    "The elements of the child accessor's result must be dereferencable to the decayed node type");
  }
};

/**
 * @brief DynamicTreeIteration is a structure for traversal and processing of tree structures which may
 *        be modified ONLY BY ADDITION OF NODES during processing.
 */
template < typename Identifier, typename ChildAccessor = DefaultChildAccessor >
struct DynamicTreeIteration
{
public:
  /**
   * @brief Processes the tree dynamically and applies the provided lambda function to each node.
   *
   * This static function repeatedly processes all nodes in the tree until no new nodes are detected.
   *
   * @param root The root node of the tree to start processing.
   * @param lambda A lambda function to process each node.
   * @param identifier A function that extracts a unique identifier from a node.
   */
  template< typename NODE, typename LAMBDA >
  static void processTree( NODE & root, LAMBDA && lambda )
  {
    using T = decltype( getIdentifier( root ) );
    static_checks< NODE, LAMBDA, T >( );
    std::set< T > visited;
    bool foundNew = true;
    std::function < void ( NODE & ) > processSubtree;
    processSubtree = [&processSubtree,&lambda,&visited,&foundNew] ( NODE & node )
    {
      decltype(auto) id = getIdentifier( node );
      if ( visited.find( id ) == visited.end() )
      {
        lambda( dereference( node ) );
        visited.insert( id );
        foundNew = true;
      }
      for ( auto & child : getChildren( node ) )
      {
        processSubtree( dereference( child ) );
      }
    };

    while( foundNew )
    {
      foundNew = false;
      processSubtree( root );
    }
  }

private:
  template < typename NODE >
  static auto getIdentifier( NODE & node ) -> decltype(auto)
  {
    static Identifier identifier;
    return identifier( dereference( node ) );
  }

  template < typename NODE >
  static auto getChildren( NODE & node ) -> decltype(auto)
  {
    static ChildAccessor accessor;
    return accessor( dereference( node ) );
  }
  template < typename NODE, typename LAMBDA , typename T >
  static constexpr void static_checks( )
  {
    // TODO (wrt/c++20) replace with concepts
    using DecayedNode = nested_dereference_t< NODE >;
    // Validate that the lambda can be invoked with a node
    static_assert( std::is_invocable_v< LAMBDA, DecayedNode & >,
                   "Lambda must be invokable with a reference to the node type" );
    static_assert(std::is_invocable_r_v< T, Identifier, DecayedNode & >,
                  "Identifier must be invokable with a reference to the node type and must return a value");
    // Check if the type T can be compared using less-than operator
    static_assert( std::is_invocable_r_v<bool, std::less<>, T, T >,
                   "The return type of the identity extractor must be comparable using the less-than operator" );

    using ChildrenType = decltype( std::declval< ChildAccessor >()( std::declval< DecayedNode & >() ) );
    using IteratorType = decltype( std::begin( std::declval< ChildrenType & >() ) );
    using DereferencedIteratorType = decltype( *std::declval< IteratorType >() );
    using ElementOfChildrenType = nested_dereference_t< DereferencedIteratorType >;
    static_assert( std::is_same_v< std::decay_t< ElementOfChildrenType >, DecayedNode >,
                    "The elements of the child accessor's result must be dereferencable to the decayed node type");
  }
};

// Helper template to create a function type from a tuple
template < typename Tuple >
struct tuple_to_function;

template <typename... Args>
struct tuple_to_function<std::tuple<Args...>>
{
  using type = std::function<void(Args...)>;
};

template <typename FlagType, typename ArgsTuple>
class FlagCallbacks
{
public:
  using flag_type = FlagType;
  using arg_tuple_type = ArgsTuple;
  using Callback = typename tuple_to_function<ArgsTuple>::type;

  void registerCallback(FlagType flag, Callback callback)
  {
    m_callbacks[flag] = std::move(callback);
  }

  template <typename... Args>
  void executeCallbacks(const std::set<FlagType> & flags, Args... args)
  {
    for (const auto & flag : flags)
    {
      auto it = m_callbacks.find(flag);
      if (it != m_callbacks.end())
      {
        it->second( args... );
      }
    }
  }

private:
  std::unordered_map<FlagType, Callback> m_callbacks;
};

template < typename FlagType, typename... ArgsTuples >
class FlagCallbackCoordinator
{
public:
  using flag_type = FlagType;

  template < typename... Args >
  void registerCallback( FlagType flag, std::function< void(Args...) > callback)
  {
    static_assert(
      std::disjunction< std::is_same< std::tuple< Args... >, ArgsTuples >...>::value,
      "The types provided to registerCallback do not match any tuple in m_keyedCallbacks."
    );

    auto & callbacks = std::get< FlagCallbacks< FlagType, std::tuple< Args... > > >( m_keyedCallbacks );
    callbacks.registerCallback( flag, callback );
  }

  template < typename... Args >
  void executeCallbacks( const std::set< FlagType > & flags, Args... args )
  {
    static_assert(
      std::disjunction< std::is_same< std::tuple< Args... >, ArgsTuples >... >::value,
      "The types provided to executeCallbacks do not match any tuple in m_keyedCallbacks."
    );

    auto & callbacks = std::get< FlagCallbacks< FlagType, std::tuple< Args... > > >( m_keyedCallbacks );
    callbacks.executeCallbacks( flags, args... );
  }

private:
  std::tuple< FlagCallbacks< FlagType, ArgsTuples >... > m_keyedCallbacks;
};
} // namespace geos

#endif