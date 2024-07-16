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
 * @file STLUtils.hpp
 */

#ifndef GEOS_DATAREPOSITORY_STLUTILS_HPP_
#define GEOS_DATAREPOSITORY_STLUTILS_HPP_

#include <tuple>

namespace geos
{

// ItemIterator to iterate over the N-th element of iterators of tuple-like structures
template < typename Iterator, std::size_t N >
class ItemIterator
{
public:
  using iterator_category = std::forward_iterator_tag;
  using base_value_type = typename std::iterator_traits<Iterator>::value_type;
  using value_type = std::remove_reference_t<decltype(std::get<N>(std::declval<base_value_type>()))>;
  using difference_type = std::ptrdiff_t;
  using reference = decltype(std::get<N>(*std::declval<Iterator>()));
  using pointer = value_type*;

  explicit ItemIterator( Iterator it ) : m_iter(it) {}

  reference operator*() const { return std::get<N>(*m_iter); }
  pointer operator->() const { return &std::get<N>(*m_iter); }

  // Prefix increment
  ItemIterator & operator++()
  {
    ++m_iter;
    return *this;
  }

  // Postfix increment
  ItemIterator operator++(int)
  {
    ItemIterator tmp = *this;
    ++(*this);
    return tmp;
  }

  bool operator==(const ItemIterator& other) const
  {
    return m_iter == other.m_iter;
  }

  bool operator!=(const ItemIterator& other) const
  {
    return m_iter != other.m_iter;
  }

private:
  Iterator m_iter;
};

// ItemView to provide a view over a specific element of tuple-like elements in containers
template < std::size_t N, class Container >
class ItemView
{
public:
  using container_type = Container;
  using iterator = ItemIterator< typename container_type::iterator, N >;

  ItemView( container_type & container ) : m_container( container ) {}

  iterator begin() { return iterator( m_container.begin() ); }
  iterator end() { return iterator( m_container.end() ); }

private:
  container_type & m_container;
};

}

#endif