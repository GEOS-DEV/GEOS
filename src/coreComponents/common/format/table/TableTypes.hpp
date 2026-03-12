

/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TableData.hpp
 */


#ifndef GEOS_COMMON_FORMAT_TABLETYPES_HPP
#define GEOS_COMMON_FORMAT_TABLETYPES_HPP

#include "common/DataTypes.hpp"

namespace geos
{

/**
 * @brief The different type a cell can handle
 */
enum class CellType : integer
{
  Header,
  Value,
  Separator,
  MergeNext,
  Hidden
};

/**
 * @brief Class for retrieving errors in the table classes
 */
class TableErrorListing
{

public:
  /**
   * @brief Add an error that will be display at the end of the table
   * @param text The string error to display.
   */
  void addError( string_view text );

  /**
   * @return true if an error has already been added
   */
  bool hasErrors() const;

  /**
   * @brief Clear all error when calling toString()
   */
  void clear();

  /// The iterator alias for the  errors vector of string
  using Iterator = std::vector< string >::const_iterator;

  /**
   * @return An Iterator pointing to the first element of the errors vector
   */
  Iterator begin() const
  { return m_errorList.begin(); }

  /**
   * @return An Iterator pointing to the last element of the errors vector
   */
  Iterator end() const
  { return m_errorList.end(); }

  /**
   * @brief Append a vector of string to the errors vector.
   * @param errors A vector of string to append
   */
  void appendErrors( std::vector< string > const & errors )
  { m_errorList.insert( m_errorList.end(), errors.begin(), errors.end() );}

  /**
   * @return A const reference to the errors vector.
   */
  std::vector< string > const & getErrors() const
  { return m_errorList; }

  /**
   * @return A reference to the errors vector.
   */
  std::vector< string > & getErrors()
  { return m_errorList; }

private:
  /// Contain all the errors  to display at the end of the table
  std::vector< string > m_errorList;
};

inline void TableErrorListing::addError( string_view text )
{ m_errorList.emplace_back( text ); }

inline bool TableErrorListing::hasErrors() const
{ return !m_errorList.empty(); }

inline void TableErrorListing::clear()
{ m_errorList.clear(); }

}

#endif /* GEOS_COMMON_FORMAT_TABLETYPES_HPP */
