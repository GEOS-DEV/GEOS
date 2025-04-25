

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
  { return errorText.begin(); }

  /**
   * @return An Iterator pointing to the last element of the errors vector
   */
  Iterator end() const
  { return errorText.end(); }

  /**
   * @brief Append a vector of string to the errors vector.
   * @param errors A vector of string to append
   */
  void appendErrors( std::vector< string > & errors )
  { errorText.insert( errorText.end(), errors.begin(), errors.end() );}

  /**
   * @return A const reference to the errors vector.
   */
  std::vector< string > const & getErrors() const
  { return errorText; }

private:
  /// Contain all the errors  to display at the end of the table
  std::vector< string > errorText;
};

/**
 * @brief Add an error to the error vector
 * @param text The string view of the error
 */
inline void TableErrorListing::addError( string_view text )
{ errorText.emplace_back( text ); }

/**
 * @return Return true if the vector contain any errors, false otherwise
 */
inline bool TableErrorListing::hasErrors() const
{ return errorText.size() != 0; }

/**
 * @brief Clear the errors vector
 */
inline void TableErrorListing::clear()
{ errorText.clear(); }

}

#endif /* GEOS_COMMON_FORMAT_TABLETYPES_HPP */
