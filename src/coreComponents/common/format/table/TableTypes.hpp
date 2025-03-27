

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

class ErrorListing
{

public:
  /// View on cell error grouping the cell information to display at the end of the table
  //CellLayoutRows errors;   //TODO à enlever crée au last moment
  /// Contain all the errors
  std::vector< string > errorText;

  /**
   * @brief Add an error that will be display at the end of the table
   * @param text The string error to display.
   * @param nbCells The numbers cells that must be equal to the number of a CellLayoutRow
   * present in headerCellsLayout or dataCellsLayout
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
};


inline void ErrorListing::addError( string_view text )
{ errorText.emplace_back( text ); }

inline bool ErrorListing::hasErrors() const
{ return errorText.size() != 0; }

inline void ErrorListing::clear()
{ errorText.clear(); }

}

#endif /* GEOS_COMMON_FORMAT_TABLETYPES_HPP */
