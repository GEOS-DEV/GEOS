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
 * @file TableLayout.hpp
 */

#ifndef GEOS_COMMON_TABLELAYOUT_HPP
#define GEOS_COMMON_TABLELAYOUT_HPP

#include "common/DataTypes.hpp"

namespace geos
{

class TableLayout
{

public:
  enum Alignment { right, left, middle };

  enum MarginValue : integer
  {
    tiny = 0,
    small = 1,
    medium = 2,
    large = 3
  };

  /**
   * @brief Structure to set up each colum parameters.
   */
  struct ColumnParam
  {
    // A vector containing all string for a header name
    std::vector< string > headerName;
    // Alignment for a column. By default aligned to the right side
    Alignment alignment = Alignment::right;
    // A boolean to display a colummn
    bool enabled = true;
  };

  /**
   * @brief Construct a new Table object, all values in the table are centered by default
   * @param columnNames
   */
  TableLayout( std::vector< string > const & columnNames );

  /**
   * @brief Construct a new Table object by specifying value alignment and optionally their displays based to log levels
   * level
   * @param columnParameter List of structures to set up each colum parameters.
   */
  TableLayout( std::vector< ColumnParam > const & columnParameter );

  void setMargin( MarginValue marginType );

  enum Section { header, values };

  /**
   * @brief Struct for a column.
   */
  struct Column
  {
    ColumnParam parameter;
    // A vector containing all column values
    std::vector< string > columnValues;
    // The largest string in the column
    string m_maxStringSize;
  };

  std::vector< Column > m_columns;

private:
  integer borderMargin;
  integer columnMargin;

};
}

#endif
