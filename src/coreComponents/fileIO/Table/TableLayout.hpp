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
  enum Alignment { right, left, center };

  enum MarginValue : integer
  {
    tiny = 0,
    small = 1,
    medium = 2,
    large = 3
  };

  /**
   * @brief Enumeration for table sections.
   */
  enum Section { header, values };

  /**
   * @brief Structure to set up each colum parameters.
   */
  struct ColumnParam
  {
    string columnName;
    // Alignment for a column. By default aligned to the right side
    Alignment alignment = Alignment::right;
    // A boolean to display a colummn
    bool enabled = true;
    // Vector containing substring column name delimited by "\n"
    std::vector< string > splitColumnName;

    /**
     * @brief Construct a ColumnParam object with the specified name and alignment.
     * @param name The name of the column
     * @param align The alignment of the column
     */
    ColumnParam( std::string const & name, Alignment align )
      : columnName( name ), alignment( align )
    {}

    /**
     * @brief Construct a ColumnParam object with the specified name, alignment, and display flag.
     * @param name The name of the column
     * @param align The alignment of the column
     * @param display Flag indicating whether the column is enabled
     */
    ColumnParam( std::string const & name, Alignment align, bool display )
      : columnName( name ), alignment( align ), enabled( display )
    {}
  };

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

  /**
   * @brief Construct a new Table object, all values in the table are centered by default
   * @param columnNames The names of the columns
   * @param title The table name
   */
  TableLayout( std::vector< string > const & columnNames, string const & title = "" );

  /**
   * @brief Construct a new Table object by specifying value alignment and optionally their displays based to log levels
   * level
   * @param columnParameter List of structures to set up each colum parameters.
   * @param title The table name
   */
  TableLayout( std::vector< ColumnParam > const & columnParameter, string const & title = "" );

  /**
   * @return The columns vector
   */
  std::vector< Column > const & getColumns() const;

  /**
   * @return The table name
   */
  string_view getTitle() const;

  /**
   * @return The border margin
   */
  integer const & getBorderMargin() const;

  /**
   * @return The column margin
   */
  integer const & getColumnMargin() const;

  /**
   * @return The margin title
   */
  integer const & getMarginTitle() const;

private:

  /**
   * @brief Set the minimal margin width between row content and borders.
   * @param marginType The margin value
   */
  void setMargin( MarginValue marginValue );

  std::vector< Column > m_columns;
  string m_tableTitle;
  integer m_borderMargin;
  integer m_columnMargin;
  integer m_marginTitle = 2;

};
}

#endif
