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
 * @file Section.hpp
 */

#ifndef GEOS_COMMON_SECTION_HPP
#define GEOS_COMMON_SECTION_HPP

#include "common/DataTypes.hpp"

namespace geos
{

class Section
{
public:

  Section( string_view m_sectionTitle );

  /**
   * @brief Add a description to the section by concatening a description name and descriptions values.
   * @param descriptionName The description name
   * @param args Descriptions values to be aligned.
   * Descriptions values can be be any types and will be aligned
   */
  template< typename ... Args >
  void addDescription( string_view descriptionName, Args const & ... args );

  /**
   * @brief Add a description to the section
   * @param description The string value of the description
   */
  void addDescription( string_view description );

  /**
   * @brief Set the minimal width of a row
   * @param minWidth The minimal width of the table
   */
  void setMinWidth( integer const & minWidth );

  /**
   * @brief Draw the first part of the section. It include the title and optionnaly, the description(s);
   * @param os An output stream (by default, std::cout)
   */
  void beginSection( std::ostream & os = std::cout );

  /**
   * @brief Draw the last part of the section. It include the title
   * @param oss An output stream (by default, std::cout)
   */
  void endSection( std::ostream & oss = std::cout ) const;

private:

  /**
   * @brief Compute the max string size (m_rowLength) between title and the description(s)
   * @param m_sectionTitle The table title
   * @param descriptions The descriptions vector
   * @return The max row length of the section
   */
  void computeMaxWidth( string_view m_sectionTitle,
                        std::vector< string > const & descriptions );

  /**
   * @brief Build a description from the name and description values
   * @param descriptionName The decription name
   * @param decriptionsValues The description values
   */
  void formatAndInsertDescriptions( string_view descriptionName,
                                    std::vector< string > const & decriptionsValues );

  /// Vector containing all descriptions
  std::vector< string > m_descriptions;

  /// title of section
  string m_sectionTitle;
  /// section length
  integer m_rowLength;
  /// min width of section length
  integer m_rowMinWidth;

  /// description border margin
  static constexpr integer m_marginBorder = 2;
  /// numbers of character used as border
  static constexpr integer m_nbSpecialChar = 2;

};

template< typename ... Args >
void Section::addDescription( string_view descriptionName, Args const &... args )
{
  std::vector< string > descriptionsValues;
  ( [&] {
    static_assert( has_formatter_v< decltype(args) >, "Argument passed in addRow cannot be converted to string" );
    string const value = GEOS_FMT( "{}", args );
    descriptionsValues.push_back( value );
  } (), ...);

  formatAndInsertDescriptions( descriptionName, descriptionsValues );
}
}

#endif
