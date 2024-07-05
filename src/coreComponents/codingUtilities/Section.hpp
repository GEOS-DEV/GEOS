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

  Section();

  /**
   * @brief Set the name of the section
   * @param m_sectionTitle The name of the section
   */
  void setName( string_view m_sectionTitle );

  /**
   * @brief Add a description to the section and composed by a description name and variadic values.
   * Use to align variadic values of the same description
   * @param descriptionName The description name 
   * @param args Values to be aligned.
   */
  template< typename ... Args >
  void addDescription( string const & descriptionName, Args const & ... args );

  /**
   * @brief Add a description to the section
   * @param description The string value of the description
   */
  void addDescription( string const & description );

  /**
   * @brief Set the minimal width of a row
   * @param minWidth The minimal width of the table
   */
  void setMinWidth( integer const & minWidth );

  /**
   * @brief Draw the first part of the section. It include the title and optionnaly, the description(s);
   * @param os An output stream (by default, std::cout)
   */
  void begin( std::ostream & os = std::cout );

  /**
   * @brief Draw the last part of the section. It include the title
   * @param oss An output stream (by default, std::cout)
   */
  void end( std::ostream & oss = std::cout );

private:

  /**
   * @brief Compute the max string size (m_rowLength) between title and the description(s)
   * @param m_sectionTitle The table title
   * @param descriptions The descriptions vector
   */
  void computeMaxRowSize( string const & m_sectionTitle,
                          std::vector< string > const & descriptions );
  /**
   * @brief Build a description from the name and variadic values descriptions 
   */
  void buildAlignDescription();

  /**
   * @brief Cleans all buffers used in the construction of a section
   */
  void clear();

  /// Vector containing all description
  std::vector< string > m_descriptions;
  /// Used if the variadic addDescription has been called
  /// Containing all "key" description name
  std::vector< string > m_descriptionNames;
  /// Used if the variadic addDescription has been called
  /// Containing all description values
  std::vector< std::vector< string > > m_descriptionsValues;

  /// title of section
  string m_sectionTitle;
  /// section length
  integer m_rowLength;
  /// min width of section length
  integer m_rowMinWidth;

  /// description border margin
  static constexpr integer m_marginBorder = 2;
  /// character used as border
  static constexpr integer m_nbSpecialChar = 2;
  /// (Temporary ?) special char with key name. I.E =>- "name": => 3char
  static constexpr integer m_embeddingName = 4;

};

template< typename ... Args >
void Section::addDescription( string const & descriptionName, Args const &... args )
{
  std::vector< string > descriptions;
  ( [&] {
    static_assert( has_formatter_v< decltype(args) >, "Argument passed in addRow cannot be converted to string" );
    string const value = GEOS_FMT( "{}", args );
    descriptions.push_back( value );
  } (), ...);

  m_descriptionsValues.push_back( descriptions );
  m_descriptionNames.push_back( descriptionName );
}
}



#endif
