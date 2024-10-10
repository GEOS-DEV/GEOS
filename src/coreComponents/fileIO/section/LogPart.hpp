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
 * @file LogPart.hpp
 */

#ifndef GEOS_COMMON_SECTION_HPP
#define GEOS_COMMON_SECTION_HPP

#include "common/DataTypes.hpp"
#include "common/format/Format.hpp"

namespace geos
{

class LogPart
{
public:

  /**
   * @brief Construct a new LogPart
   * @param m_sectionTitle The section title
   */
  LogPart( string_view m_sectionTitle );

  /**
   * @brief Add a description to the section by concatening a description name and descriptions values.
   * @param descriptionName The description name
   * @param args Descriptions values to be concatened.
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
   * @brief Add a description to the end of the section by concatening a description name and descriptions values.
   * @param descriptionName The description name
   * @param args Descriptions values to be concatened.
   * Descriptions values can be be any types and will be aligned
   */
  template< typename ... Args >
  void addEndDescription( string_view descriptionName, Args const & ... args );

  /**
   * @brief Add a description to the end of the section
   * @param description The string value of the description
   */
  void addEndDescription( string_view description );

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
   * @brief Build a description from the name and description values
   * @param descriptionName The decription name
   * @param decriptionsValues The description values
   */
  void formatAndInsertDescriptions( std::vector< string > & descriptionContainer,
                                    string_view descriptionName,
                                    std::vector< string > const & decriptionsValues );

  /**
   * @brief Constructs a string from a vector of descriptions.
   * @param descriptions A vector of string containing the descriptions
   * @return A string containing all the formatted descriptions.
   */
  string constructDescriptionsFromVector( std::vector< string > const & descriptions ) const;

  /**
   * @brief Formats a description and adds it to an output stream
   * @param oss The output stream
   * @param description The description to be formatted
   * @param remainingLength The remaining length available for the description
   * after accounting for borders and margins
   */
  void formatDescription( std::ostringstream & oss,
                          string const & description,
                          integer const remainingLength ) const;

  /// Vector containing all descriptions
  std::vector< string > m_descriptions;
  /// title of section
  std::vector< string > m_endLogMessages;

  /// title of section
  string m_sectionTitle;
  /// section length
  integer m_sectionWidth;
  /// min width of section length
  integer m_rowMinWidth = 70;

  /// description border margin
  static constexpr integer m_marginBorder = 2;
  /// numbers of character used as border
  static constexpr integer m_nbSpecialChar = 2;

  /// String containing horizontal border
  string horizontalBorder;

  /// Start title footer string
  static string_view constexpr m_footerTitle = "End : ";

};

template< typename ... Args >
void LogPart::addDescription( string_view descriptionName, Args const &... args )
{
  std::vector< string > descriptionsValues;
  ( [&] {
    static_assert( has_formatter_v< decltype(args) >, "Argument passed in addRow cannot be converted to string" );
    string const value = GEOS_FMT( "{}", args );
    descriptionsValues.push_back( value );
  } (), ...);

  formatAndInsertDescriptions( m_descriptions, descriptionName, descriptionsValues );
}

template< typename ... Args >
void LogPart::addEndDescription( string_view descriptionName, Args const &... args )
{
  std::vector< string > descriptionsValues;
  ( [&] {
    static_assert( has_formatter_v< decltype(args) >, "Argument passed in addRow cannot be converted to string" );
    string const value = GEOS_FMT( "{}", args );
    descriptionsValues.push_back( value );
  } (), ...);

  formatAndInsertDescriptions( m_endLogMessages, descriptionName, descriptionsValues );
}
}

#endif
