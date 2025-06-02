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

#ifndef GEOS_COMMON_FORMAT_LOGPART_HPP
#define GEOS_COMMON_FORMAT_LOGPART_HPP

#include "common/DataTypes.hpp"
#include "common/format/Format.hpp"
#include "common/format/StringUtilities.hpp"

namespace geos
{

using namespace stringutilities;
/**
 * @brief Class for displaying section for different steps of simulation
 */
class LogPart
{
public:

  /**
   * @brief Initialize a LogPart given a title
   * @param logPartTitle The title who will be used for top and bottom LogPart
   * @param enableOutput Boolean to activate or not csv output
   */
  LogPart( string_view logPartTitle, bool enableOutput );

  /**
   * @brief Add a description to the top LogPart
   * @param name The first part of the description
   * @param args The remaining part of the description,all remaining values will be concatened and aligned
   * @note Descriptions values can be be any formatted types.
   */
  template< typename ... Args >
  void addDescription( string_view name, Args const & ... args );

  /**
   * @brief Add a description to the top logPart. No specific format will be apply the this description
   * @param description The string value of the description.
   */
  void addDescription( string_view description );

  /**
   * @brief Add a description to the bottom logPart by concatening a description name and descriptions values.
   * @param name The first part of the description
   * @param args The remaining part of the description, all remaining values will be concatened and aligned
   * @note Descriptions values can be be any formatted types.
   */
  template< typename ... Args >
  void addEndDescription( string_view name, Args const & ... args );

  /**
   * @brief Add a description to the top logPart. No specific format will be apply the this description
   * @param description The string value of the description
   */
  void addEndDescription( string_view description );

  /**
   * @brief Set the minimal width of the LogPart
   * @param minWidth The minimal width to apply
   */
  void setMinWidth( size_t const & minWidth );

  /**
   * @brief Set the maximal width of the LogPart
   * @param maxWidth The maximal width to apply
   */
  void setMaxWidth( size_t const & maxWidth );

  /**
   * @brief Draw the first part of the logPart. It include the title and optionnaly the description(s).
   * @param os An output stream (by default, std::cout)
   */
  void begin( std::ostream & os = std::cout );

  /**
   * @brief Draw the last part of the logPart. It include the title and optionnaly the end description(s).
   * @param oss An output stream (by default, std::cout)
   */
  void end( std::ostream & oss = std::cout );

  /**
   * @brief Toggles the CSV output feature.
   * @param enabled Boolean to activate csv output
   */
  void enableOutput( bool enabled )
  { m_enableOutput = enabled; }

private:

  /**
   * @brief Structure containing all information needed in order to construct a top or bottom LogPart.
   * All these variables can be adjusted.
   */
  struct Description
  {
    /// Name of the description (first part of a description), it can be splited by \\n
    std::vector< std::vector< string > > m_names;
    /// Values in the description (remaining part of a description),
    /// each vector of values is associated with one name
    std::vector< std::vector< string > > m_values;
  };

  /**
   * @brief Structure containing all the information needed to display a logPart.
   */
  struct FormattedDescription
  {
    /// Log part title
    string m_title;
    /// Vector containing the descriptions formatted by formatDescriptions()
    std::vector< string > m_lines;
    /// max length name (first part of a description) of a logPart
    size_t m_maxNameWidth;
    /// max length name (remaining part of a description) of a logPart
    size_t m_maxValueWidth;
  };

  Description m_startDescription = { {}, {} };
  Description m_endDescription  = { {}, {} };

  FormattedDescription m_formattedStartDescription = {  "", {}, 0, 0 };
  FormattedDescription m_formattedEndDescription  = { "", {}, 0, 0 };

  /// logPart default length
  size_t m_width = 100;
  /// minimal length of a log part
  size_t m_minWidth = 100;
  /// maximal length of a log part
  size_t m_maxWidth = SIZE_MAX;
  /// margin (left and right) between all descriptions and the log part borders
  static constexpr size_t m_borderMargin = 2;
  /// numbers of character used for the border
  static constexpr size_t m_nbBorderChar = 2;
  /// character used for border
  char const m_borderCharacter = '#';
  /// prefix to append to the title of bottom section
  static constexpr string_view m_prefixEndTitle = "End of ";
  /// string used to separate the name/description
  static constexpr string_view m_delimiter = " : ";
  /// Active the LogPart output
  bool m_enableOutput = true;

  /**
   * @brief Add a description to a specific section (top or bottom)
   * @param description Structure containing all the information (name, values, length) needed for building a logPart
   * @param name The first part of the description
   * @param args The remaining part of the description, all remaining values will be concatened and aligned
   * @note Descriptions values can be be any formatted types. Values will be aligned altogether.
   */
  template< typename ... Args >
  void addDescriptionBySection( Description & description, FormattedDescription & formattedDescription,
                                string_view name, Args const &... args );

  /**
   * @brief Construct a formatted description from the all the descriptions set in input
   * @param description Structure containing all the information (name, values, length) needed for building a logPart
   * @param formattedDescription Structure containing the formatted description
   */
  void formatDescriptions( LogPart::Description & description,
                           FormattedDescription & formattedDescription );

  /**
   * @param formattedDescription Structure containing the formatted description
   * @return A formatted string containing the log part title
   */
  string outputTitle( FormattedDescription & formattedDescription );

  /**
   * @param formattedDescription Structure containing the formatted description
   * @return A formatted string containing the log part descriptions
   */
  string outputDescription( FormattedDescription & formattedDescription );
};

template< typename ... Args >
void LogPart::addDescriptionBySection( Description & description, FormattedDescription & formattedDescription,
                                       string_view name, Args const &... args )
{
  std::vector< string > values;
  size_t & maxValueSize = formattedDescription.m_maxValueWidth;
  size_t & maxNameSize = formattedDescription.m_maxNameWidth;
  ( [&] {
    static_assert( has_formatter_v< decltype(args) >,
                   "Argument passed cannot be converted to string" );
    string const value = GEOS_FMT( "{}", args );

    std::vector< string_view > splitValues =  divideLines< string_view >( maxValueSize, value );
    values.insert( values.end(), splitValues.begin(), splitValues.end() );
  } (), ...);

  description.m_values.push_back( values );

  size_t lineWidth = 0;
  std::vector< string > nameDivided = divideLines< string >( lineWidth, name );
  if( lineWidth == 0 )
    lineWidth = name.size();
  maxNameSize = std::max( maxNameSize, lineWidth );

  description.m_names.push_back( nameDivided );

  size_t const formattingCharSize = m_nbBorderChar * 2 + m_borderMargin * 2;
  size_t const currentTotalWidth =  maxNameSize + maxValueSize + formattingCharSize;
  m_width = std::max( m_width, currentTotalWidth );
  m_width = std::max( m_width, formattedDescription.m_title.size());
}

template< typename ... Args >
void LogPart::addDescription( string_view name, Args const &... args )
{
  addDescriptionBySection( m_startDescription, m_formattedStartDescription, name, args ... );
}

template< typename ... Args >
void LogPart::addEndDescription( string_view name, Args const &... args )
{
  addDescriptionBySection( m_endDescription, m_formattedEndDescription, name, args ... );
}

}

#endif
