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
 * @file LogPart.cpp
 */

#include "LogPart.hpp"
#include "common/format/StringUtilities.hpp"
#include "common/logger/ErrorHandling.hpp"
#include <algorithm>

using namespace geos::stringutilities;
namespace geos
{

LogPart::LogPart( string_view logpartName, bool enableOutput )
{
  m_formattedStartDescription.m_title = logpartName;
  m_formattedEndDescription.m_title = GEOS_FMT( "{}{}", m_prefixEndTitle, logpartName );

  m_enableOutput = enableOutput;
}

void LogPart::addDescription( string_view description )
{
  size_t compareWidth = m_width;
  m_startDescription.m_names.push_back( stringutilities::divideLines< string >( compareWidth, description ) );
  m_startDescription.m_values.push_back( stdVector< string >() );
}

void LogPart::addEndDescription( string_view description )
{
  size_t compareWidth = m_width;
  m_endDescription.m_names.push_back( stringutilities::divideLines< string >( compareWidth, description ) );
  m_endDescription.m_values.push_back( stdVector< string >() );

}


void LogPart::setMinWidth( size_t const & minWidth )
{
  m_minWidth = minWidth;
}

void LogPart::setMaxWidth( size_t const & maxWidth )
{
  m_maxWidth = maxWidth;
}


double clamp( double v, double min, double max )
{
  return LvArray::math::min( max, LvArray::math::max( min, v ));
}

void LogPart::formatDescriptions( LogPart::Description & description,
                                  FormattedDescription & formattedDescription )
{
  stdVector< string > & formattedLines = formattedDescription.m_lines;
  size_t const borderSpaceWidth = m_nbBorderChar * 2 + m_borderMargin * 2;

  size_t const decorationWidth = borderSpaceWidth;
  size_t & maxNameSize = formattedDescription.m_maxNameWidth;
  size_t & maxValueSize = formattedDescription.m_maxValueWidth;

  formattedLines.reserve( description.m_names.size() * 2 );

  m_width = clamp( m_width, m_minWidth, m_maxWidth );

  for( size_t idxName = 0; idxName < description.m_names.size(); idxName++ )
  {
    auto const & rawName =  description.m_names[idxName];
    auto const & rawValues =  description.m_values[idxName];

    // Format name with no values associated
    if( rawValues.empty())
    {
      size_t availableWidth = m_width - borderSpaceWidth;
      auto wrappedNames = stringutilities::wrapTextToMaxLength( rawName, availableWidth );

      for( auto & name : wrappedNames )
      {
        auto const currMaxNameSize = std::max( name.size(), maxNameSize );
        if( currMaxNameSize + decorationWidth < m_width )
        {
          // append space at the end of name if needed
          name.reserve( m_width - borderSpaceWidth );
          name.append( std::string( m_width - currMaxNameSize - decorationWidth, ' ' ));
        }
        formattedLines.push_back( name );
      }
      continue;
    }

    // Format name with values assiociated
    size_t availableWidthForValues = m_width - maxNameSize - decorationWidth - m_delimiter.size();
    auto wrappedValues = stringutilities::wrapTextToMaxLength( rawValues, availableWidthForValues );

    // format name
    stdVector< string > formatNames {rawName};
    for( size_t idxSubName = 0; idxSubName < formatNames.size(); idxSubName++ )
    {
      size_t const paddingNeeded = idxSubName < wrappedValues.size() ?
                                   maxNameSize  - formatNames[idxSubName].size() :
                                   m_width - formatNames[idxSubName].size() - decorationWidth;
      // append space at the end of name if needed
      formatNames[idxSubName].reserve( formatNames[idxSubName].size() + paddingNeeded );
      formatNames[idxSubName].append( paddingNeeded, ' ' );
    }

    size_t const lineCount = std::max( formatNames.size(), wrappedValues.size());

    // format values
    size_t const minValueSizeRequired = m_width - maxNameSize - decorationWidth - m_delimiter.size();
    for( auto & wrappedValue : wrappedValues )
    {
      wrappedValue.reserve( minValueSizeRequired );
      wrappedValue.append( minValueSizeRequired - wrappedValue.size(), ' ' );
    }
    maxValueSize = std::max( maxValueSize,
                             (std::max_element( wrappedValues.begin(), wrappedValues.end() ))->size() );

    // add the first line
    string firstLine;
    firstLine.reserve( formatNames.front().size() + m_delimiter.size() + wrappedValues.front().size());
    firstLine.append( formatNames.front()).append( m_delimiter ).append( wrappedValues.front());
    formattedLines.push_back( firstLine );

    // combination name + value
    for( size_t idxLine = 1; idxLine < lineCount; ++idxLine )
    {
      if( idxLine < formatNames.size() && idxLine < wrappedValues.size())
      { // name + value
        formattedLines.push_back( GEOS_FMT( "{}{}{}", formatNames[idxLine], m_delimiter, wrappedValues[idxLine] ));
      }
      else if( idxLine < formatNames.size())
      { // subnames remaining
        formattedLines.push_back( formatNames[idxLine] );
      }
      else if( idxLine < wrappedValues.size())
      { // subvalues remaining
        size_t const spaceAvailable = maxNameSize + wrappedValues[idxLine].size() + m_delimiter.size();
        formattedLines.push_back( GEOS_FMT( "{:>{}}", wrappedValues[idxLine], spaceAvailable ));
      }
    }
  }
}

string LogPart::outputDescription( FormattedDescription & formattedDescription )
{
  std::ostringstream oss;
  string const borderCharacters = string( m_nbBorderChar, m_borderCharacter );
  string const borderSpaces = string( m_borderMargin, ' ' );

  for( auto const & line : formattedDescription.m_lines )
  {
    oss << borderCharacters;
    oss << borderSpaces << line <<  borderSpaces;
    oss << borderCharacters << '\n';
  }
  return oss.str();
}

string LogPart::outputTitle( LogPart::FormattedDescription & formattedDescription )
{
  size_t const titleRowLength =  m_width;
  string const borderCharacters =  string( m_nbBorderChar, m_borderCharacter );

  return GEOS_FMT( "\n{}{:^{}}{}\n",
                   borderCharacters,
                   formattedDescription.m_title,
                   titleRowLength  - 4,
                   borderCharacters );
}

void LogPart::begin( std::ostream & os )
{
  ErrorLogger::global().setCurrentLogPart( m_formattedStartDescription.m_title );

  if( !m_enableOutput )
    return;

  if( !m_startDescription.m_names.empty())
  {
    formatDescriptions( m_startDescription, m_formattedStartDescription );
  }

  string const line = string( m_width, m_borderCharacter );
  os << '\n' << line;
  os << outputTitle( m_formattedStartDescription );
  os << line << '\n';
  os << outputDescription( m_formattedStartDescription ) << '\n';
}

void LogPart::end( std::ostream & os )
{
  if( !m_enableOutput )
    return;

  formatDescriptions( m_endDescription, m_formattedEndDescription );

  string const line =  string( m_width, m_borderCharacter );
  if( !m_endDescription.m_names.empty() )
  {
    os << '\n';
    os << outputDescription( m_formattedEndDescription );
    os << line;
  }
  os << outputTitle( m_formattedEndDescription );
  os << line << "\n\n";
}

}
