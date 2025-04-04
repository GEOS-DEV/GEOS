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
#include <algorithm>

using namespace geos::stringutilities;
namespace geos
{

LogPart::LogPart( string_view logPartTitle )
{
  m_formattedStartDescription.m_title = logPartTitle;
  m_formattedEndDescription.m_title = GEOS_FMT( "{}{}", m_prefixEndTitle, logPartTitle );
}

void LogPart::addDescription( string_view description )
{
  size_t compareWidth = m_width;
  m_startDescription.m_names.push_back( stringutilities::divideLines< string >( compareWidth, description ) );
  m_startDescription.m_values.push_back( std::vector< string >() );
  m_width = std::max( compareWidth, m_width );
}

void LogPart::addEndDescription( string_view description )
{
  size_t compareWidth = m_width;
  m_endDescription.m_names.push_back( stringutilities::divideLines< string >( compareWidth, description ) );
  m_endDescription.m_values.push_back( std::vector< string >() );
  m_width = std::max( compareWidth, m_width );

}


void LogPart::setMinWidth( size_t const & minWidth )
{
  m_minWidth = minWidth;
}

void LogPart::setMaxWidth( size_t const & maxWidth )
{
  m_maxWidth = maxWidth;
}

void LogPart::formatDescriptions( LogPart::Description & description,
                                  FormattedDescription & formattedDescription )
{
  std::vector< string > & formattedLines = formattedDescription.m_lines;
  size_t const borderSpaceWidth = m_nbBorderChar * 2 + m_borderMargin * 2;
  size_t const formattingCharSize = borderSpaceWidth + m_delimiter.size();

  m_width = std::min( m_width, m_maxWidth );
  m_width = std::max( m_width, m_minWidth );

  for( size_t idxName = 0; idxName < description.m_names.size(); idxName++ )
  {
    std::vector< string > const & nonFormattedNames =  description.m_names[idxName];
    std::vector< string > const & nonFormattedValues =  description.m_values[idxName];

    if( nonFormattedValues.empty())
    {
      size_t maxLineLength = m_maxWidth - borderSpaceWidth;
      auto wrappedNames = stringutilities::wrapTextToMaxLength( nonFormattedNames, maxLineLength );
      formattedLines.insert( formattedLines.end(), wrappedNames.begin(), wrappedNames.end());
      continue;
    }

    //format name
    std::vector< string > formatNames {nonFormattedNames};
    size_t formattedNameWidth = formattedDescription.m_maxNameWidth;
    for( auto it=formatNames.begin(); it<formatNames.end(); it++ )
    {
      // found a name too large
      if( it->size() > formattedNameWidth )
      {
        std::vector< string > nameToResize{*it};
        auto resizedNames = stringutilities::wrapTextToMaxLength( nameToResize, formattedNameWidth );
        formatNames.insert( it, resizedNames.begin(), resizedNames.end());
        formatNames.pop_back();
        it=formatNames.begin();
      }
    }

    for( auto & name : formatNames )
    {
      string const spaces = std::string( formattedNameWidth  - name.size(), ' ' );
      name = GEOS_FMT( "{}{}", name, spaces );
    }

    //format values
    size_t spacesForValues = m_maxWidth - formattedNameWidth - formattingCharSize;
    auto wrappedValues = stringutilities::wrapTextToMaxLength( nonFormattedValues, spacesForValues );

    size_t const lineCount = std::max( formatNames.size(), wrappedValues.size());
    formattedLines.push_back( GEOS_FMT( "{}{}{}", formatNames.front(), m_delimiter, wrappedValues.front() ));
    if( formattedLines.front().size() < formattedNameWidth )
    {
      string const spaces = std::string( formattedNameWidth- formattedLines.front().size(), ' ' );
      formattedLines.front() = GEOS_FMT( "{}{}", formattedLines.front(), spaces );
    }
    for( size_t idxLine = 1; idxLine < lineCount; ++idxLine )
    {
      if( idxLine < formatNames.size() && idxLine < wrappedValues.size())
      {
        formattedLines.push_back( GEOS_FMT( "{}{}{}", formatNames[idxLine], m_delimiter, wrappedValues[idxLine] ));
      }
      else if( idxLine < formatNames.size())
      {
        formattedLines.push_back( formatNames[idxLine] );
      }
      else if( idxLine < wrappedValues.size())
      {
        size_t const spaceAvailable = formattedNameWidth + wrappedValues[idxLine].size() + m_delimiter.size();
        formattedLines.push_back( GEOS_FMT( "{:>{}}", wrappedValues[idxLine], spaceAvailable ));
      }

      if( formattedLines[idxLine].size() < formattedNameWidth )
      {
        string const spaces = std::string( formattedNameWidth - formattedLines[idxLine].size(), ' ' );
        formattedLines[idxLine] = GEOS_FMT( "{}{}", formattedLines[idxLine], spaces );
      }
    }
  }
}

string LogPart::outputDescription( FormattedDescription & formattedDescription )
{
  std::ostringstream oss;
  string const borderCharacters = string( m_nbBorderChar, m_borderCharacter );

  for( auto const & line : formattedDescription.m_lines )
  {
    oss << borderCharacters;
    oss << GEOS_FMT( "{:<{}}{:<{}}", " ", m_borderMargin,
                     line, m_width - m_nbBorderChar * 2 - m_borderMargin );
    oss << borderCharacters << '\n';
  }
  return oss.str();
}

string LogPart::outputTitle( LogPart::FormattedDescription & formattedDescription )
{
  size_t const titleRowLength = m_width;
  string const borderCharacters =  string( m_nbBorderChar, m_borderCharacter );

  return GEOS_FMT( "\n{}{:^{}}{}\n",
                   borderCharacters,
                   formattedDescription.m_title,
                   titleRowLength  - 4,
                   borderCharacters );
}

void LogPart::begin( std::ostream & os )
{
  if( MpiWrapper::commRank() != 0 )
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
  if( MpiWrapper::commRank() != 0 )
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
