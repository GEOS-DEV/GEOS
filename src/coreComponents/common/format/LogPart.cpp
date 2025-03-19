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
  m_startDescription.m_names.push_back( stringutilities::divideLines< string >( description ) );
  m_startDescription.m_values.push_back( std::vector< string >() );
}

void LogPart::addEndDescription( string_view description )
{
  m_endDescription.m_names.push_back( stringutilities::divideLines< string >( description ) );
  m_endDescription.m_values.push_back( std::vector< string >() );
}


void LogPart::setMinWidth( size_t const & minWidth )
{
  m_startDescription.m_minWidth = minWidth;
  m_endDescription.m_minWidth = minWidth;
}

void LogPart::setMaxWidth( size_t const & maxWidth )
{
  m_startDescription.m_maxWidth = maxWidth;
  m_endDescription.m_maxWidth = maxWidth;
}

void LogPart::formatDescriptions( LogPart::Description & description,
                                  FormattedDescription & formattedDescription )
{
  std::vector< string > & formattedLines = formattedDescription.m_lines;
  size_t const borderSpaceWidth = m_nbBorderChar * 2 + m_borderMargin * 2;
  for( size_t idxName = 0; idxName < description.m_names.size(); idxName++ )
  {
    std::vector< string > & names =  description.m_names[idxName];
    std::vector< string > & values =  description.m_values[idxName];

    // if no values process only the names
    if( values.empty())
    {
      auto wrappedNames = wrapTextToMaxLength( names, description.m_maxWidth - borderSpaceWidth );
      formattedLines.insert( formattedLines.end(), wrappedNames.begin(), wrappedNames.end());
      continue;
    }

    //format name
    for( auto & name: names )
    {
      string const spaces = std::string( formattedDescription.m_maxNameWidth- name.size(), ' ' );
      name = GEOS_FMT( "{}{}", name, spaces );
    }

    //format values
    size_t const valueSpaceAvailable = description.m_maxWidth - formattedDescription.m_maxNameWidth - borderSpaceWidth - m_delimiter.size();
    auto wrappedValues = wrapTextToMaxLength( values, valueSpaceAvailable );

    size_t maxValueWidth = (*std::max_element( wrappedValues.begin(), wrappedValues.end(),
                                               []( string_view a, string_view b ) {
      return a.length() < b.length();
    } )).size();
    // maxValueWidth = 0;
    size_t const totalLineWidth = formattedDescription.m_maxNameWidth + maxValueWidth + borderSpaceWidth + m_delimiter.size();
    m_width = std::max( m_width, totalLineWidth );

    //2.5 merge both name and values
    size_t const lineCount = std::max( names.size(), wrappedValues.size());
    formattedLines.push_back( GEOS_FMT( "{}{}{}", names[0], m_delimiter, wrappedValues[0] ));
    if( formattedLines[0].size() < formattedDescription.m_maxNameWidth )
    {
      string const spaces = std::string( formattedDescription.m_maxNameWidth- formattedLines[0].size(), ' ' );
      formattedLines[0] = GEOS_FMT( "{}{}", formattedLines[0], spaces );
    }
    for( size_t idxLine = 1; idxLine < lineCount; ++idxLine )
    {
      if( idxLine < names.size() && idxLine < wrappedValues.size())
      {
        formattedLines.push_back( GEOS_FMT( "{}{}{}", names[idxLine], m_delimiter, wrappedValues[idxLine] ));
      }
      else if( idxLine < names.size())
      {
        formattedLines.push_back( names[idxLine] );
      }
      else if( idxLine < wrappedValues.size())
      {
        size_t const spaceAvailable = formattedDescription.m_maxNameWidth + wrappedValues[idxLine].size() + m_delimiter.size();
        formattedLines.push_back( GEOS_FMT( "{:>{}}", wrappedValues[idxLine], spaceAvailable ));
      }

      if( formattedLines[idxLine].size() < formattedDescription.m_maxNameWidth )
      {
        string const spaces = std::string( formattedDescription.m_maxNameWidth - formattedLines[idxLine].size(), ' ' );
        formattedLines[idxLine] = GEOS_FMT( "{}{}", formattedLines[idxLine], spaces );
      }
      m_width = std::max( m_width, formattedLines[idxLine].size() );
    }

  }
  m_width = std::max( m_width, description.m_minWidth );
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

void LogPart::computeInitialLogWidth( LogPart::Description & description,
                                      FormattedDescription & formattedDescription )
{

  if( !description.m_names.empty() )
  {
    size_t maxStringSize = 0;
    for( size_t idxDescription = 0; idxDescription  < description.m_names.size(); idxDescription++ )
    {
      if( !description.m_values[idxDescription].empty())
      {
        for( auto const & name : description.m_names[idxDescription] )
        {
          maxStringSize = std::max( maxStringSize, name.size());
        }
      }
    }

    formattedDescription.m_maxNameWidth= maxStringSize;
  }
}

void LogPart::begin( std::ostream & os )
{
  if( MpiWrapper::commRank() != 0 )
    return;

  computeInitialLogWidth( m_startDescription, m_formattedStartDescription );

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

  computeInitialLogWidth( m_endDescription, m_formattedEndDescription );

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
