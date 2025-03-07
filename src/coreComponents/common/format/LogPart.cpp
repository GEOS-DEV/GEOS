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

LogPart::LogPart( string_view logPartTitle, integer commRank )
{
  m_startDesc.m_title = logPartTitle;
  m_endDesc.m_title = GEOS_FMT( "{}{}", m_prefixEndTitle, logPartTitle );
  m_rank = commRank;
}

void LogPart::addDescription( string const & description )
{
  m_startDesc.m_names.push_back( splitStringByNewLine( description ) );
  m_startDesc.m_values.push_back( std::vector< string >() );
}

void LogPart::addEndDescription( string const & description )
{
  m_endDesc.m_names.push_back( splitStringByNewLine( description ) );
  m_endDesc.m_values.push_back( std::vector< string >() );
}


void LogPart::setMinWidth( size_t const & minWidth )
{
  m_startDesc.m_logPartMinWidth = minWidth;
  m_endDesc.m_logPartMinWidth = minWidth;
}

void LogPart::setMaxWidth( size_t const & maxWidth )
{
  m_startDesc.m_logPartMaxWidth = maxWidth;
  m_endDesc.m_logPartMaxWidth = maxWidth;
}

std::vector< std::string > LogPart::splitStringByNewLine( std::string const & str )
{
  auto result = std::vector< std::string >{};
  auto ss = std::stringstream{str};

  for( std::string line; std::getline( ss, line, '\n' ); )
    result.push_back( line );

  return result;
}

void LogPart::formatDescriptions( LogPart::Description & description )
{
  size_t const logPartMaxWidth = description.m_logPartMaxWidth;
  size_t const & logPartMinWidth = description.m_logPartMinWidth;
  size_t const & logPartMaxNameWidth =  description.m_logPartMaxNameWidth;
  std::vector< string > & formattedLines = description.m_formattedDescriptionLines;
  size_t const borderSpaceWidth = m_nbBorderChar * 2 + m_borderMargin * 2;
  for( size_t idxName = 0; idxName < description.m_names.size(); idxName++ )
  {
    std::vector< string > & names =  description.m_names[idxName];
    std::vector< string > & values =  description.m_values[idxName];

    // if no values process only the names
    if( values.empty())
    {
      auto wrappedNames = wrapTextToMaxLength( names, logPartMaxWidth - borderSpaceWidth );
      formattedLines.insert( formattedLines.end(), wrappedNames.begin(), wrappedNames.end());
      continue;
    }

    //format name
    for( auto & name: names )
    {
      string const spaces = std::string( logPartMaxNameWidth- name.size(), ' ' );
      name = GEOS_FMT( "{}{}", name, spaces );
    }

    //format values
    size_t const valueSpaceAvailable = logPartMaxWidth - logPartMaxNameWidth - borderSpaceWidth - m_delimiter.size();
    auto wrappedValues = wrapTextToMaxLength( values, valueSpaceAvailable );

    size_t const maxValueWidth = (*max_element( wrappedValues.begin(), wrappedValues.end(),
                                                []( string const & a, string const & b ) {
      return a.length() < b.length();
    } )).size();
    size_t const totalLineWidth = logPartMaxNameWidth + maxValueWidth + borderSpaceWidth + m_delimiter.size();
    description.m_logPartWidth = std::max( description.m_logPartWidth, totalLineWidth );

    //2.5 merge both names and values
    size_t const lineCount = std::max( names.size(), wrappedValues.size());
    formattedLines.push_back( GEOS_FMT( "{}{}{}", names[0], m_delimiter, wrappedValues[0] ));
    if( formattedLines[0].size() < logPartMaxNameWidth )
    {
      string const spaces = std::string( logPartMaxNameWidth- formattedLines[0].size(), ' ' );
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
        size_t const spaceAvailable = logPartMaxNameWidth + wrappedValues[idxLine].size() + m_delimiter.size();
        formattedLines.push_back( GEOS_FMT( "{:>{}}", wrappedValues[idxLine], spaceAvailable ));
      }

      if( formattedLines[idxLine].size() < logPartMaxNameWidth )
      {
        string const spaces = std::string( logPartMaxNameWidth - formattedLines[idxLine].size(), ' ' );
        formattedLines[idxLine] = GEOS_FMT( "{}{}", formattedLines[idxLine], spaces );
      }
      description.m_logPartWidth = std::max( description.m_logPartWidth, formattedLines[idxLine].size() );
    }

  }
  description.m_logPartWidth = std::max( description.m_logPartWidth, logPartMinWidth );
}

string LogPart::buildDescriptionPart( LogPart::Description const & description )
{
  std::ostringstream oss;
  for( auto const & formattedDescription : description.m_formattedDescriptionLines )
  {
    string const borderCharacters = string( m_nbBorderChar, m_borderCharacter );
    oss << borderCharacters;
    oss << GEOS_FMT( "{:<{}}{:<{}}", " ", m_borderMargin,
                     formattedDescription, description.m_logPartWidth - m_nbBorderChar * 2 - m_borderMargin );
    oss << borderCharacters << '\n';
  }
  return oss.str();
}

string LogPart::buildTitlePart( LogPart::Description const & description )
{
  std::ostringstream oss;
  size_t const titleRowLength = description.m_logPartWidth;
  string const borderCharacters =  string( m_nbBorderChar, m_borderCharacter );
  oss << GEOS_FMT( "{}{:^{}}{}\n",
                   borderCharacters,
                   description.m_title,
                   titleRowLength  - 4,
                   borderCharacters );
  return oss.str();
}

void LogPart::computeInitialLogWidth( LogPart::Description & description,
                                      std::vector< std::vector< string > > const & descriptionNames,
                                      std::vector< std::vector< string > > m_values )
{

  if( !descriptionNames.empty() )
  {
    size_t maxStringSize = 0;
    for( size_t idxDescription = 0; idxDescription  < descriptionNames.size(); idxDescription++ )
    {
      if( !m_values[idxDescription].empty())
      {
        maxStringSize = std::max( maxStringSize,
                                  (*max_element( descriptionNames[idxDescription].begin(), descriptionNames[idxDescription].end(),
                                                 []( string const & a, string const & b ) {
          return a.length() < b.length();
        } )).size());
      }
    }

    description.m_logPartMaxNameWidth= maxStringSize;
  }
}

void LogPart::begin( std::ostream & os )
{
  if( m_rank != -1 && MpiWrapper::commRank() != m_rank )
    return;

  string bottomPart = "";
  auto & descriptionNames = m_startDesc.m_names;

  computeInitialLogWidth( m_startDesc, descriptionNames, m_startDesc.m_values );

  if( !descriptionNames.empty())
  {
    formatDescriptions( m_startDesc );
  }

  bottomPart = buildDescriptionPart( m_startDesc );

  string const line = string( m_startDesc.m_logPartWidth, m_borderCharacter );
  string const topPart =  GEOS_FMT( "{}\n{}{}\n",
                                    line,
                                    buildTitlePart( m_startDesc ),
                                    line );
  os << GEOS_FMT( "\n{}{}\n", topPart, bottomPart );
}

void LogPart::end( std::ostream & os )
{
  if( m_rank != -1 && MpiWrapper::commRank() != m_rank )
    return;

  string topPart = "";
  auto & descriptionNames = m_endDesc.m_names;

  computeInitialLogWidth( m_endDesc, descriptionNames, m_endDesc.m_values );

  formatDescriptions( m_endDesc );

  string const line =  string( m_endDesc.m_logPartWidth, m_borderCharacter );
  if( !descriptionNames.empty() )
  {
    topPart = GEOS_FMT( "{}{}\n", buildDescriptionPart( m_endDesc ), line );
  }

  string const bottomPart = GEOS_FMT( "{}{}\n", buildTitlePart( m_endDesc ), line );
  os << GEOS_FMT( "\n{}{}\n", topPart, bottomPart );
}

}
