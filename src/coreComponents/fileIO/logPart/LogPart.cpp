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
#include <algorithm>

namespace geos
{

LogPart::LogPart( string_view sectionTitle ):
  m_sectionTitle( string( sectionTitle ) ),
  m_sectionWidth( sectionTitle.length() )
{
  m_footerTitle = GEOS_FMT( "End : {}", m_sectionTitle );
}

void LogPart::addDescription( string_view description )
{
  m_descriptions.push_back( string( description ) );
}

void LogPart::addEndDescription( string_view description )
{
  m_endLogMessages.push_back( string( description ) );
}

void LogPart::setMinWidth( integer const & minWidth )
{
  m_rowMinWidth = minWidth;
}

void LogPart::formatAndInsertDescriptions( std::vector< string > & descriptionContainer,
                                           string_view name,
                                           std::vector< string > const & values )
{
  string const nameFormatted =  GEOS_FMT( "- {}: ", string( name ));
  string const descriptionFormatted = GEOS_FMT( "{}{}", nameFormatted, values[0] );
  integer const spacesFromBorder = m_marginBorder * 2 + m_nbBorderChar * 2;
  integer const completeDescriptionLength = descriptionFormatted.length() + spacesFromBorder;
  descriptionContainer.push_back( descriptionFormatted );

  m_sectionWidth = std::max( completeDescriptionLength, m_sectionWidth );

  for( size_t idxValue = 1; idxValue < values.size(); idxValue++ )
  {
    size_t const spaces = values[idxValue].length() + nameFormatted.length();
    descriptionContainer.push_back( GEOS_FMT( "{:>{}}", values[idxValue], spaces ) );
  }
}

string LogPart::buildDescriptionPart( std::vector< string > const & descriptions ) const
{
  std::ostringstream oss;
  for( auto const & description : descriptions )
  {
    integer const remainingLength = m_sectionWidth - m_nbBorderChar * 2 - m_marginBorder;
    string const borderCharacters = GEOS_FMT( "{:#<{}}", "", m_nbBorderChar );
    oss << borderCharacters;
    oss << GEOS_FMT( "{:<{}}{:<{}}", " ", m_marginBorder, description, remainingLength );
    oss << borderCharacters << '\n';
  }
  return oss.str();
}

string LogPart::buildTitlePart( string_view title ) const
{
  std::ostringstream oss;
  integer const titleRowLength = m_sectionWidth - m_nbBorderChar * 2;
  string const borderCharacters =  GEOS_FMT( "{:#<{}}", "", m_nbBorderChar );
  oss <<  GEOS_FMT( "{}{:^{}}{}\n",
                    borderCharacters,
                    title,
                    titleRowLength,
                    borderCharacters );
  return oss.str();
}

void LogPart::begin( std::ostream & os )
{
  m_sectionWidth = std::max( m_sectionWidth, m_rowMinWidth );

  string bottomPart;
  if( !m_descriptions.empty())
  {
    bottomPart = buildDescriptionPart( m_descriptions );
  }

  m_horizontalBorder =  GEOS_FMT( "{:#>{}}\n", "", m_sectionWidth );
  string topPart =  GEOS_FMT( "{}{}{}", m_horizontalBorder,
                              buildTitlePart( m_sectionTitle ),
                              m_horizontalBorder );
  os << GEOS_FMT( "\n{}{}\n", topPart, bottomPart );
}

void LogPart::end( std::ostream & os ) const
{
  string topPart;
  if( !m_endLogMessages.empty() )
  {
    topPart =  GEOS_FMT( "{}{}", buildDescriptionPart( m_endLogMessages ), m_horizontalBorder );
  }
  string bottomPart = GEOS_FMT( "{}{}", buildTitlePart( m_footerTitle ), m_horizontalBorder );
  os << GEOS_FMT( "\n{}{}\n", topPart, bottomPart );
}

}
