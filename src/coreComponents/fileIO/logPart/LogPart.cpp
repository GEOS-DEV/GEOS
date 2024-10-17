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
{}

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
                                           string_view descriptionName,
                                           std::vector< string > const & descriptionValues )
{
  string const descriptionNameFormatted =  GEOS_FMT( "- {}: ", string( descriptionName ));
  string const descriptionFormatted = GEOS_FMT( "{}{}", descriptionNameFormatted, descriptionValues[0] );
  integer const spacesFromBorder = m_marginBorder * 2 + m_nbBorderChar * 2;
  integer const completeDescriptionLength = descriptionFormatted.length() + spacesFromBorder;
  descriptionContainer.push_back( descriptionFormatted );

  m_sectionWidth = std::max( completeDescriptionLength, m_sectionWidth );

  for( size_t idxValue = 1; idxValue < descriptionValues.size(); idxValue++ )
  {
    size_t const spaces = descriptionValues[idxValue].length() + descriptionNameFormatted.length();
    descriptionContainer.push_back( GEOS_FMT( "{:>{}}", descriptionValues[idxValue], spaces ) );
  }
}

string LogPart::constructDescriptionsFromVector( std::vector< string > const & descriptions ) const
{
  std::ostringstream oss;
  for( const auto & description : descriptions )
  {
    integer const rowDescriptionLength = m_sectionWidth - m_nbBorderChar * 2 - m_marginBorder;
    formatDescription( oss, description, rowDescriptionLength );
  }
  return oss.str();
}

void LogPart::formatDescription( std::ostringstream & oss,
                                 string_view description,
                                 integer const remainingLength ) const
{
  string const borderCharacters = GEOS_FMT( "{:#<{}}", "", m_nbBorderChar );
  oss << borderCharacters;
  oss << GEOS_FMT( "{:<{}}{:<{}}", " ", m_marginBorder, description, remainingLength );
  oss << borderCharacters << '\n';
}

void LogPart::begin( std::ostream & os )
{
  m_sectionWidth = std::max( m_sectionWidth, m_rowMinWidth );
  integer const titleRowLength = m_sectionWidth - m_nbBorderChar * 2;

  string descriptions;
  if( !m_descriptions.empty())
  {
    descriptions = constructDescriptionsFromVector( m_descriptions );
  }
  string const borderCharacters =  GEOS_FMT( "{:#<{}}", "", m_nbBorderChar );
  string titleRowFormatted =  GEOS_FMT( "{}{:^{}}{}\n",
                                        borderCharacters,
                                        m_sectionTitle,
                                        titleRowLength,
                                        borderCharacters );

  m_horizontalBorder =  GEOS_FMT( "{:#>{}}\n", "", m_sectionWidth );
  string const topPart = GEOS_FMT( "{}{}{}", m_horizontalBorder, titleRowFormatted, m_horizontalBorder );
  string const bottomPart = descriptions;
  os << GEOS_FMT( "\n{}{}\n", topPart, bottomPart );
}

void LogPart::end( std::ostream & os ) const
{
  string topPart;
  if( !m_endLogMessages.empty() )
  {
    topPart =  GEOS_FMT( "{}{}", constructDescriptionsFromVector( m_endLogMessages ), m_horizontalBorder );
  }

  string const borderCharacters =  GEOS_FMT( "{:#<{}}", "", m_nbBorderChar );
  string const footerTitle = GEOS_FMT( "{}{}", m_footerTitle, m_sectionTitle );
  integer const titleRowLength = m_sectionWidth - m_nbBorderChar * 2;
  string titleRowFormatted = GEOS_FMT( "{}{:^{}}{}\n",
                                       borderCharacters,
                                       footerTitle,
                                       titleRowLength,
                                       borderCharacters );
  string const bottomPart = GEOS_FMT( "{}{}", titleRowFormatted, m_horizontalBorder );
  os << GEOS_FMT( "\n{}{}\n", topPart, bottomPart );
}

}
