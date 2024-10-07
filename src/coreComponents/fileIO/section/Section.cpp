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
 * @file Section.cpp
 */

#include "Section.hpp"
#include <algorithm>

namespace geos
{

Section::Section( string_view sectionTitle ):
  m_sectionTitle( string( sectionTitle ) ),
  m_sectionWidth( sectionTitle.length() )
{}

void Section::addDescription( string_view description )
{
  m_descriptions.push_back( string( description ) );
}

void Section::addEndDescription( string_view description )
{
  m_endLogMessages.push_back( string( description ) );
}

void Section::setMinWidth( integer const & minWidth )
{
  m_rowMinWidth = minWidth;
}

void Section::formatAndInsertDescriptions( std::vector< string > & descriptionContainer,
                                           string_view descriptionName,
                                           std::vector< string > const & descriptionValues )
{
  string descriptionNameFormatted =  GEOS_FMT( "- {}: ", string( descriptionName ));
  string const descriptionFormatted = GEOS_FMT( "{}{}", descriptionNameFormatted, descriptionValues[0] );
  integer spacesFromBorder = m_marginBorder * 2 + m_nbSpecialChar * 2;
  integer const completeDescriptionLength = descriptionFormatted.length() + spacesFromBorder;
  descriptionContainer.push_back( descriptionFormatted );

  m_sectionWidth = std::max( completeDescriptionLength, m_sectionWidth );

  for( size_t idxValue = 1; idxValue < descriptionValues.size(); idxValue++ )
  {
    size_t const spaces = descriptionValues[idxValue].length() + descriptionNameFormatted.length();
    descriptionContainer.push_back( GEOS_FMT( "{:>{}}", descriptionValues[idxValue], spaces ) );
  }
}

string Section::constructDescriptionsWithContainer( std::vector< string > const & descriptions ) const
{
  std::ostringstream oss;
  for( const auto & description : descriptions )
  {
    integer const rowDescriptionLength = m_sectionWidth - m_nbSpecialChar * 2 - m_marginBorder;
    constructDescription( oss, description, rowDescriptionLength );
  }
  return oss.str();
}

void Section::constructDescription( std::ostringstream & oss,
                                    string const & description,
                                    integer const remainingLength ) const
{
  string borderCharacters = GEOS_FMT( "{:#<{}}", "", m_nbSpecialChar );
  oss << borderCharacters;
  oss << GEOS_FMT( "{:<{}}{:<{}}", " ", m_marginBorder, description, remainingLength );
  oss << borderCharacters << '\n';
}

void Section::outputSection( std::ostream & oss, string_view topPart, string_view bottomPart ) const
{
  oss << '\n';
  oss << topPart;
  oss << bottomPart;
  oss << '\n';
}

void Section::beginSection( std::ostream & oss )
{
  m_sectionWidth = std::max( m_sectionWidth, m_rowMinWidth );
  m_lineSection =  GEOS_FMT( "{:#>{}}\n", "", m_sectionWidth );
  string const borderCharacters =  GEOS_FMT( "{:#<{}}", "", m_nbSpecialChar );

  integer const titleLength = m_sectionWidth - m_nbSpecialChar * 2;
  string descriptions;

  if( !m_descriptions.empty())
  {
    descriptions = constructDescriptionsWithContainer( m_descriptions );
  }

  string titleRowFormatted =  GEOS_FMT( "{}{:^{}}{}\n",
                                        borderCharacters,
                                        m_sectionTitle,
                                        titleLength,
                                        borderCharacters );

  string topPart = GEOS_FMT( "{}{}{}", m_lineSection, titleRowFormatted, m_lineSection );
  string bottomPart = descriptions;
  outputSection( oss, topPart, bottomPart );
}

void Section::endSection( std::ostream & os ) const
{
  string const footerTitle = GEOS_FMT( "{}{}", m_footerTitle, m_sectionTitle );
  string const lineSection = GEOS_FMT( "{:#^{}}\n", "", m_sectionWidth );
  string const borderCharacters =  GEOS_FMT( "{:#<{}}", "", m_nbSpecialChar );
  integer const titleLength = m_sectionWidth - m_nbSpecialChar * 2;
  string topPart;

  if( !m_endLogMessages.empty() )
  {
    topPart =  GEOS_FMT( "{}{}", constructDescriptionsWithContainer( m_endLogMessages ), m_lineSection );
  }

  string titleRowFormatted = GEOS_FMT( "{}{:^{}}{}\n",
                                       borderCharacters,
                                       footerTitle,
                                       titleLength,
                                       borderCharacters );
  string bottomPart = GEOS_FMT( "{}{}", titleRowFormatted, m_lineSection );
  outputSection( os, topPart, bottomPart );
}
}
