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
  m_rowMinWidth( 70 )
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

void Section::computeWidth()
{
  integer const titleLength = m_footerTitle.length() + m_sectionTitle.length();
  integer maxDescriptionLength = titleLength;

  if( !m_descriptions.empty())
  { 
    auto it = std::max_element( m_descriptions.begin(), m_descriptions.end(),
                                []( auto const & a, auto const & b ) {
      return a.size() < b.size();
    } );
    string const maxDescriptionSize = *it;
    maxDescriptionLength = std::max( (integer)maxDescriptionSize.length(), maxDescriptionLength );
  }

  if( !m_endLogMessages.empty())
  {
    auto it = std::max_element( m_endLogMessages.begin(), m_endLogMessages.end(),
                                []( auto const & a, auto const & b ) {
      return a.size() < b.size();
    } );
    string const maxDescriptionSize = *it;
    maxDescriptionLength = std::max( (integer)maxDescriptionSize.length(), maxDescriptionLength );
  }

  m_sectionWidth = maxDescriptionLength + m_marginBorder * 2 + m_nbSpecialChar * 2;
}

void Section::formatAndInsertDescriptions( std::vector< string > & descriptionContainer,
                                           string_view descriptionName,
                                           std::vector< string > const & descriptionValues )
{
  string const descNameFormatted =  GEOS_FMT( "- {}: ", string( descriptionName ));
  integer const descNameLength = descNameFormatted.length();
  descriptionContainer.push_back( GEOS_FMT( "{}{}", descNameFormatted, descriptionValues[0] ) );
  for( size_t idxValue = 1; idxValue < descriptionValues.size(); idxValue++ )
  {
    descriptionContainer.push_back( GEOS_FMT( "{:>{}}{}", " ",
                                              descNameLength,
                                              descriptionValues[idxValue] ) );
  }
}

string Section::constructDescriptionsWithContainer( std::vector< string > const & descriptions, integer const length ) const
{
  std::ostringstream oss;
  for( const auto & description : descriptions )
  {
    constructDescription( oss, description, length );
  }
  return oss.str();
}

void Section::constructDescription( std::ostringstream & oss, string const & description, integer length ) const
{
  oss << m_borderSpaces;
  oss << GEOS_FMT( "{:<{}}{:<{}}", " ", m_marginBorder, description, length );
  oss << m_borderSpaces << '\n';
}

void Section::beginSection( std::ostream & oss )
{
  computeWidth();

  m_lineSection =  GEOS_FMT( "{:#>{}}\n", "", m_sectionWidth );
  m_borderSpaces =  GEOS_FMT( "{:#<{}}", "", m_nbSpecialChar );

  integer const titleLength = m_sectionWidth - m_nbSpecialChar * 2;
  integer const descriptionLength = m_sectionWidth - m_nbSpecialChar * 2 - m_marginBorder;
  string descriptions;

  if( !m_descriptions.empty())
  {
    descriptions = constructDescriptionsWithContainer( m_descriptions, descriptionLength );
  }

  //TODO fonction output (args)
  oss << '\n';
  oss << m_lineSection;
  oss << GEOS_FMT( "{}{:^{}}{}\n", m_borderSpaces, m_sectionTitle, titleLength, m_borderSpaces );//TODO refacto here
  oss << m_lineSection;
  oss << descriptions;
  oss << '\n';
}

void Section::endSection( std::ostream & oss ) const
{
  string const footerTitle = GEOS_FMT( "{}{}", m_footerTitle, m_sectionTitle );
  string const lineSection = GEOS_FMT( "{:#^{}}\n", "", m_sectionWidth );
  string const horizontalChars =  GEOS_FMT( "{:#<{}}", "", m_nbSpecialChar );

  integer const titleLength = m_sectionWidth - m_nbSpecialChar * 2;
  integer const descriptionLength = m_sectionWidth - m_nbSpecialChar * 2 - m_marginBorder;
  string endDescriptions;

  if( !m_endLogMessages.empty() )
  {
    endDescriptions = constructDescriptionsWithContainer( m_endLogMessages, descriptionLength );
  }

  //TODO fonction output (args)
  oss << '\n';
  oss << GEOS_FMT( "{}{:^{}}{}\n", horizontalChars, footerTitle, titleLength, horizontalChars );
  oss << lineSection;
  oss << endDescriptions;
  oss << '\n';
}
}
