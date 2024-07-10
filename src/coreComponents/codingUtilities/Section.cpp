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
  m_rowMinWidth( 20 )
{}

void Section::addDescription( string_view description )
{
  m_descriptions.push_back( string( description ) );
}

void Section::setMinWidth( integer const & minWidth )
{
  m_rowMinWidth = minWidth;
}

void Section::computeMaxWidth( string_view title,
                               std::vector< string > const & rowsDescription )
{
  string const footerStr = "End : ";
  integer const titleLength = footerStr.length() + title.length() + m_marginBorder * 2 + m_nbSpecialChar * 2;

  m_rowLength = std::max( m_rowMinWidth, titleLength );

  if( rowsDescription.size() == 0 )
  {
    return;
  }

  auto it = std::max_element( rowsDescription.begin(), rowsDescription.end(),
                              []( auto const & a, auto const & b ) {
    return a.size() < b.size();
  } );

  string const maxDescriptionSize = *it;

  integer const maxDescriptionLength = integer( maxDescriptionSize.length()) +
                                       m_marginBorder * 2 + m_nbSpecialChar * 2;

  m_rowLength = std::max( maxDescriptionLength, m_rowLength );
}

void Section::formatAndInsertDescriptions( string_view descriptionName,
                                           std::vector< string > const & descriptionValues )
{
  string const descNameFormatted =  GEOS_FMT( "- {}: ", string( descriptionName ));
  integer const descNameLength = descNameFormatted.length();
  m_descriptions.push_back( GEOS_FMT( "{}{}", descNameFormatted, descriptionValues[0] ) );
  for( size_t idxValue = 1; idxValue < descriptionValues.size(); idxValue++ )
  {
    m_descriptions.push_back( GEOS_FMT( "{:>{}}{}", " ",
                                        descNameLength,
                                        descriptionValues[idxValue] ) );
  }
}

void Section::beginSection( std::ostream & oss )
{
  computeMaxWidth( m_sectionTitle, m_descriptions );

  string const lineSection =  GEOS_FMT( "{:#>{}}\n", "", m_rowLength );
  string const horizontalChars =  GEOS_FMT( "{:#<{}}", "", m_nbSpecialChar );
  integer const titleLength = m_rowLength - m_nbSpecialChar * 2;
  integer const descriptionLength = m_rowLength - m_nbSpecialChar * 2 - m_marginBorder;

  //build section title
  oss << '\n';
  oss << lineSection;
  oss << GEOS_FMT( "{}{:^{}}{}\n", horizontalChars, m_sectionTitle, titleLength, horizontalChars );

  //build section descriptions
  oss << lineSection;

  if( !m_descriptions.empty())
  {
    for( string const & description : m_descriptions )
    {
      oss << horizontalChars;
      oss << GEOS_FMT( "{:<{}}{:<{}}", " ", m_marginBorder, description, descriptionLength );
      oss << horizontalChars;
      oss << '\n';
    }
  }
  oss << "\n";
}

void Section::endSection( std::ostream & oss ) const
{
  string const footerTitle = GEOS_FMT( "{} {}", "End :", m_sectionTitle );
  string const lineSection = GEOS_FMT( "{:#^{}}\n", "", m_rowLength );
  string const horizontalChars =  GEOS_FMT( "{:#<{}}", "", m_nbSpecialChar );
  integer const titleLength = m_rowLength - m_nbSpecialChar * 2;

  oss << '\n';
  oss << GEOS_FMT( "{}{:^{}}{}\n", horizontalChars, footerTitle, titleLength, horizontalChars );
  oss << lineSection;
  oss << '\n';
}
}
