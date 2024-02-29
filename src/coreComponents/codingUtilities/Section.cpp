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

namespace geos
{

void Section::setName( string_view title )
{
  m_sectionTitle = title;
}

void Section::addDescription( string const & description )
{
  m_vDescriptions.push_back( description );
}

void Section::computeMaxRowSize( string const & title,
                                 std::vector< string > const & rowsDescription )
{
  integer marginBorder = 2;
  integer nbSpecialChar = 2;
  integer maxDescriptionLength = 0;
  integer titleLength = title.length() + marginBorder * 2 + nbSpecialChar * 2;

  if( rowsDescription.size() == 0 )
  {
    m_rowLength = titleLength;
    return;
  }

  auto it = std::max_element( rowsDescription.begin(),
                              rowsDescription.end(),
                              []( const auto & a, const auto & b ) {
    return a.size() < b.size();
  } );

  string maxDescriptionSize = *it;

  maxDescriptionLength = integer( maxDescriptionSize.length()) + marginBorder * 2 + nbSpecialChar * 2;

  m_rowLength = maxDescriptionLength > titleLength ? maxDescriptionLength : titleLength;
}

void Section::buildLineSection( string & lineSection )
{
  lineSection =  GEOS_FMT( "{:#>{}}\n", "", m_rowLength );
}

void Section::addTitleRow( string & sectionToBeBuilt, string const & title )
{
  sectionToBeBuilt += GEOS_FMT( "##{:^{}}##\n", title, m_rowLength - 4 );
}

void Section::addEndSectionRow( string & sectionToBeBuilt, string const & title )
{
  sectionToBeBuilt += GEOS_FMT( "##{:^{}}##\n", title, m_rowLength - 4 );
}

void Section::addDescriptionRows( string & sectionToBeBuilt, std::vector< string > const & rowValues )
{
  for( string rowValue : rowValues )
  {
    sectionToBeBuilt += GEOS_FMT( "##  {:<{}}##\n", rowValue, m_rowLength - 6 );
  }
}

void Section::begin( std::ostream & oss )
{
  string lineSection;
  string sectionToBeBuilt;
  string titleToAdd = "Section : " + m_sectionTitle;

  computeMaxRowSize( titleToAdd, m_vDescriptions );
  buildLineSection( lineSection );

  sectionToBeBuilt += '\n' + lineSection;
  addTitleRow( sectionToBeBuilt, titleToAdd );
  sectionToBeBuilt += lineSection;
  addDescriptionRows( sectionToBeBuilt, m_vDescriptions );
  sectionToBeBuilt += '\n';

  oss << sectionToBeBuilt;
}

void Section::end( std::ostream & oss )
{
  string lineSection;
  string sectionToBeBuilt;
  string titleToAdd = "End : " + m_sectionTitle;

  buildLineSection( lineSection );

  sectionToBeBuilt += '\n';
  addTitleRow( sectionToBeBuilt, titleToAdd );
  sectionToBeBuilt += lineSection;
  sectionToBeBuilt += '\n';

  oss << sectionToBeBuilt;
}
}
