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

Section::Section():
  m_rowMinWidth( 20 )
{}

void Section::setName( string_view title )
{
  m_sectionTitle = title;
}

void Section::addDescription( string const & description )
{
  m_descriptions.push_back( description );
}

void Section::addDescription( string const & descriptionName, std::vector< string > const & descriptionValues )
{
  m_descriptionsValues.push_back( descriptionValues );
  m_descriptionNames.push_back( descriptionName );
}

void Section::setMinWidth( integer const & minWidth )
{
  m_rowMinWidth = minWidth;
}

void Section::computeMaxRowSize( string const & title,
                                 std::vector< string > const & rowsDescription )
{
  static constexpr integer marginBorder = 2;
  static constexpr integer nbSpecialChar = 2;
  integer maxDescriptionLength = 0;
  integer titleLength = title.length() + marginBorder * 2 + nbSpecialChar * 2;

  m_rowLength = std::max( m_rowMinWidth, titleLength );

  if( rowsDescription.size() == 0 )
  {
    return;
  }

  auto it = std::max_element( rowsDescription.begin(),
                              rowsDescription.end(),
                              []( const auto & a, const auto & b ) {
    return a.size() < b.size();
  } );

  string maxDescriptionSize = *it;

  maxDescriptionLength = integer( maxDescriptionSize.length()) + marginBorder * 2 + nbSpecialChar * 2;

  m_rowLength = std::max( maxDescriptionLength, m_rowLength );
}

void Section::buildLineSection( string & lineSection )
{
  lineSection =  GEOS_FMT( "{:#>{}}\n", "", m_rowLength );
}

void Section::addTitleRow( string & sectionToBeBuilt, string_view title )
{
  sectionToBeBuilt += GEOS_FMT( "##{:^{}}##\n", title, m_rowLength - 4 );
}

void Section::addEndSectionRow( string & sectionToBeBuilt, string_view title )
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

void Section::clear()
{
  m_descriptions.clear();
  m_sectionTitle.clear();
}

void Section::begin( std::ostream & oss )
{
  string lineSection;
  string sectionToBeBuilt;
  string const titleToDisplay = "Section : " + m_sectionTitle;

  //TODO function and test and refacto below
  if( m_descriptionsValues.empty())
  {
    int maxLenName = 0;
    computeMaxRowSize( "", m_descriptionNames );
    maxLenName = m_rowLength;

    if( m_descriptionsValues.length == 1 )
    {
      m_descriptions.push_back( GEOS_FMT( " - {}:{}", m_descriptionNames[0], m_descriptionsValues[0] ));
    }
    else
    {
      int i = 0;
      for( std::vector< string > descriptionsValues : m_descriptionsValues )
      {
        string description = GEOS_FMT( " - {}:", m_descriptionNames[i] );
        for( string values: descriptionsValues )
        {
          description += GEOS_FMT( "{:-<{}}", values, m_rowLength );
          m_descriptions.push_back();

        }
        i++;
      }
    }
  }

  // check if descvalues empty
  // split descValues into m_descriptions and set marginValues by default 0
  // marginValues is used to construct m_description so just here (?)

  //and back to normal... just need to format m_descriptions !
  //don't forget to rename variable and functions et to remove unused  function !
  computeMaxRowSize( titleToDisplay, m_descriptions );
  buildLineSection( lineSection );

  sectionToBeBuilt += '\n' + lineSection;
  addTitleRow( sectionToBeBuilt, titleToDisplay );
  sectionToBeBuilt += lineSection;
  addDescriptionRows( sectionToBeBuilt, m_descriptions );
  sectionToBeBuilt += '\n';

  oss << sectionToBeBuilt;
}

void Section::end( std::ostream & oss )
{
  string lineSection;
  string sectionToBeBuilt;
  string titleToDisplay = "End : " + m_sectionTitle;

  buildLineSection( lineSection );

  sectionToBeBuilt += '\n';
  addTitleRow( sectionToBeBuilt, titleToDisplay );
  sectionToBeBuilt += lineSection;
  sectionToBeBuilt += '\n';

  oss << sectionToBeBuilt;

  clear();
}
}
