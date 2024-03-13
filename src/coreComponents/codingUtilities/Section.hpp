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
 * @file Section.hpp
 */

#ifndef GEOS_COMMON_SECTION_HPP
#define GEOS_COMMON_SECTION_HPP

#include "common/DataTypes.hpp"

namespace geos
{

class Section
{
public:

  Section();

  /**
   * @brief Set the name of the section
   * @param m_sectionTitle The name of the section
   */
  void setName( string_view m_sectionTitle );

  /**
   * @brief Add a description to the section
   * @param description The string value of the description
   */
  void addDescription( string const & description );

  /**
   * @brief Set the minimal width of a row
   * @param minWidth the minimal width of the table
   */
  void setMinWidth( integer const & minWidth );

  /**
   * @brief Draw the first part of the section. It include the title and optionnaly, the description(s);
   * @param os An output stream (by default, std::cout)
   */
  void begin( std::ostream & os = std::cout );

  /**
   * @brief Draw the last part of the section. It include the title
   * @param oss An output stream (by default, std::cout)
   */
  void end( std::ostream & oss = std::cout );
private:

  /**
   * @brief Compute the max string size (m_rowLength) between title and the description(s)
   * @param m_sectionTitle The title of the table
   * @param vDescriptions The vector of descriptions
   */
  void computeMaxRowSize( string const & m_sectionTitle,
                          std::vector< string > const & vDescriptions );

  /**
   * @brief Build the line section in order to build the section
   * @param lineSection An empty string
   */
  void buildLineSection( string & lineSection );

  /**
   * @brief Build and add the title to the first part of the section
   * @param sectionToBeBuilt The current section being built
   * @param title The section name
   */
  void addTitleRow( string & sectionToBeBuilt, string const & title );

  /**
   * @brief Build and add the title to the last part of the section
   * @param sectionToBeBuilt The current section being built
   * @param title The section name
   */
  void addEndSectionRow( string & sectionToBeBuilt, string const & title );

  /**
   * @brief Build and add the descriptions to the first part of the section
   * @param sectionToBeBuilt The current section being built
   * @param rowsValue The vector of descriptions
   */
  void addDescriptionRows( string & sectionToBeBuilt, std::vector< string > const & rowsValue );

  std::vector< string > m_vDescriptions;

  string m_sectionTitle;
  integer m_rowLength;
  integer m_rowMinWidth;
};
}



#endif
