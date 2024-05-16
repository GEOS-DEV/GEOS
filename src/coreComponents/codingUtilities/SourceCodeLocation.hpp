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
 * @file SourceCodeLocation.hpp
 *
 * Collection of utilities to facilitate I/O of enumeration types.
 * Provides a macro definition that allows associating string names
 * with enumeration constants and a set of functions that make use
 * of these strings, like stream insertion/extraction operators.
 */

#ifndef GEOS_SOURCECODELOCATION_HPP
#define GEOS_SOURCECODELOCATION_HPP

#include <string_view>
#include <iostream>
#include <type_traits>
#include <algorithm>

#if __cplusplus >= 202002L
#include <source_location>
#endif

namespace geos
{


/**
 * @brief Instanciated with GEOS_SRCLOC macro, this class behave as c++20 std::source_location,
 * but has a limitation: It cannot be used as a default parameter of a function as std::source_location
 * can be.
 */
class SourceCodeLocation
{
public:
  constexpr SourceCodeLocation( const std::string_view fileName, const size_t line ) noexcept:
    m_fileName( fileName ),
    m_line( line )
  {}

  constexpr std::string_view getFileName() const noexcept
  { return m_fileName; }

  constexpr size_t getLine() const noexcept
  { return m_line; }
  
  friend std::ostream & operator<<( std::ostream & os, const SourceCodeLocation & loc )
  { os << loc.m_fileName << ":" << loc.m_line; return os; }

private:
  std::string_view m_fileName;
  size_t m_line;
};


#if __cplusplus < 202002L

/**
 * @brief Instanciate a SourceCodeLocation using compiler macros.
 */
#define GEOS_SRCLOC() SourceCodeLocation( __FILE__, __LINE__ )

#else

/**
 * @brief Instanciate a SourceCodeLocation using c++20 std::source_location.
 * @warning To prevent returning a wrong line number, do not insert line breaks in the macro.
 */
#define GEOS_SRCLOC() SourceCodeLocation( std::source_location::current().file_name(), std::source_location::current().line() )

#endif


} // namespace geos

#endif /* GEOS_SOURCECODELOCATION_HPP */
