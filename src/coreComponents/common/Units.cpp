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
 * @file Units.cpp
 */

#include "Units.hpp"

namespace geos
{

namespace units
{


string formatLongDuration( std::chrono::system_clock::duration duration )
{
  using namespace std::chrono;

  const auto hms = duration_cast< seconds >( duration );
  const int64_t totalSeconds = int64_t( hms.count() );
  const int64_t microsecondOnly = duration_cast< microseconds >( duration - hms ).count();

  const int64_t totalHours = totalSeconds / ( 60 * 60 );
  static constexpr int64_t oneYearHours = int64_t( 365.25 * 24 );
  //    - `yDiv.quot` is the years count,
  //    - `yDiv.rem` the year hours count.
  //    - `dDiv.quot` is the year days count,
  //    - `dDiv.rem` the day hours count.
  const auto yDiv = std::div( totalHours, oneYearHours );
  const auto dDiv = std::div( yDiv.rem, int64_t( 24 ));

  std::ostringstream oss;

  if( yDiv.quot != 0 )
  {
    oss << yDiv.quot << "y, ";
  }
  if( dDiv.quot != 0 )
  {
    oss << dDiv.quot << "d, ";
  }
  oss << GEOS_FMT( "{}h{:%Mm%Ss} ({}.{:0>6}s)", dDiv.rem, hms, totalSeconds, microsecondOnly );

  return oss.str();
}


} // end namespace units

} // end namespace geos