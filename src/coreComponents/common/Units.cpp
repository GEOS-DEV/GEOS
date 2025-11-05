/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file Units.cpp
 */

#include "Units.hpp"
#include <string>
#include "logger/Logger.hpp"

namespace geos
{

namespace units
{


TimeFormatInfo::TimeFormatInfo( double const totalSeconds, int const years, int const days,
                                int const hours, int const minutes, int const seconds ):
  m_totalSeconds( totalSeconds ),
  m_years( years ),
  m_days( days ),
  m_hours( hours ),
  m_minutes( minutes ),
  m_seconds( seconds )
{}

string TimeFormatInfo::toString() const
{
  return GEOS_FMT( "{} ({})", toUnfoldedString(), toSecondsString() );
}
string TimeFormatInfo::toUnfoldedString() const
{
  std::ostringstream oss;
  if( m_totalSeconds < 0.0 )
  {
    oss << "-(";
  }
  if( m_years != 0 )
  {
    oss << m_years << "y, " << m_days << "d, ";
  }
  else if( m_days != 0 )
  {
    oss << m_days << "d, ";
  }
  oss << GEOS_FMT( "{:0>2}h{:0>2}m{:0>2}s", m_hours, m_minutes, m_seconds );
  if( m_totalSeconds < 0.0 )
  {
    oss << ")";
  }
  return oss.str();
}
string TimeFormatInfo::toSecondsString() const
{
  return GEOS_FMT( "{} s", m_totalSeconds );
}

std::ostream & operator<<( std::ostream & os, TimeFormatInfo const & info )
{
  os << info.toString();
  return os;
}

TimeFormatInfo TimeFormatInfo::fromSeconds( double const seconds )
{
  GEOS_LOG( "*** seconds = "+std::to_string(seconds) );
  double remainingSeconds = seconds < 0.0 ? -seconds : seconds;
  GEOS_LOG( "*** remainingSeconds = "+std::to_string(remainingSeconds) );
  GEOS_LOG( "*** YearSeconds = "+std::to_string(YearSeconds) );
  int const totalYears = int( remainingSeconds / YearSeconds );
  GEOS_LOG( "*** totalYears = "+std::to_string(totalYears) );
  remainingSeconds -= totalYears * YearSeconds;
  GEOS_LOG( "*** remainingSeconds = "+std::to_string(remainingSeconds) );
  GEOS_LOG( "*** DaySeconds = "+std::to_string(DaySeconds) );
  int const daysOut = int( remainingSeconds / DaySeconds );
  GEOS_LOG( "*** daysOut = "+std::to_string(daysOut) );
  remainingSeconds -= daysOut * DaySeconds;
  GEOS_LOG( "*** remainingSeconds = "+std::to_string(remainingSeconds) );
  GEOS_LOG( "*** HourSeconds = "+std::to_string(HourSeconds) );
  int const hoursOut = int( remainingSeconds / HourSeconds );
  GEOS_LOG( "*** hoursOut = "+std::to_string(hoursOut) );
  remainingSeconds -= hoursOut * HourSeconds;
  GEOS_LOG( "*** remainingSeconds = "+std::to_string(remainingSeconds) );
  GEOS_LOG( "*** MinuteSeconds = "+std::to_string(MinuteSeconds) );
  int const minutesOut = int( remainingSeconds / MinuteSeconds );
  GEOS_LOG( "*** minutesOut = "+std::to_string(minutesOut) );
  remainingSeconds -= minutesOut * MinuteSeconds;
  GEOS_LOG( "*** remainingSeconds = "+std::to_string(remainingSeconds) );
  int const secondsOut = int( remainingSeconds );
  GEOS_LOG( "*** secondsOut = "+std::to_string(secondsOut) );

  return TimeFormatInfo( seconds, totalYears, daysOut, hoursOut, minutesOut, secondsOut );
}


} // end namespace units

} // end namespace geos
