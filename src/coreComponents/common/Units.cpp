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
  std::ostringstream oss;
  if( m_years != 0 )
  {
    oss << m_years << "y, " << m_days << "d, ";
  }
  else if( m_days != 0 )
  {
    oss << m_days << "d, ";
  }
  oss << GEOS_FMT( "{:0>2}h{:0>2}m{:0>2}s ({} s)",
                   m_hours, m_minutes, m_seconds, m_totalSeconds );
  return oss.str();
}

std::ostream & operator<<( std::ostream & os, TimeFormatInfo const & info )
{
  os << info.toString();
  return os;
}


template< typename Duration >
TimeFormatInfo TimeFormatInfo::fromDuration( Duration const value )
{
  using namespace std::chrono;

  auto const totalYears = duration_cast< units::Years >( value );
  auto const daysOut = duration_cast< units::Days >( value - totalYears );
  auto const hoursOut = duration_cast< hours >( value - totalYears - daysOut );
  auto const minutesOut = duration_cast< minutes >( value - totalYears - daysOut - hoursOut );
  auto const secondsOut = duration_cast< seconds >( value - totalYears - daysOut - hoursOut - minutesOut );

  return TimeFormatInfo( duration< double >( value ).count(), int( totalYears.count() ),
                         int( daysOut.count() ), int( hoursOut.count() ),
                         int( minutesOut.count() ), int( secondsOut.count() ) );
}
// available specializations
template TimeFormatInfo TimeFormatInfo::fromDuration< SystemClock::duration >( SystemClock::duration duration );

TimeFormatInfo TimeFormatInfo::fromSeconds( double const seconds )
{
  int totalYears = int(   seconds / YearSeconds );
  int daysOut = int(    ( seconds - totalYears * YearSeconds ) / DaySeconds );
  int hoursOut = int(   ( seconds - totalYears * YearSeconds - daysOut * DaySeconds ) / HourSeconds );
  int minutesOut = int( ( seconds - totalYears * YearSeconds - daysOut * DaySeconds - hoursOut * HourSeconds ) / MinuteSeconds );
  int secondsOut = int(   seconds - totalYears * YearSeconds - daysOut * DaySeconds - hoursOut * HourSeconds - minutesOut * MinuteSeconds );

  return TimeFormatInfo( seconds, totalYears, daysOut, hoursOut, minutesOut, secondsOut );
}


} // end namespace units

} // end namespace geos
