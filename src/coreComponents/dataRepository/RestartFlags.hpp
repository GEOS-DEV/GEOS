/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef RESTARTFLAGS_H_
#define RESTARTFLAGS_H_

namespace geosx
{
namespace dataRepository
{

enum class RestartFlags : unsigned char
{
  NO_WRITE,
  WRITE,
  WRITE_AND_READ
};

enum class PlotLevel : int
{
  LEVEL_0,
  LEVEL_1,
  LEVEL_2,
  LEVEL_3,
  NOPLOT
};

inline PlotLevel IntToPlotLevel( int const val )
{
  PlotLevel rval = PlotLevel::NOPLOT;
  switch( val )
  {
    case 0:
    {
      rval = PlotLevel::LEVEL_0;
      break;
    }
    case 1:
    {
      rval = PlotLevel::LEVEL_1;
      break;
    }
    case 2:
    {
      rval = PlotLevel::LEVEL_2;
      break;
    }
    case 3:
    {
      rval = PlotLevel::LEVEL_3;
      break;
    }
    default:
    {
      break;
    }
  }
  return rval;
}
inline bool operator<( PlotLevel const left, PlotLevel const right)
{
  return static_cast<int>(left) < static_cast<int>(right);
}

inline bool operator>( PlotLevel const left, PlotLevel const right)
{
  return static_cast<int>(left) > static_cast<int>(right);
}

inline bool operator==( PlotLevel const left, PlotLevel const right)
{
  return static_cast<int>(left) == static_cast<int>(right);
}

inline bool operator!=( PlotLevel const left, PlotLevel const right)
{
  return static_cast<int>(left) != static_cast<int>(right);
}




}   /* namespace dataRepository */
}   /* namespace geosx */

#endif  /* RESTARTFLAGS_H_ */
