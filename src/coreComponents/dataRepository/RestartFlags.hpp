/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_DATAREPOSITORY_RESTARTFLAGS_HPP_
#define GEOSX_DATAREPOSITORY_RESTARTFLAGS_HPP_

/**
 * @file RestartFlags.hpp
 */

namespace geosx
{
namespace dataRepository
{

/**
 * @enum RestartFlags
 *
 * A scoped enum for the restart options.
 */
enum class RestartFlags : unsigned char
{
  NO_WRITE,      ///< Doe no write into restart
  WRITE,         ///< Write into restart
  WRITE_AND_READ ///< Write and read from restart
};

/**
 * @enum PlotLevel
 *
 * A scoped enum for the Plot options.
 */
enum class PlotLevel : int
{
  LEVEL_0, ///< Write to plot always
  LEVEL_1, ///< Write to plot when plotLevel>=1 is specified in input
  LEVEL_2, ///< Write to plot when plotLevel>=2 is specified in input
  LEVEL_3, ///< Write to plot when plotLevel>=3 is specified in input
  NOPLOT   ///< Do not ever write to plot file
};

/**
 * @brief Function to get a PlotLevel enum from an int
 * @param val int that represents the PlotLevel
 * @return The PlotLevel that corresponds to the input
 */
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

/**
 * @brief boolean operator< for PlotLevel
 * @param[in] left
 * @param[in] right
 * @return boolean result of left<right
 */
inline bool operator<( PlotLevel const left, PlotLevel const right )
{
  return static_cast< int >(left) < static_cast< int >(right);
}

/**
 * @brief boolean operator> for PlotLevel
 * @param[in] left
 * @param[in] right
 * @return boolean result of left>right
 */
inline bool operator>( PlotLevel const left, PlotLevel const right )
{
  return static_cast< int >(left) > static_cast< int >(right);
}

/**
 * @brief boolean operator== for PlotLevel
 * @param[in] left
 * @param[in] right
 * @return boolean result of left==right
 */
inline bool operator==( PlotLevel const left, PlotLevel const right )
{
  return static_cast< int >(left) == static_cast< int >(right);
}

/**
 * @brief boolean operator!= for PlotLevel
 * @param[in] left
 * @param[in] right
 * @return boolean result of left!=right
 */
inline bool operator!=( PlotLevel const left, PlotLevel const right )
{
  return static_cast< int >(left) != static_cast< int >(right);
}



}   /* namespace dataRepository */
}   /* namespace geosx */

#endif  /* GEOSX_DATAREPOSITORY_RESTARTFLAGS_HPP_ */
