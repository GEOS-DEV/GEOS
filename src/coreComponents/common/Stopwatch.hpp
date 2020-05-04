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

/**
 * @file Stopwatch.hpp
 */

#ifndef GEOSX_STOPWATCH_HPP
#define GEOSX_STOPWATCH_HPP

#include <chrono>

namespace geosx
{
/**
 * @class Stopwatch
 * @brief Class defining a simple stopwatch for interval timing.
 * @note Current implementation relies on std::chrono::steady_clock
 */
class Stopwatch
{
public:
  /**
   * @brief Constructor.  Calls zero().
   */
  Stopwatch()
  {
    zero();
  }

  /**
   * @brief Zero out the timer.
   */
  void zero()
  {
    start = std::chrono::steady_clock::now();
  }

  /**
   * @brief Return elapsed time in seconds since zero() was last called.
   */
  real64 elapsedTime()
  {
    end = std::chrono::steady_clock::now();
    diff = end-start;
    return diff.count();
  }

private:
  std::chrono::steady_clock::time_point start;
  std::chrono::steady_clock::time_point end;
  std::chrono::duration< real64 > diff;
};

} // end geosx

#endif
