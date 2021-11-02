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

#ifndef GEOSX_COMMON_FORMAT_HPP_
#define GEOSX_COMMON_FORMAT_HPP_

#if __cplusplus < 202002L
#define GEOSX_USE_FMT
#endif

#ifdef GEOSX_USE_FMT
#include <fmt/core.h>
#include <fmt/chrono.h>
#define GEOSX_FMT_NS ::fmt
#else // use C++20's <format>
#include <format>
#define GEOSX_FMT_NS ::std
#endif

/**
 * @brief Interpolate arguments into a message format string.
 * @param msg the message format string, must be a constant expression
 */
#define GEOSX_FMT( msg, ... ) GEOSX_FMT_NS::format( msg, __VA_ARGS__ )

/**
 * @brief Interpolate arguments into a message format string and write into an output iterator.
 * @param iter the output iterator to write to
 * @param size maximum number of characters to write
 * @param msg the message format string, must be a constant expression
 * @note Ensures the output buffer is zero-terminated (std::format_to_n doesn't)
 */
#define GEOSX_FMT_TO( iter, size, msg, ... ) *GEOSX_FMT_NS::format_to_n( iter, size - 1, msg, __VA_ARGS__ ).out = '\0'

#endif //GEOSX_COMMON_FORMAT_HPP_
