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
 * @file InputFlags.hpp
 */

#ifndef GEOSX_DATAREPOSITORY_INPUTFLAGS_HPP_
#define GEOSX_DATAREPOSITORY_INPUTFLAGS_HPP_

#include "common/Logger.hpp"

namespace geosx
{

namespace dataRepository
{

/**
 * @enum InputFlags
 *
 * Enumeration of flags that control reading XML input and schema generation.
 */
enum class InputFlags : int
{
  INVALID,            ///< Invalid value
  FALSE,              ///< Not read from input
  OPTIONAL,           ///< Optional in input
  OPTIONAL_NONUNIQUE, ///< Optional in input, may be repeated
  REQUIRED,           ///< Required in input
  REQUIRED_NONUNIQUE, ///< Required in input, may be repeated
  PROBLEM_ROOT,       ///< Root of the hierarchy
};

/**
 * @brief Convert integer value to InputFlags
 * @param val value to convert
 * @return converted enumeration
 */
inline InputFlags IntToInputFlag( int const val )
{
  InputFlags rval = InputFlags::INVALID;
  switch( val )
  {
    case 0:
    {
      rval = InputFlags::FALSE;
      break;
    }
    case 1:
    {
      rval = InputFlags::OPTIONAL;
      break;
    }
    case 2:
    {
      rval = InputFlags::REQUIRED;
      break;
    }
    default:
    {
      GEOSX_ERROR( "Invalid integer conversion to InputFlag" );
    }
  }
  return rval;
}

/**
 * @brief Convert InputFlags to int
 * @param val value to convert
 * @return converted integer
 */
inline int InputFlagToInt( InputFlags const val )
{
  return static_cast< int >(val);
}

/**
 * @brief Comparison operator for InputFlags enumeration.
 * @param left  lhs value
 * @param right rhs value
 * @return comparison result
 */
inline bool operator==( InputFlags const left, InputFlags const right )
{
  return static_cast< int >(left) == static_cast< int >(right);
}

/**
 * @copydoc operator==(InputFlags const, InputFlags const)
 */
inline bool operator!=( InputFlags const left, InputFlags const right )
{
  return static_cast< int >(left) != static_cast< int >(right);
}

/**
 * @copydoc operator==(InputFlags const, InputFlags const)
 */
inline bool operator<( InputFlags const left, InputFlags const right )
{
  return static_cast< int >(left) < static_cast< int >(right);
}

/**
 * @copydoc operator==(InputFlags const, InputFlags const)
 */
inline bool operator>( InputFlags const left, InputFlags const right )
{
  return static_cast< int >(left) > static_cast< int >(right);
}

/**
 * @copydoc operator==(InputFlags const, InputFlags const)
 */
inline bool operator<=( InputFlags const left, InputFlags const right )
{
  return static_cast< int >(left) <= static_cast< int >(right);
}

/**
 * @copydoc operator==(InputFlags const, InputFlags const)
 */
inline bool operator>=( InputFlags const left, InputFlags const right )
{
  return static_cast< int >(left) >= static_cast< int >(right);
}
}


} /* namespace geosx */



#endif /* GEOSX_DATAREPOSITORY_INPUTFLAGS_HPP_ */
