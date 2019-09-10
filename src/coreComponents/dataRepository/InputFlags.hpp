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

#ifndef CORECOMPONENTS_DATAREPOSITORY_INPUTFLAGS_HPP_
#define CORECOMPONENTS_DATAREPOSITORY_INPUTFLAGS_HPP_

namespace geosx
{

namespace dataRepository
{

enum class InputFlags : int
{
  INVALID,
  FALSE,
  OPTIONAL,
  OPTIONAL_NONUNIQUE,
  REQUIRED,
  REQUIRED_NONUNIQUE,
  PROBLEM_ROOT,
};

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
      GEOS_ERROR( "Invalid integer conversion to InputFlag" );
    }
  }
  return rval;
}

inline int InputFlagToInt( InputFlags const val )
{
  return static_cast< int >(val);
}

inline bool operator==( InputFlags const left, InputFlags const right )
{
  return static_cast< int >(left) == static_cast< int >(right);
}

inline bool operator!=( InputFlags const left, InputFlags const right )
{
  return static_cast< int >(left) != static_cast< int >(right);
}

inline bool operator<( InputFlags const left, InputFlags const right )
{
  return static_cast< int >(left) < static_cast< int >(right);
}

inline bool operator>( InputFlags const left, InputFlags const right )
{
  return static_cast< int >(left) > static_cast< int >(right);
}

inline bool operator<=( InputFlags const left, InputFlags const right )
{
  return static_cast< int >(left) <= static_cast< int >(right);
}

inline bool operator>=( InputFlags const left, InputFlags const right )
{
  return static_cast< int >(left) >= static_cast< int >(right);
}
}


} /* namespace geosx */



#endif /* CORECOMPONENTS_DATAREPOSITORY_INPUTFLAGS_HPP_ */
