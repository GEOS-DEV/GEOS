/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

/*
 * InputFlags.hpp
 *
 *  Created on: Dec 17, 2018
 *      Author: settgast
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
  return static_cast<int>(val);
}

inline bool operator==( InputFlags const left, InputFlags const right)
{
  return static_cast<int>(left) == static_cast<int>(right);
}

inline bool operator!=( InputFlags const left, InputFlags const right)
{
  return static_cast<int>(left) != static_cast<int>(right);
}

inline bool operator<( InputFlags const left, InputFlags const right)
{
  return static_cast<int>(left) < static_cast<int>(right);
}

inline bool operator>( InputFlags const left, InputFlags const right)
{
  return static_cast<int>(left) > static_cast<int>(right);
}

inline bool operator<=( InputFlags const left, InputFlags const right)
{
  return static_cast<int>(left) <= static_cast<int>(right);
}

inline bool operator>=( InputFlags const left, InputFlags const right)
{
  return static_cast<int>(left) >= static_cast<int>(right);
}
}


} /* namespace geosx */



#endif /* CORECOMPONENTS_DATAREPOSITORY_INPUTFLAGS_HPP_ */
