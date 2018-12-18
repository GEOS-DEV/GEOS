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
  REQUIRED,
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


}


} /* namespace geosx */



#endif /* CORECOMPONENTS_DATAREPOSITORY_INPUTFLAGS_HPP_ */
