/*
 * SchemaFlags.hpp
 *
 *  Created on: Jan 4, 2018
 *      Author: sherman
 */

#ifndef CORECOMPONENTS_FILEIO_SCHEMA_HPP_
#define CORECOMPONENTS_FILEIO_SCHEMA_HPP_

namespace geosx
{

namespace dataRepository
{

enum class SchemaFlags : int
{
  IGNORE,
  ROOT,
  NODE,
  REQUIRED_NODE,
  UNIQUE_NODE,
  REQUIRED_UNIQUE_NODE,
};

inline SchemaFlags IntToSchemaFlag( int const val )
{
  SchemaFlags rval = SchemaFlags::IGNORE;
  switch( val )
  {
    case 0:
    {
      rval = SchemaFlags::IGNORE;
      break;
    }
    case 1:
    {
      rval = SchemaFlags::ROOT;
      break;
    }
    case 2:
    {
      rval = SchemaFlags::NODE;
      break;
    }
    case 3:
    {
      rval = SchemaFlags::REQUIRED_NODE;
      break;
    }
    case 4:
    {
      rval = SchemaFlags::UNIQUE_NODE;
      break;
    }
    case 5:
    {
      rval = SchemaFlags::REQUIRED_UNIQUE_NODE;
      break;
    }
    default:
    {
      GEOS_ERROR( "Invalid integer conversion to InputFlag" );
    }
  }
  return rval;
}

inline int SchemaFlagToInt( SchemaFlags const val )
{
  return static_cast<int>(val);
}

inline bool operator==( SchemaFlags const left, SchemaFlags const right)
{
  return static_cast<int>(left) == static_cast<int>(right);
}

inline bool operator!=( SchemaFlags const left, SchemaFlags const right)
{
  return static_cast<int>(left) != static_cast<int>(right);
}

inline bool operator<( SchemaFlags const left, SchemaFlags const right)
{
  return static_cast<int>(left) < static_cast<int>(right);
}

inline bool operator>( SchemaFlags const left, SchemaFlags const right)
{
  return static_cast<int>(left) > static_cast<int>(right);
}

inline bool operator<=( SchemaFlags const left, SchemaFlags const right)
{
  return static_cast<int>(left) <= static_cast<int>(right);
}

inline bool operator>=( SchemaFlags const left, SchemaFlags const right)
{
  return static_cast<int>(left) >= static_cast<int>(right);
}
}


} /* namespace geosx */



#endif /* CORECOMPONENTS_FILEIO_SCHEMA_HPP_ */
