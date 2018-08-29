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

#include "BufferOps.hpp"

namespace geosx
{

//localIndex CommBufferOps::PackSize( string const & var )
//{
//  string::size_type sizeOfPackedChars = 0;
//
//  sizeOfPackedChars += PackSize( sizeOfPackedChars );
//  sizeOfPackedChars += var.size();
//
//  return integer_conversion<localIndex>(sizeOfPackedChars);
//}


namespace bufferOps
{
localIndex Unpack( char const *& buffer, string& var )
{
  localIndex sizeOfUnpackedChars = 0;
  string::size_type stringsize = 0;

  sizeOfUnpackedChars += Unpack( buffer, stringsize );

  var.resize(stringsize);
  var.assign( buffer, stringsize );
  buffer += stringsize;
  sizeOfUnpackedChars += stringsize;

  return sizeOfUnpackedChars;
}

}

//localIndex CommBufferOps::PackSize( const string_array& container )
//{
//  localIndex sizeOfPackedChars = 0;
//  const localIndex arrayLength = container.size();
//
//  sizeOfPackedChars += PackSize( arrayLength );
//  for( string_array::const_iterator i=container.begin() ; i!=container.end() ; ++i )
//  {
//    sizeOfPackedChars += PackSize(*i);
//  }
//
//  return sizeOfPackedChars;
//}




//localIndex CommBufferOps::Unpack( const char*& buffer,
//                               string_array& array )
//{
//  array.clear();
//  localIndex sizeOfUnpackedChars = 0;
//
//  localIndex arrayLength;
//  sizeOfUnpackedChars += Unpack( buffer, arrayLength );
//  array.resize(arrayLength);
//
//  for( auto i=0 ; i<arrayLength ; ++i )
//  {
//    sizeOfUnpackedChars += Unpack( buffer, array[i] );
//  }
//
//  return sizeOfUnpackedChars;
//}

}
