#include "CommBufferOps.hpp"

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



localIndex CommBufferOps::Unpack( char const *& buffer, string& var )
{
  localIndex sizeOfUnpackedChars = 0;
  string::size_type stringsize = 0;

  sizeOfUnpackedChars += Unpack( buffer, stringsize );

  var.assign( buffer, stringsize );
  buffer += stringsize;
  sizeOfUnpackedChars += stringsize;

  return sizeOfUnpackedChars;
}



//localIndex CommBufferOps::PackSize( const array<string>& container )
//{
//  localIndex sizeOfPackedChars = 0;
//  const localIndex arrayLength = container.size();
//
//  sizeOfPackedChars += PackSize( arrayLength );
//  for( array<string>::const_iterator i=container.begin() ; i!=container.end() ; ++i )
//  {
//    sizeOfPackedChars += PackSize(*i);
//  }
//
//  return sizeOfPackedChars;
//}




//localIndex CommBufferOps::Unpack( const char*& buffer,
//                               array<string>& array )
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
