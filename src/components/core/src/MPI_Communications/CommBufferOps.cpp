#include "CommBufferOps.hpp"

namespace geosx
{
localIndex CommBufferOps::Pack( array<char> & buffer, string const & var )
{
  string::size_type sizeOfPackedChars = var.size();

  Pack( buffer, sizeOfPackedChars );

  const char* cvar = var.data();
  for( string::size_type i=0 ; i<sizeOfPackedChars ; ++i )
  {
    buffer.push_back(*(cvar++));
  }

  sizeOfPackedChars += sizeof( localIndex );
  return integer_conversion<localIndex>(sizeOfPackedChars);
}

localIndex CommBufferOps::Pack( char*& buffer,  const std::string& var )
{
  string::size_type sizeOfPackedChars = var.size();

  Pack( buffer, sizeOfPackedChars );

  for( string::size_type i=0 ; i<sizeOfPackedChars ; ++i )
  {
    *buffer = var[i];
    buffer++;
  }

  sizeOfPackedChars += sizeof( localIndex );
  return integer_conversion<localIndex>(sizeOfPackedChars);
}

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


localIndex CommBufferOps::Pack( array<char> & buffer,
                                const array<string>& container )
{
  localIndex sizeOfPackedChars = 0;
  const localIndex arrayLength = container.size();

  sizeOfPackedChars += Pack( buffer, arrayLength );
  for( array<string>::const_iterator i=container.begin() ; i!=container.end() ; ++i )
  {
    sizeOfPackedChars += Pack(buffer,*i);
  }

  return sizeOfPackedChars;
}

localIndex CommBufferOps::Pack( char*& buffer,
                                const array<string>& container )
{
  localIndex sizeOfPackedChars = 0;
  const localIndex arrayLength = container.size();

  sizeOfPackedChars += Pack( buffer, arrayLength );
  for( array<string>::const_iterator i=container.begin() ; i!=container.end() ; ++i )
  {
    sizeOfPackedChars += Pack(buffer,*i);
  }

  return sizeOfPackedChars;
}

localIndex CommBufferOps::Unpack( const char*& buffer,
                               array<string>& array )
{
  array.clear();
  localIndex sizeOfUnpackedChars = 0;

  localIndex arrayLength;
  sizeOfUnpackedChars += Unpack( buffer, arrayLength );
  array.resize(arrayLength);

  for( auto i=0 ; i<arrayLength ; ++i )
  {
    sizeOfUnpackedChars += Unpack( buffer, array[i] );
  }

  return sizeOfUnpackedChars;
}

}
