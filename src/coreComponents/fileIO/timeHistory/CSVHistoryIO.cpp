#include <sstream>

#include "CSVHistoryIO.hpp"
#include "ASCIIFile.hpp"
#include "common/MpiWrapper.hpp"

namespace geos
{

namespace detail
{

template <typename T>
void inline printAdvance(std::ostream & os, buffer_unit_type * & data)
{
  os << *reinterpret_cast<T*>(data);
  data += sizeof(T);
}

void inline printAndAdvance(std::ostream & os, std::type_index const & type, buffer_unit_type * & data)
{
  static const std::unordered_map<std::type_index, std::function<void(std::ostream&, buffer_unit_type*&)>> actions =
  {
    { std::type_index(typeid(char)), &printAdvance<char> },
    { std::type_index(typeid(signed char)), &printAdvance<signed char> },
    { std::type_index(typeid(real32)), &printAdvance<real32> },
    { std::type_index(typeid(real64)), &printAdvance<real64> },
    { std::type_index(typeid(integer)), &printAdvance<integer> },
    { std::type_index(typeid(localIndex)), &printAdvance<localIndex> },
    { std::type_index(typeid(globalIndex)), &printAdvance<globalIndex> }
  };

  auto it = actions.find(type);
  if(it != actions.end())
  {
    it->second(os, data);
  }
  else
  {
    GEOS_ERROR( "Error: unsupported type for CSV Time History output!" );
  }
}

void inline printRankedData( std::ostream & outFile,
                              buffer_unit_type *& dataBuffer,
                              size_t currentDim,
                              size_t dataRank,
                              const std::vector< size_t > & dataDims,
                              std::type_index atomDatatype )
{
  if ( currentDim < dataRank - 1 )
  {
    for ( size_t ii = 0; ii < dataDims[currentDim]; ++ii )
    {
      printRankedData(outFile, dataBuffer, currentDim + 1, dataRank, dataDims, atomDatatype);
      if (currentDim == dataRank - 2)
      {
        outFile << " ";
      }
    }
  }
  else
  {
    printAndAdvance( outFile, atomDatatype, dataBuffer );
  }
}

void inline printBufferedDataRows( std::ostream & outFile,
                                   buffer_unit_type * & dataBuffer,
                                   size_t rowCount,
                                   const std::vector< size_t > & colCounts,
                                   size_t dataRank,
                                   const std::vector< size_t > & dataDims,
                                   std::type_index atomDatatype )
{
  for ( size_t row = 0; row < rowCount; ++row )
  {
    for( size_t col = 0; col < colCounts[row]; ++col )
    {
      printRankedData( outFile, dataBuffer, 0, dataRank, dataDims, atomDatatype);
      outFile << ", ";
    }
    outFile << "\n";
  }
}

}

CSVHistoryIO::CSVHistoryIO( string const & filename,
                            localIndex dataRank,
                            std::vector< localIndex > const & dims,
                            std::type_index typeIdx,
                            localIndex writeHead,
                            MPI_Comm comm ):
  BufferedHistoryIO( typeIdx, dataRank, dims ),
  m_filename( filename ),
  m_writeHead( writeHead ),
  m_comm( comm )
{
  int rank = MpiWrapper::commRank( m_comm );
  if (m_filename.substr(m_filename.size() - 4) == ".csv")
  {
      m_filename = m_filename.substr(0, m_filename.size() - 4)
                + "." + std::to_string(rank) + ".csv";
  }
  else
  {
      m_filename += "." + std::to_string(rank) + ".csv";
  }
}

void CSVHistoryIO::init( bool existsOkay )
{
  ASCIIFile outfile( m_filename, false, existsOkay );
  GEOS_ERROR_IF( outfile.countLinesInFile( ) != m_writeHead,
                 "Error: Unexpected number of lines in existing CSV time history output!" );
}

void CSVHistoryIO::write()
{
  if(m_bufferedCount > 0)
  {
    buffer_unit_type * dataBuffer = nullptr;
    if(m_dataBuffer.size() > 0)
    {
      dataBuffer = &m_dataBuffer[0];
    }
    else
    {
      GEOS_ERROR("Error: Attempting to write buffered TimeHistory data but no data is buffered." );
    }

    // Use a stringstream to buffer the output
    std::stringstream bufferStream;

    detail::printBufferedDataRows(bufferStream, dataBuffer, m_bufferedCount, m_bufferedLocalIdxCounts, m_rank, m_dims, m_typeIdx);

    ASCIIFile asciiFile(m_filename, false, true);
    // Append the buffered data to the ASCII file
    asciiFile.append(bufferStream.str());
  }
  emptyBuffer();
}

}
