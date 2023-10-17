/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_FILEIO_TIMEHISTORY_ASCIIFILE_HPP_
#define GEOS_FILEIO_TIMEHISTORY_ASCIIFILE_HPP_

#include <fstream>
#include <string>
#include <filesystem>
#include <algorithm>
#include <iterator>

namespace geos
{

class ASCIIFile
{
public:
  // Constructor: Opens the ASCII file. If the file already exists and `existsOkay` is false, throws an error.
  ASCIIFile( const std::string & filename,
             bool deleteExisting,
             bool existsOkay) :
    m_filename(filename),
    m_filestream()
  {
    bool exists = std::filesystem::exists(m_filename);
    if ( exists )
    {
      if ( !existsOkay )
      {
        GEOS_ERROR("File already exists: " + m_filename);
      }
      else if ( deleteExisting )
      {
        m_filestream.open(m_filename,std::ios::trunc);
      }
    }
    else
    {
      m_filestream.open(m_filename,std::ios::app);
      if (!m_filestream)
      {
        GEOS_ERROR("Failed to create file: " + m_filename);
      }
      m_filestream.close();
    }
  }

  // Destructor: Closes the file
  ~ASCIIFile()
  {
    if (m_filestream.is_open())
    {
      m_filestream.close();
    }
  }

  // Append data to the ASCII file
  void append( const std::string &data )
  {
    ensureOpenForAppend();
    m_filestream << data;
  }

  // Counts the number of lines in the file
  int countLinesInFile()
  {
    std::ifstream inFile(m_filename);
    if (!inFile)
    {
      GEOS_ERROR( "Error: Cannot open the file " + m_filename );
    }

    return std::count(std::istreambuf_iterator<char>(inFile),
                      std::istreambuf_iterator<char>(), '\n');
  }

private:
  // Ensures the file stream is open for appending data
  void ensureOpenForAppend()
  {
    if (!m_filestream.is_open())
    {
      m_filestream.open(m_filename, std::ios::app);
      if (!m_filestream)
      {
        GEOS_ERROR("Failed to open file for appending: " + m_filename);
      }
    }
  }

  std::string m_filename;
  std::ofstream m_filestream;
};


}

#endif