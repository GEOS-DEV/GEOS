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
 * @file IOUtilities.h
 * @author walsh24
 * @date Oct 13, 2011
 */

#ifndef IOUTILITIES_H_
#define IOUTILITIES_H_


#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

namespace geosx
{

class IOUtilities
{
public:
  IOUtilities();
  virtual ~IOUtilities();

  template< typename T >
  static void parse_file( array1d<T> & target, std::string filename, char delimiter );
};

template< typename T >
void IOUtilities::parse_file( array1d<T> & target, std::string filename, char delimiter )
{
  std::ifstream inputStream(filename.c_str());
  std::string lineString;
  T value;

  if (inputStream)
  {
    while (!inputStream.eof())
    {
      std::getline(inputStream, lineString);
      std::istringstream ss( lineString );

      while(ss.peek() == delimiter || ss.peek() == ' ')
      {
        ss.ignore();
      }
      while( ss>>value )
      {
        target.push_back( value );
        while(ss.peek() == delimiter || ss.peek() == ' ')
        {
          ss.ignore();
        }
      }
    }

    inputStream.close();
  }
  else
  {
    GEOS_ERROR("Could not read input file!");
  }
}

}


/*
   /// Read from a character deliminated file into a vector
   template<class ARRAY>
   void dlmreadVector(const std::string& filename, ARRAY& values, char delim = '
      ', int skipLines = 0)
   {

   std::ifstream inputStream(filename.c_str());
   if (inputStream)
   {
    std::string lineString;

    for (int i = 0; i < skipLines; ++i)
      std::getline(inputStream, lineString);

    // read every line from the stream
    while (std::getline(inputStream, lineString))
    {
      std::istringstream lineStream(lineString);
      std::string cellString;
      // read every element from the line separated by the delimiter
      while (std::getline(lineStream, cellString, delim))
      {
        values.push_back(geosx::stringutilities::fromString<typename
           ARRAY::value_type>(cellString));
      }
    }
    inputStream.close();
   }
   else
   {
    GEOS_ERROR("dlmreadVector: Failed to load file:" + filename + " \n");
   }
   }


   ///////////////////////////////////////////////////
   //
   ///// Read from a character deliminated file into a vector of vectors
   template<class ARRAY>
   void dlmreadArray(const std::string& filename, ARRAY& values, char delim = '
      ', int skipLines = 0)
   {

   std::ifstream inputStream(filename.c_str());
   if (inputStream)
   {
    std::string lineString;

    for (int i = 0; i < skipLines; ++i)
      std::getline(inputStream, lineString);

    // read every line from the stream
    while (std::getline(inputStream, lineString))
    {
      std::istringstream lineStream(lineString);
      std::vector < typename ARRAY::value_type > row;
      std::string cellString;
      // read every element from the line separated by the delimiter
      while (std::getline(lineStream, cellString, delim))
      {
        row.push_back(geosx::stringutilities::fromString<typename
           ARRAY::value_type>(cellString));
      }
      values.push_back(row);
    }
    inputStream.close();
   }
   else
   {
    GEOS_ERROR("dlmreadArray: Failed to load file:" + filename + " \n");
   }
   }

   ///// Read from a character deliminated file into a vector of vectors and
      return the transpose
   ///// TODO: Could be improved - not particularly fast or memory efficient
   template<class ARRAY>
   void dlmreadArrayTranspose(const std::string& filename, ARRAY& values, char
      delim = ' ', int skipLines = 0)
   {

   std::ifstream inputStream(filename.c_str());
   std::string lineString;

   if (inputStream)
   {
    for (int i = 0; i < skipLines; ++i)
      std::getline(inputStream, lineString);

    std::vector < std::vector<typename ARRAY::value_type> > valuesT;

    // read every line from the stream
    while (std::getline(inputStream, lineString))
    {
      std::istringstream lineStream(lineString);
      std::vector < typename ARRAY::value_type > row;
      std::string cellString;
      // read every element from the line separated by the delimiter
      while (std::getline(lineStream, cellString, delim))
      {
        row.push_back(geosx::stringutilities::fromString<typename
           ARRAY::value_type>(cellString));
      }
      valuesT.push_back(row);
    }
    inputStream.close();

    localIndex n = valuesT.size();
    localIndex m = valuesT[0].size();
    values.resize(m);
    for (localIndex i = 0; i < m; ++i)
    {
      values[i].resize(n);
      for (localIndex j = 0; j < n; ++j)
        values[i][j] = valuesT[j][i];
    }
   }
   else
   {
    GEOS_ERROR("dlmreadArrayTranspose: Failed to load file:" + filename + "
       \n");
   }
   }
 */

#endif /*IOUTILITIES_H_*/
