/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file EQ36Database.cpp
 */

#include "EQ36Database.hpp"

using namespace std;

namespace geosx
{
using namespace stringutilities;

namespace constitutive
{
EQ36Database::EQ36Database( path const & databaseFileName, string_array const primarySpeciesNames, string_array const secondarySpeciesNames, real64 const temperature ):
  ThermoDatabaseBase( databaseFileName )
{
  CreateChemicalSystem( primarySpeciesNames, secondarySpeciesNames, temperature );
}

void EQ36Database::CreateChemicalSystem( string_array const & primarySpeciesNames, 
                                         string_array const & secondarySpeciesNames )
{
// Read activity coefficient model related enteries
  std::ifstream m_EQ36File ( m_databaseFileName );
  constexpr std::streamsize m_bufferSize = 256;
  char fileLine[m_bufferSize];
  string nextBlockString

  while( m_EQ36File.getline( fileLine, m_bufferSize ))
  {
    std::string fileLineString ( fileLine );
    nextBlockString = "Temperature grid"
    auto found = fileLineString.find( nextBlockString );
    if( found!=std::string::npos )
      break;
  }

  // read and store temperature information
  nextBlockString = "Pressure grid"
  readActivityCoefficientBlock( nextBlockString, m_actCoeffParameters.temperatures ) 
 
  // read and ignore pressure and pressure envelope information
  // As of now, we don't make use of pressure and hence I am not storing it
  nextBlockString = "Pressure envelope"
  readAndIgnoreActivityCoefficientBlock( nextBlockString ) 
  nextBlockString = "Debye-Huckel A_gamma"
  readAndIgnoreActivityCoefficientBlock( nextBlockString ) 
 
  // read and store Debye-Huckel A information
  nextBlockString = "Debye-Huckel A_H"
  readActivityCoefficientBlock( nextBlockString, m_actCoeffParameters.DebyeHuckelAs ) 
 
  // read and ignore Debye-Huckel A_H information
  nextBlockString = "Debye-Huckel B_gamma"
  readAndIgnoreActivityCoefficientBlock( nextBlockString ) 
 
  // read and store Debye-Huckel B information and convert it in the correct units
  nextBlockString = "Debye-Huckel B_H"
  readActivityCoefficientBlock( nextBlockString, m_actCoeffParameters.DebyeHuckelBs ) 
  m_actCoeffParameters.DebyeHuckelBs = m_actCoeffParameters.DebyeHuckelBs*1e8
 
  // read and ignore Debye-Huckel B_H information
  nextBlockString = "B-dot"
  readAndIgnoreActivityCoefficientBlock( nextBlockString ) 
 
  // read and store WATEQB-Dot information
  nextBlockString = "B-dot_H"
  readActivityCoefficientBlock( nextBlockString, m_actCoeffParameters.WATEQBDots ) 

  // if the size of the variables to be interpolated are not the same, print an error statement
  GEOSX_ERROR_IF(
    m_actCoeffParameters.temperatures.size() != m_actCoeffParameters.DebyeHuckelAs.size() || m_actCoeffParameters.temperatures.size() != m_actCoeffParameters.DebyeHuckelBs.size() || m_actCoeffParameters.temperatures.size() != m_actCoeffParameters.WATEQBDots.size(),
    "Internal error when reading database" );

  // Piecewise linear interpolation for the parameter values
  int loopIndex = 0
  while ( temperature < m_actCoeffParameters.temperatures[loopIndex] )
    loopIndex++
  if ( loopIndex < m_actCoeffParameters.temperatures.size() & loopIndex > 0 )
  {
    real64 DebyeHuckelA = m_actCoeffParameters.DebyeHuckelAs[loopIndex-1] +  
                          (m_actCoeffParameters.DebyeHuckelAs[loopIndex] - m_actCoeffParameters.DebyeHuckelAs[loopIndex-1]) / 
                          (m_actCoeffParameters.temperatures[loopIndex] - m_actCoeffParameters.temperatures[loopIndex-1]) * 
                          (temperature - m_actCoeffParameters.temperatures[loopIndex-1])
    real64 DebyeHuckelB = m_actCoeffParameters.DebyeHuckelBs[loopIndex-1] +  
                          (m_actCoeffParameters.DebyeHuckelBs[loopIndex] - m_actCoeffParameters.DebyeHuckelBs[loopIndex-1]) / 
                          (m_actCoeffParameters.temperatures[loopIndex] - m_actCoeffParameters.temperatures[loopIndex-1]) * 
                          (temperature - m_actCoeffParameters.temperatures[loopIndex-1])
    real64 WATEQBDot = m_actCoeffParameters.WATEQBDots[loopIndex-1] +  
                          (m_actCoeffParameters.WATEQBDots[loopIndex] - m_actCoeffParameters.WATEQBDots[loopIndex-1]) / 
                          (m_actCoeffParameters.temperatures[loopIndex] - m_actCoeffParameters.temperatures[loopIndex-1]) * 
                          (temperature - m_actCoeffParameters.temperatures[loopIndex-1])
  }
  GEOSX_ERROR_IF( loopIndex == 0 || loopIndex == m_actCoeffParameters.temperatures.size(), "Input temperature not in permitted range." )

  // Read primary species related entries

  // read and ignore everything till you reach the basis species block
  readAndIgnoreActivityCoefficientBlock( EQ36File, buffer_size, "basis species") 

  // Resize variables that will store the read information
  integer numPrimarySpecies = primarySpeciesNames.size();
  m_MWPrimary.resize( numPrimarySpecies )
  m_ionSizePrimary.resize( numPrimarySpecies )
  m_chargePrimary.resize( numPrimarySpecies )

  // read and store primary species information
  readPrimarySpeciesBlock( EQ36File, buffer_size, "auxiliary basis species", primarySpeciesNames, m_MWPrimary, m_ionSizePrimary, m_chargePrimary )

  // Read secondary species related entries
  // flag to determine is user has provided a list of secondary species to track
  bool inputSecondarySpeciesFlag = secondarySpeciesNames.size() == 0 ? 1 : 0; 

  // We probably don't need this. Will evaluate in next round. 
  unordered_map< string, int > secondarySpeciesMap;
  for( int ic = 0; ic < secondarySpeciesNames.size(); ic++ )
    secondarySpeciesMap[ secondarySpeciesNames[ic] ] = ic;

  // Have to figure out how to size the different variables as we may not know how many secondary species would we need. 
  // read auxiliary basis species information
  readSecondarySpeciesBlock( EQ36File, buffer_size, "aqueous species", secondarySpeciesNames, m_MWSec, m_ionSizeSec, m_chargeSec, m_stoichMatrix, m_log10EqConstTempFunction )

  // read auxiliary aqueous species information
  readSecondarySpeciesBlock( EQ36File, buffer_size, "solids", secondarySpeciesNames, m_MWSec, m_ionSizeSec, m_chargeSec, m_stoichMatrix, m_log10EqConstTempFunction )

  // read solid species information
  readSecondarySpeciesBlock( EQ36File, buffer_size, "liquids", secondarySpeciesNames, m_MWSec, m_ionSizeSec, m_chargeSec, m_stoichMatrix, m_log10EqConstTempFunction )

  // read liquid species information
  readSecondarySpeciesBlock( EQ36File, buffer_size, "gases", secondarySpeciesNames, m_MWSec, m_ionSizeSec, m_chargeSec, m_stoichMatrix, m_log10EqConstTempFunction )

  // read gases species information
  readSecondarySpeciesBlock( EQ36File, buffer_size, "solid solutions", secondarySpeciesNames, m_MWSec, m_ionSizeSec, m_chargeSec, m_stoichMatrix, m_log10EqConstTempFunction )

  is.close();
}


void EQ36Database::readActivityCoefficientBlock( std::ifstream const EQ36File, 
                                                 std::streamsize const bufferSize,
                                                 string const nextBlockString,
                                                 arraySlice1d< real64 > const & readVariable ) const

{

  // This function reads a block from the EQ3/6 Database file related to the B-Dot activity coefficient model and returns the read entries in a vector.

  char fileLine[bufferSize];
  // I am not sure if information on which line was read is preserved when passing information through a function. 
  // That is critical because if the file is going to be read from the first line every single time then this won't work
  while( EQ36File.getline( fileLine, bufferSize ))         //read the line
  {
    std::string fileLineString ( fileLine );
    auto found = fileLineString.find( nextBlockString );   // can you find the next block string string in the line
    if( found!=std::string::npos )                         // if yes, you have reached the next block and should exit
      break;

    string_array lineEntries = Tokenize( fileLineString, " " );      // separate the different numeric enteries and store them
    for( int i = 0; i < lineEntries.size(); i++ )
    {
      readVariable.emplace_back( std::stod( lineEntries[i] ));       // do we need to size the variable beforehand?
    }
  }

}

void EQ36Database::readAndIgnoreActivityCoefficientBlock( std::ifstream const EQ36File, 
                                                          std::streamsize const bufferSize,
                                                          string const nextBlockString ) const

{

  // This function reads a block from the EQ3/6 Database file and ignores it as it isn't currently useful
  char fileLine[bufferSize];
  while( EQ36File.getline( fileLine, bufferSize ))         //read the line
  {
    std::string fileLineString ( fileLine );
    auto found = fileLineString.find( nextBlockString );   // can you find the next block string string in the line
    if( found!=std::string::npos )                         // if yes, you have reached the next block and should exit
      break;

  }

}

void EQ36Database::readPrimarySpeciesBlock( std::ifstream const EQ36File, 
                                            std::streamsize const bufferSize,
                                            string const nextBlockString,
                                            string_array const & primarySpeciesNames,
                                            arraySlice1d< real64 > const & MWPrimary,
                                            arraySlice1d< real64 > const & ionSizePrimary,
                                            arraySlice1d< integer > const & chargePrimary ) const
{

  // We probably don't need this. Left it as is for now, will evaluate later
  unordered_map< string, int > primarySpeciesMap;
  // populate primarySpeciesMap with the input primary species
  for( int ic = 0; ic < primarySpeciesNames.size(); ic++ )
  {
    primarySpeciesMap[primarySpeciesNames[ic] ] = -1;
  }
  int H2OIndex = -1, O2gIndex = -1;

  // define variables to store things
  real64 MW
  real64 ionSize
  integer charge
  string speciesName;

  int count = 0;
  while( EQ36File.getline( fileLine, bufferSize ))          // read the line
  {
    std::string fileLineString( fileLine );
    {
      auto found = fileLineString.find( "+-------" );       // if one finds "+-----" in the line, it means that the next line is either the start of a new block or new primary species
      if( found!=std::string::npos )                        // Either you have information for the last species or you are just starting
      {
        // auto it = primarySpeciesMap.find( speciesName );
        // if( it != primarySpeciesMap.end() || speciesName == "H2O" ||speciesName == "O2(g)" )
        // Only include the primary species mentioned in the input file. 
        // For the time being ignore H2O and O2(g) because I don't know why they need special treatment
        // Not a 100 % sure if the map variable can be replaced by the string array variable while still using the find() and end() functions. 
        auto it = primarySpeciesNames.find( speciesName );          
        if( it != primarySpeciesNames.end() )
        {
          // Ignoring speciesType information since it doesn't make a difference at this point
          MWPrimary[speciesIndex] = MW;
          ionSizePrimary[speciesIndex] = ionSize;
          chargePrimary[speciesIndex] = charge;

          // Have to figure out how are we going to find out the species index to replace this line. Ask Matteo
          primarySpeciesMap[speciesName] = count;
          count++;
        }

        EQ36File.getline( fileLine, buffer_size );          // read next line
        speciesName = fileLine;

        auto found2 = speciesName.find( nextBlockString );  // if the next block has started, exit the loop
        if( found2 != std::string::npos )
        {
          break;
        }
      }
    }

    // molecular weight section
    {
      auto found = fileLineString.find( "mol.wt." );
      if( found != std::string::npos )
      {
        string_array lineEntries = Tokenize( fileLineString, " " );
        MW = std::stod( lineEntries[3] ) * 0.001;
      }
    }

    // ion size section
    {
      auto found = fileLineString.find( "DHazero" );
      if( found != std::string::npos )
      {
        string_array lineEntries = Tokenize( fileLineString, " " );
        ionSize = std::stod( lineEntries[3] ) * 1e-8;
      }
    }
    // charge section
    {
      auto found = fileLineString.find( "charge" );
      if( found != std::string::npos )
      {
        string_array lineEntries = Tokenize( fileLineString, " " );
        charge = std::stoi( lineEntries[2] );
      }
    }
  }

  // Likely we don't need any of this. Will evaluate in the next round. 
  integer idx;
  integer numPrimarySpecies = primarySpeciesNames.size();
  m_primarySpeciesIndices.resize( numPrimarySpecies );
  for( int ic = 0; ic < numPrimarySpecies; ic++ )
  {
    idx = primarySpeciesMap[primarySpeciesNames[ic] ];
    m_primarySpeciesIndices[ic] = idx;
    primarySpeciesMap[primarySpeciesNames[ic] ] = int(ic);
  }
}

void EQ36Database::readSecondarySpeciesBlock( std::ifstream const EQ36File, 
                                              std::streamsize const bufferSize,
                                              string const nextBlockString,
                                              string_array const & secondarySpeciesNames,
                                              arraySlice1d< real64 > const & MWSec,
                                              arraySlice1d< real64 > const & ionSizeSec,
                                              arraySlice1d< integer > const & chargeSec, 
                                              arraySlice2d< real64 > const & stoichMatrix, 
                                              arraySlice1d< real64 > const & log10EqConstTempFunction ) const
{

  // variables to store read information
  string_array speciesNames;
  array1d< integer > speciesIndices;
  array1d< real64 > stoichs;
  array1d< real64 > logKs;

  count = 0;

  while( EQ36File.getline( fileLine, bufferSize ))          // read from file
  {
    std::string fileLineString( fileLine );
    {
      auto found = fileLineString.find( "+-------" );       // If +----- is found, it means you are at the beginning of the next species or the beginning of a new block
      if( found!=std::string::npos )
      {
        // record these entries only when all the products in the reaction associated with this species are listed as a primary species in the input file
        // Have to update how that is determined in the bottom. Will do in next round
        if( primarySpeciesIndices.size() > 1 )
        {
          // Find if the species is listed as a secondary species or if no secondary species are provided in the input file
          // Have to check if find() and end() can be used with the varaible secondarySpeciesNames
          auto iter = secondarySpeciesNames.find( speciesName );
          if( iter != secondarySpeciesNames.end() || inputSecondarySpeciesFlag )
          {
            // Have to figure out how to get the speciesIndex
            if( inputSecondarySpeciesFlag )
            {
              secondarySpeciesNames.emplace_back( speciesName )
            }
            // Ignoring speciesType information since I don't know what difference does it make
            MWSec[speciesIndex] = MW;
            ionSizeSec[speciesIndex] = ionSize;
            chargeSec[speciesIndex] = charge;
            log10EqConstTempFunction[speciesIndex][:] = logKs;         // rectify syntax
            stoichMatrix[speciesIndex][:] = stoichs;                   // rectify syntax

            // This will probably get replaced by whatever lines we use to determine the species index
            count++;
          }

          speciesIndices.clear();
          stoichs.clear();
          logKs.clear();
          speciesNames.clear();
        }

        EQ36File.getline( fileLine, bufferSize );                       // read the next line
        std::string fileLineString2( fileLine );
        string_array lineEntries2 = Tokenize( fileLineString2, " " );
        speciesName = lineEntries2[0];

        auto found2 = lineEntries2.find( "aqueous species" );           // if reached next block, exit the loop
        if( found2 != std::string::npos )
        {
          break;
        }
      }
    }

    // Molecular weight section
    {
      auto found = fileLineString.find( "mol.wt." );
      if( found != std::string::npos )
      {
        string_array lineEntries = Tokenize( fileLineString, " " );
        MW = std::stod( lineEntries[3] ) * 0.001;
      }
    }

    // Ion size section
    {
      auto found = fileLineString.find( "DHazero" );
      if( found != std::string::npos )
      {
        string_array lineEntries = Tokenize( fileLineString, " " );
        ionSize = std::stod( lineEntries[3] ) * 1e-8;
      }
    }

    // charge section
    {
      auto found = fileLineString.find( "charge" );
      if( found != std::string::npos )
      {
        string_array lineEntries = Tokenize( fileLineString, " " );
        charge = std::stoi( lineEntries[2] );
      }
    }

    // stoichiometric section
    {
      auto found = fileLineString.find( "aqueous dissociation reaction" );
      if( found != std::string::npos )
      {
        string_array lineEntries = Tokenize( fileLineString, " " );
        integer numSpecies =  std::stoi( lineEntries[0] );       // number of species that participate in the reaction
        speciesNames.clear();
        stoichs.clear();

        while( EQ36File.getline( fileLine, bufferSize ))
        {
          std::string fileLineString2( fileLine );
          auto found2 = fileLineString2.find( "* Log K" );      // reached the end of this block so exit loop
          if( found2 != std::string::npos )
            break;

          string_array lineEntries2 = Tokenize( fileLineString2, " " );
          for( int i = 0; i < lineEntries2.size(); i++ )
          {
            if( i % 2 == 0 )
              stoichs.emplace_back( std::stod( lineEntries2[i] ));
            else
              // have to figure out a way to find the index corresponding to this primary species. Needs to be done in the next round
              // Will have to treat H2O as special case
              speciesNames.emplace_back( lineEntries2[i] );
          }
        }

        GEOSX_ERROR_IF( numSpecies != speciesNames.size() || numSpecies != stoichs.size() || speciesName != speciesNames[0], "Internal error when reading database" );

        // Find if all the  species listed in the dissociation reaction are listed as primary species in the input file 
        // Will have to write something special for H2O
        // This might get combined in the piece of code we need to write to tie the primary species in the reaction to its corresponding index
        bool notFound = 0;
        speciesIndices.resize( numSpecies );
        speciesIndices[0] = count;

        for( int i = 1; i < numSpecies; i++ )
        {
          auto it = primarySpeciesNames.find( speciesNames[i] );
          if( it != primarySpeciesNames.end())
          {
            speciesIndices[i] = it->second;
          }
          else
          {
            notFound = 1;
            break;
          }
        }

        // If even one of the species in the reaction is not listed as a primary species in the input file, do not include this secondary species
        // Else, read the equilibrium constant section.         
        // This may also get combined in the piece of code we need to write to tie the primary species in the reaction to its corresponding index
        if( notFound )
        {
          speciesIndices.clear();
        }
        else
        {
          logKs.clear();
          while( EQ36File.getline( fileLine, bufferSize ))
          {
            std::string fileLineString2( fileLine );
            auto found2 = fileLineString2.find( "*" );          // if reached next section, exit the loop
            if( found2 != std::string::npos )
              break;

            string_array lineEntries2 = Tokenize( fileLineString2, " " );
            for( int i = 0; i < lineEntries2.size(); i++ )
              logKs.emplace_back( std::stod( lineEntries2[i] ));
          }
        }
      }
    }
  }
}

// Not sure what is this for
REGISTER_CATALOG_ENTRY( ThermoDatabaseBase,
                        EQ36Database,
                        const Path &,
                        const string_array &,
                        const string_array & )

}

}
