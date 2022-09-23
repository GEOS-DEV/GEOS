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
  nextBlockString = "basis species"
  readAndIgnoreActivityCoefficientBlock( nextBlockString ) 

  // Resize variables that will store the read information
  integer numPrimarySpecies = primarySpeciesNames.size();
  MWPrimary.resize( numPrimarySpecies )
  ionSizePrimary.resize( numPrimarySpecies )
  chargePrimary.resize( numPrimarySpecies )

  // read and store primary species information
  nextBlockString = "auxiliary basis species"
  readPrimarySpeciesBlock( nextBlockString, primarySpeciesNames, MWPrimary, ionSizePrimary, chargePrimary )

  // Read secondary species related entries
  // flag to determine is user has provided a list of secondary species to track
  // Code is currently not set up to work when user does not list secondary species due to difficulty in working with unknown sizes
  // To allow for automatic population of secondary species based on the database a function similar to readSeconarySpeciesBlock would 
  // have to be written that just counts the number of secondary species and creates the secondarySpeciesNames variable
  bool inputSecondarySpeciesFlag = secondarySpeciesNames.size() == 0 ? 1 : 0; 
  GEOSX_ERROR_IF( inputSecondarySpeciesFlag, "Currently user is expected to supply secondary species names as well." )

  // create a map to get the index numbers of the secondary species
  unordered_map< string, int > secondarySpeciesMap;
  for( int ic = 0; ic < secondarySpeciesNames.size(); ic++ )
    secondarySpeciesMap[ secondarySpeciesNames[ic] ] = ic;

  // read auxiliary basis species information
  nextBlockString = "aqueous species"
  readSecondarySpeciesBlock( nextBlockString, secondarySpeciesNames, MWSec, ionSizeSec, chargeSec, stoichMatrix, log10EqConstTempFunction )

  // read auxiliary aqueous species information
  nextBlockString = "solids"
  readSecondarySpeciesBlock( nextBlockString, secondarySpeciesNames, MWSec, ionSizeSec, chargeSec, stoichMatrix, log10EqConstTempFunction )

  // read solid species information
  nextBlockString = "liquids"
  readSecondarySpeciesBlock( nextBlockString, secondarySpeciesNames, MWSec, ionSizeSec, chargeSec, stoichMatrix, log10EqConstTempFunction )

  // read liquid species information
  nextBlockString = "gases"
  readSecondarySpeciesBlock( nextBlockString, secondarySpeciesNames, MWSec, ionSizeSec, chargeSec, stoichMatrix, log10EqConstTempFunction )

  // read gases species information
  nextBlockString = "solid solutions"
  readSecondarySpeciesBlock( nextBlockString, secondarySpeciesNames, MWSec, ionSizeSec, chargeSec, stoichMatrix, log10EqConstTempFunction )

  is.close();
}


void EQ36Database::readActivityCoefficientBlock( string const nextBlockString,
                                                 arraySlice1d< real64 > const & readVariable ) const

{

  // This function reads a block from the EQ3/6 Database file related to the B-Dot activity coefficient model and returns the read entries.

  char fileLine[m_bufferSize];
  while( m_EQ36File.getline( fileLine, m_bufferSize ))         //read the line
  {
    std::string fileLineString ( fileLine );
    auto found = fileLineString.find( nextBlockString );   // can you find the next block string string in the line
    if( found!=std::string::npos )                         // if yes, you have reached the next block and should exit
      break;

    string_array lineEntries = Tokenize( fileLineString, " " );      // separate the different numeric enteries and store them
    for( int i = 0; i < lineEntries.size(); i++ )
    {
      // Ask Matteo if this is acceptable as readVariable hasn't been resized anywhere
      readVariable.emplace_back( std::stod( lineEntries[i] ));       // do we need to size the variable beforehand?
    }
  }

}

void EQ36Database::readAndIgnoreActivityCoefficientBlock( string const nextBlockString ) const

{

  // This function reads a block from the EQ3/6 Database file and ignores it as it isn't currently useful
  char fileLine[m_bufferSize];
  while( m_EQ36File.getline( fileLine, m_bufferSize ))         //read the line
  {
    std::string fileLineString ( fileLine );
    auto found = fileLineString.find( nextBlockString );   // can you find the next block string string in the line
    if( found!=std::string::npos )                         // if yes, you have reached the next block and should exit
      break;

  }

}

void EQ36Database::readPrimarySpeciesBlock( string const nextBlockString,
                                            string_array const & primarySpeciesNames,
                                            arraySlice1d< real64 > const & MWPrimary,
                                            arraySlice1d< real64 > const & ionSizePrimary,
                                            arraySlice1d< integer > const & chargePrimary ) const
{

  // create a map to get the index numbers of the primary species
  unordered_map< string, int > primarySpeciesMap;
  for( int ic = 0; ic < primarySpeciesNames.size(); ic++ )
  {
    primarySpeciesMap[primarySpeciesNames[ic] ] = ic;
  }

  // define variables to store things
  real64 MW
  real64 ionSize
  integer charge
  string speciesName = " ";

  int countPrimary = 0;
  while( m_EQ36File.getline( fileLine, m_bufferSize ))          // read the line
  {
    std::string fileLineString( fileLine );
    {
      auto found = fileLineString.find( "+-------" );       // if one finds "+-----" in the line, it means that the next line is either the start of a new block or new primary species
      if( found!=std::string::npos )                        // Either you have all the information for the previous species or you are just starting
      {
        auto it = primarySpeciesMap.find( speciesName );          
        if( it != primarySpeciesMap.end() )
        {
          // Ignoring speciesType information since it doesn't make a difference at this point
          int speciesIndex = primarySpeciesMap[ speciesName ] 
          MWPrimary[speciesIndex] = MW;
          ionSizePrimary[speciesIndex] = ionSize;
          chargePrimary[speciesIndex] = charge;
          countPrimary++
        }

        m_EQ36File.getline( fileLine, m_bufferSize );          // read next line
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

  GEOSX_ERROR_IF( countPrimary != primarySpeciesNames.size(), "Could not find all input primary species." )
}

void EQ36Database::readSecondarySpeciesBlock( string const nextBlockString,
                                              unordered_map const & secondarySpeciesMap,
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

  countSecondary = 0;
  bool primarySpeciesFound = 0

  while( m_EQ36File.getline( fileLine, m_bufferSize ))          // read from file
  {
    std::string fileLineString( fileLine );
    {
      auto found = fileLineString.find( "+-------" );       // If +----- is found, it means you are at the beginning of the next species or the beginning of a new block
      if( found!=std::string::npos )
      {
        // record these entries only when all the products in the reaction associated with this species are listed as a primary species in the input file
        // Have to update how that is determined in the bottom. Will do in next round
        if( primarySpeciesFound )
        {
          // Find if the species is listed as a secondary species 
          auto iter = secondarySpeciesMap.find( speciesName );
          if( iter != secondarySpeciesMap.end() )
          {

            int secSpeciesIndex = secondarySpeciesMap[speciesName]

            // Ignoring speciesType information since I don't know what difference does it make
            MWSec[secSpeciesIndex] = MW;
            ionSizeSec[secSpeciesIndex] = ionSize;
            chargeSec[secSpeciesIndex] = charge;
            log10EqConstTempFunction[secSpeciesIndex][:] = logKs;         // rectify syntax
            stoichMatrix[secSpeciesIndex][:] = stoichs;                   // rectify syntax

            countSecondary++;
          }

          speciesIndices.clear();
          stoichs.clear();
          logKs.clear();
          speciesNames.clear();
        }

        m_EQ36File.getline( fileLine, m_bufferSize );                       // read the next line
        std::string fileLineString2( fileLine );
        string_array lineEntries2 = Tokenize( fileLineString2, " " );
        speciesName = lineEntries2[0];

        auto found2 = lineEntries2.find( nextBlockString );           // if reached next block, exit the loop
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
        integer numRxnSpecies =  std::stoi( lineEntries[0] );       // number of species that participate in the reaction
        reactionPrimarySpeciesNames.clear();    
        reactionPrimarySpeciesIndices.clear();    
        stoichs.clear();
        reactionPrimarySpeciesNames.resize( numRxnSpecies-1 )
        reactionPrimarySpeciesIndices.resize( numRxnSpecies-1 )
        stoichs.resize( numRxnSpecies-1 )

        int reactionSpeciesCount = 0
        while( m_EQ36File.getline( fileLine, m_bufferSize ))
        {
          std::string fileLineString2( fileLine );
          auto found2 = fileLineString2.find( "* Log K" );      // reached the end of this block so exit loop
          if( found2 != std::string::npos )
            break;

          string_array lineEntries2 = Tokenize( fileLineString2, " " );
          for( int i = 0; i < lineEntries2.size(); i = i+2 )
          {
            GEOSX_ERROR_IF( reactionSpeciesCount == 0 & lineEntries2[0] != speciesName, "Internal error when reading database" );
            GEOSX_ERROR_IF( reactionSpeciesCount == 0 & lineEntries2[1] != -1, "First species in the reaction is expected to have a stoichiometry of -1" );
            reactionSpeciesCount++
            if( reactionSpeciesCount > 1 )
            {
              auto it = primarySpeciesMap.find( lineEntries2[i+1] );
              // Note that no special treatment for H2O or O2 gas has been implemented. 
              // If these are present in one's problem, they need to be listed as input primary species. 
              if ( it != primarySpeciesMap.end() )
              {
                reactionPrimarySpeciesNames[reactionSpeciesCount-1] = lineEntries2[i+1];
                stoichs[reactionSpeciesCount-1] = lineEntries2[i];
                reactionPrimarySpeciesIndices[reactionSpeciesCount-1] = primarySpeciesMap[ lineEntries2[i+1] ];
                primarySpeciesFound = 1
              }
              else
              {
                primarySpeciesFound = 0
                break
              }
            }
          }
        }

        GEOSX_ERROR_IF( reactionSpeciesCount != numRxnSpecies, "Internal error when reading database" );

        // If even one of the species in the reaction is not listed as a primary species in the input file, do not include this secondary species
        // Else, read the equilibrium constant section.         
        // This may also get combined in the piece of code we need to write to tie the primary species in the reaction to its corresponding index
        if( primarySpeciesFound )
        {
          logKs.clear();
          while( m_EQ36File.getline( fileLine, m_bufferSize ))
          {
            std::string fileLineString2( fileLine );
            auto found2 = fileLineString2.find( "*" );          // if reached next section, exit the loop
            if( found2 != std::string::npos )
              break;

            string_array lineEntries2 = Tokenize( fileLineString2, " " );
            for( int i = 0; i < lineEntries2.size(); i++ )
              logKs.emplace_back( std::stod( lineEntries2[i] ));      // Check with Matteo what to do about this since we don't know the size of the vector
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
