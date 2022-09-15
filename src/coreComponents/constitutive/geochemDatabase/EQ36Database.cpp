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
EQ36Database::EQ36Database( path const databaseFileName, string_array const primarySpeciesNames, string_array const secondarySpeciesNames ):
  ThermoDatabaseBase( databaseFileName )
{
  CreateChemicalSystem( primarySpeciesNames, secondarySpeciesNames );
}

void EQ36Database::CreateChemicalSystem( string_array const primarySpeciesNames, 
                                         string_array const secondarySpeciesNames )
{
// Read activity coefficient model related enteries
  std::ifstream is( m_databaseFileName );
  constexpr std::streamsize buf_size = 256;
  char buf[buf_size];

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    auto found = str.find( "Temperature grid" );
    if( found!=std::string::npos )
      break;
  }

// Lines 54 - 67 repeat multiple times with the only difference being in 
// lines 57 and 65 (the term it is searching for to determine the end of 
// the block and the variable in which to store the enteries (including unit conversion))  
  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    auto found = str.find( "Pressure grid" );
    if( found!=std::string::npos )
      break;

    string_array strs = Tokenize( str, " " );
    for( int i = 0; i < strs.size(); i++ )
    {
      // store the temperature grid
      m_actCoeffParameters.temperatures.emplace_back( std::stod( strs[i] ));
    }
  }

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    auto found = str.find( "Pressure envelope" );
    if( found!=std::string::npos )
      break;

    string_array strs = Tokenize( str, " " );
    for( int i = 0; i < strs.size(); i++ )
    {
      // store the pressure grid
      // I think this is the saturation pressure of water at the temperatures in the temperature grid
      m_actCoeffParameters.pressures.emplace_back( std::stod( strs[i] ));
    }
  }

// Lines 93 - 99 repeat multiple times with the only difference being in lines 96
// (the term it is searching for to determine the end of the block)  
// This block is used when the database enteries are not necessary
  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    auto found = str.find( "Debye-Huckel A_gamma" );
    if( found!=std::string::npos )
      break;
  }

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    auto found = str.find( "Debye-Huckel A_H" );
    if( found!=std::string::npos )
      break;

    string_array strs = Tokenize( str, " " );
    for( int i = 0; i < strs.size(); i++ )
    {
      m_actCoeffParameters.DebyeHuckelAs.emplace_back( std::stod( strs[i] ));
    }
  }

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    auto found = str.find( "Debye-Huckel B_gamma" );
    if( found!=std::string::npos )
      break;
  }

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    auto found = str.find( "Debye-Huckel B_H" );
    if( found!=std::string::npos )
      break;

    string_array strs = Tokenize( str, " " );
    for( int i = 0; i < strs.size(); i++ )
    {
      m_actCoefParameters.DebyeHuckelBs.emplace_back( std::stod( strs[i] )*1e8 );
    }
  }

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    auto found = str.find( "B-dot" );
    if( found!=std::string::npos )
      break;
  }

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    auto found = str.find( "B-dot_H" );
    if( found!=std::string::npos )
      break;

    string_array strs = Tokenize( str, " " );
    for( int i = 0; i < strs.size(); i++ )
    {
      m_actCoeffParameters.WATEQBDots.emplace_back( std::stod( strs[i] ));
    }
  }

  GEOSX_ERROR_IF(
    m_actCoeffParameters.temperatures.size() != m_actCoeffParameters.pressures.size() || m_actCoeffParameters.temperatures.size() != m_actCoeffParameters.DebyeHuckelAs.size() || m_actCoeffParameters.temperatures.size() != m_actCoeffParameters.DebyeHuckelBs.size() || m_actCoeffParameters.temperatures.size() != m_actCoeffParameters.WATEQBDots.size(),
    "Internal error when reading database" );


// Read primary species related enteries

  unordered_map< string, int > primarySpeciesMap;
// populate primarySpeciesMap with the input primary species
  for( int ic = 0; ic < primarySpeciesNames.size(); ic++ )
  {
    primarySpeciesMap[primarySpeciesNames[ic] ] = -1;
  }

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    auto found = str.find( "basis species" );
    if( found!=std::string::npos )
      break;
  }

  string speciesName;
  real64 MW = 0;
  integer charge = 0;
  real64 ionSize = 0;
  int H2OIndex = -1, O2gIndex = -1;

  int count = 0;
  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    {
      auto found = str.find( "+-------" );
      if( found!=std::string::npos )
      {
        auto it = primarySpeciesMap.find( speciesName );
        // Only include the primary species mentioned in the input file, H2O and O2(g)
        if( it != primarySpeciesMap.end() || speciesName == "H2O" ||speciesName == "O2(g)" )
        {
          speciesProperties entry;
          entry.name = speciesName;
          entry.type = speciesType::aqueous;
          entry.MW = MW;
          entry.ionSize = ionSize;
          entry.charge = charge;
          m_primarySpecies.emplace_back( entry );

          if( speciesName == "H2O" )
            H2OIndex = count;
          else if( speciesName == "O2(g)" )
            O2gIndex = count;
          else
            primarySpeciesMap[speciesName] = count;

          count++;
        }
        is.getline( buf, buf_size );
        speciesName = buf;
        auto found2 = speciesName.find( "auxiliary basis species" );

        if( found2 != std::string::npos )
        {
          break;
        }
      }
    }

    {
      auto found = str.find( "mol.wt." );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        MW = std::stod( strs[3] ) * 0.001;
      }
    }

    {
      auto found = str.find( "DHazero" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        ionSize = std::stod( strs[3] ) * 1e-8;
      }
    }

    {
      auto found = str.find( "charge" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        charge = std::stoi( strs[2] );
      }
    }
  }

  integer idx;
  integer numPrimarySpecies = primarySpeciesNames.size();
  m_primarySpeciesIndices.resize( numPrimarySpecies + 2 );
  for( int ic = 0; ic < numPrimarySpecies; ic++ )
  {
    idx = primarySpeciesMap[primarySpeciesNames[ic] ];
    m_primarySpeciesIndices[ic] = idx;
    primarySpeciesMap[primarySpeciesNames[ic] ] = int(ic);
  }

  m_primarySpeciesIndices[numPrimarySpecies] = H2OIndex;
  primarySpeciesMap[m_primarySpecies[H2OIndex].name] = int(numPrimarySpecies);

  m_primarySpeciesIndices[numPrimarySpecies + 1] = O2gIndex;
  primarySpeciesMap[m_primarySpecies[O2gIndex].name ] = int(numPrimarySpecies + 1);


// Read secondary species related enteries
  /* read aux basis species */
// Not sure why some species are defined as auxilliary basis species
  unordered_map< string, int > secondarySpeciesMap;
  for( int ic = 0; ic < secondarySpeciesNames.size(); ic++ )
    secondarySpeciesMap[ secondarySpeciesNames[ic] ] = ic;
 // I am guessing this is to indicate whether the user has provided a list of names as secondary species. Not 100 % sure
  bool inputSecondarySpeciesFlag = secondarySpeciesNames.size() == 0 ? 1 : 0; 
  
// Lines 278 - 430 repeat multiple times with the only difference being in lines 297 and 324 
// (the species type and the term it is searching for to determine the end of the block)  
  string_array speciesNames;
  array1d< integer > speciesIndices;
  array1d< real64 > stoichs;
  array1d< real64 > logKs;

  count = 0;

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    {
      auto found = str.find( "+-------" );
      if( found!=std::string::npos )
      {
// speciesIndices has elements only when all the products in the reaction associated with this species are listed as a primary species in the input file
        if( speciesIndices.size() > 1 )
        {
          speciesProperties entry;
          entry.name = speciesName;
          entry.type = speciesType::aqueous;
          entry.MW = MW;
          entry.ionSize = ionSize;
          entry.charge = charge;
          entry.speciesIndices = speciesIndices;
          entry.stoichs = stoichs;
          entry.logKs = logKs;

// Find if the auxiliary species is listed as a secondary species or if no secondary species are provided in the input file
          auto iter = secondarySpeciesMap.find( speciesName );
          if( iter != secondarySpeciesMap.end() || inputSecondarySpeciesFlag )
          {
            m_secondarySpecies.emplace_back( entry );
            count++;
          }

          speciesIndices.clear();
          stoichs.clear();
          logKs.clear();
          speciesNames.clear();
        }

        is.getline( buf, buf_size );
        std::string str2( buf );
        string_array strs2 = Tokenize( str2, " " );
        speciesName = strs2[0];

        auto found2 = str2.find( "aqueous species" );
        if( found2 != std::string::npos )
        {
          break;
        }
      }
    }

    {
      auto found = str.find( "mol.wt." );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        MW = std::stod( strs[3] ) * 0.001;
      }
    }

    {
      auto found = str.find( "DHazero" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        ionSize = std::stod( strs[3] ) * 1e-8;
      }
    }

    {
      auto found = str.find( "charge" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        charge = std::stoi( strs[2] );
      }
    }

    {
      auto found = str.find( "aqueous dissociation reaction" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        integer num1 =  std::stoi( strs[0] );
        speciesNames.clear();
        stoichs.clear();

        while( is.getline( buf, buf_size ))
        {
          std::string str2( buf );
          auto found2 = str2.find( "* Log K" );
          if( found2 != std::string::npos )
            break;


          string_array strs2 = Tokenize( str2, " " );
          for( int i = 0; i < strs2.size(); i++ )
          {
            if( i % 2 == 0 )
              stoichs.emplace_back( std::stod( strs2[i] ));
            else
              speciesNames.emplace_back( strs2[i] );
          }
        }

        GEOSX_ERROR_IF( num1 != speciesNames.size() || num1 != stoichs.size() || speciesName != speciesNames[0], "Internal error when reading database" );

// Find if all the  species listed in the dissociation reaction are listed as primary species in the input file 
        bool notFound = 0;
        speciesIndices.resize( num1 );
        speciesIndices[0] = count;

        for( int i = 1; i < num1; i++ )
        {
          auto it = primarySpeciesMap.find( speciesNames[i] );
          if( it != primarySpeciesMap.end())
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
// Else, read the equilibrium constant block.         
        if( notFound )
        {
          speciesIndices.clear();
        }
        else
        {
          logKs.clear();
          while( is.getline( buf, buf_size ))
          {
            std::string str2( buf );
            auto found2 = str2.find( "*" );
            if( found2 != std::string::npos )
              break;

            string_array strs2 = Tokenize( str2, " " );
            for( int i = 0; i < strs2.size(); i++ )
              logKs.emplace_back( std::stod( strs2[i] ));
          }
        }
      }
    }
  }

  /* read aux aqueous species */
  speciesIndices.clear();
  stoichs.clear();
  logKs.clear();
  speciesNames.clear();

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    {
      auto found = str.find( "+-------" );
      if( found!=std::string::npos )
      {
        if( speciesIndices.size() > 1 )
        {
          speciesProperties entry;
          entry.name = speciesName;
          entry.type = speciesType::aqueous;
          entry.MW = MW;
          entry.ionSize = ionSize;
          entry.charge = charge;
          entry.speciesIndices = speciesIndices;
          entry.stoichs = stoichs;
          entry.logKs = logKs;

          auto iter = secondarySpeciesMap.find( speciesName );
          if( iter != secondarySpeciesMap.end() || inputSecondarySpeciesFlag )
          {
            m_secondarySpecies.emplace_back( entry );
            count++;
          }

          speciesIndices.clear();
          stoichs.clear();
          logKs.clear();
          speciesNames.clear();
        }

        is.getline( buf, buf_size );
        std::string str2( buf );
        string_array strs2 = Tokenize( str2, " " );
        speciesName = strs2[0];

        auto found2 = str2.find( "solids" );
        if( found2 != std::string::npos )
        {
          break;
        }
      }
    }

    {
      auto found = str.find( "mol.wt." );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        MW = std::stod( strs[3] ) * 0.001;
      }
    }

    {
      auto found = str.find( "DHazero" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        ionSize = std::stod( strs[3] ) * 1e-8;
      }
    }

    {
      auto found = str.find( "charge" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        charge = std::stoi( strs[2] );
      }
    }

    {
      auto found = str.find( "aqueous dissociation reaction" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        integer num1 = std::stoi( strs[0] );
        speciesNames.clear();
        stoichs.clear();

        while( is.getline( buf, buf_size ))
        {
          std::string str2( buf );
          auto found2 = str2.find( "* Log K" );
          if( found2 != std::string::npos )
            break;

          string_array strs2 = Tokenize( str2, " " );
          for( int i = 0; i < strs2.size(); i++ )
          {
            if( i % 2 == 0 )
              stoichs.emplace_back( std::stod( strs2[i] ));
            else
              speciesNames.emplace_back( strs2[i] );
          }
        }

        GEOSX_ERROR_IF( num1 != speciesNames.size() || num1 != stoichs.size() || speciesName != speciesNames[0], "Internal error when reading database" );

        bool notFound = 0;
        speciesIndices.resize( num1 );
        speciesIndices[0] = count;

        for( int i = 1; i < num1; i++ )
        {
          auto it = primarySpeciesMap.find( speciesNames[i] );
          if( it != primarySpeciesMap.end())
          {
            speciesIndices[i] = it->second;
          }
          else
          {
            notFound = 1;
            break;
          }
        }

        if( notFound )
        {
          speciesIndices.clear();
        }
        else
        {
          logKs.clear();
          while( is.getline( buf, buf_size ))
          {

            std::string str2( buf );
            auto found2 = str2.find( "*" );
            if( found2 != std::string::npos )
              break;

            string_array strs2 = Tokenize( str2, " " );
            for( int i = 0; i < strs2.size(); i++ )
              logKs.emplace_back( std::stod( strs2[i] ));

          }
        }
      }
    }
  }

  /* read solid species */
  speciesIndices.clear();
  stoichs.clear();
  logKs.clear();
  speciesNames.clear();

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    {
      auto found = str.find( "+-------" );
      if( found!=std::string::npos )
      {
        if( speciesIndices.size() > 1 )
        {
          speciesProperties entry;
          entry.name = speciesName;
          entry.type = speciesType::solid;
          entry.MW = MW;
          entry.ionSize = 0;
          entry.charge = 0;
          entry.speciesIndices = speciesIndices;
          entry.stoichs = stoichs;
          entry.logKs = logKs;

          auto iter = secondarySpeciesMap.find( speciesName );
          if( iter != secondarySpeciesMap.end() || inputSecondarySpeciesFlag )
          {
            m_secondarySpecies.emplace_back( entry );
            count++;
          }

          speciesIndices.clear();
          stoichs.clear();
          logKs.clear();
          speciesNames.clear();

        }

        is.getline( buf, buf_size );
        std::string str2( buf );
        string_array strs2 = Tokenize( str2, " " );
        speciesName = strs2[0];

        auto found2 = str2.find( "liquids" );
        if( found2 != std::string::npos )
        {
          break;
        }
      }
    }

    {
      auto found = str.find( "mol.wt." );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        MW = std::stod( strs[3] ) * 0.001;
      }
    }

    {
      auto found = str.find( "DHazero" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        ionSize = std::stod( strs[3] ) * 1e-8;
      }
    }

    {
      auto found = str.find( "charge" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        charge = std::stoi( strs[2] );
      }
    }

    {
      auto found = str.find( "aqueous dissociation reaction" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        integer num1 = std::stoi( strs[0] );

        speciesNames.clear();
        stoichs.clear();

        while( is.getline( buf, buf_size ))
        {
          std::string str2( buf );
          auto found2 = str2.find( "* Log K" );
          if( found2 != std::string::npos )
            break;

          string_array strs2 = Tokenize( str2, " " );
          for( int i = 0; i < strs2.size(); i++ )
          {
            if( i % 2 == 0 )
              stoichs.emplace_back( std::stod( strs2[i] ));
            else
              speciesNames.emplace_back( strs2[i] );
          }
        }

        GEOSX_ERROR_IF( num1 != speciesNames.size() || num1 != stoichs.size() || speciesName != speciesNames[0], "Internal error when reading database" );

        bool notFound = 0;

        speciesIndices.resize( num1 );
        speciesIndices[0] = count;

        for( int i = 1; i < num1; i++ )
        {
          auto it = primarySpeciesMap.find( speciesNames[i] );
          if( it != primarySpeciesMap.end())
          {
            speciesIndices[i] = it->second;
          }
          else
          {
            notFound = 1;
            break;
          }
        }

        if( notFound )
        {
          speciesIndices.clear();
        }
        else
        {
          logKs.clear();
          while( is.getline( buf, buf_size ))
          {
            std::string str2( buf );
            auto found2 = str2.find( "*" );
            if( found2 != std::string::npos )
              break;

            string_array strs2 = Tokenize( str2, " " );
            for( int i = 0; i < strs2.size(); i++ )
              logKs.emplace_back( std::stod( strs2[i] ));

          }
        }
      }
    }
  }


  /* read liquid species */
  speciesIndices.clear();
  stoichs.clear();
  logKs.clear();
  speciesNames.clear();

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    {
      auto found = str.find( "+-------" );
      if( found!=std::string::npos )
      {
        if( speciesIndices.size() > 1 )
        {
          speciesProperties entry;
          entry.name = speciesName;
          entry.type = speciesType::liquid;
          entry.MW = MW;
          entry.ionSize = 0;
          entry.charge = 0;
          entry.speciesIndices = speciesIndices;
          entry.stoichs = stoichs;
          entry.logKs = logKs;

// Not sure why the lines about dependentSpeciesMap is not included here
          m_secondarySpecies.emplace_back( entry );

          speciesIndices.clear();
          stoichs.clear();
          logKs.clear();
          speciesNames.clear();

          count++;
        }

        is.getline( buf, buf_size );
        std::string str2( buf );
        string_array strs2 = Tokenize( str2, " " );
        speciesName = strs2[0];

        auto found2 = str2.find( "gases" );

        if( found2 != std::string::npos )
        {
          break;
        }
      }
    }

    {
      auto found = str.find( "mol.wt." );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        MW = std::stod( strs[3] ) * 0.001;
      }
    }

    {
      auto found = str.find( "DHazero" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        ionSize = std::stod( strs[3] ) * 1e-8;
      }
    }

    {
      auto found = str.find( "charge" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        charge = std::stoi( strs[2] );
      }
    }

    {
      auto found = str.find( "aqueous dissociation reaction" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        integer num1 = std::stoi( strs[0] );

        speciesNames.clear();
        stoichs.clear();

        while( is.getline( buf, buf_size ))
        {
          std::string str2( buf );
          auto found2 = str2.find( "* Log K" );
          if( found2 != std::string::npos )
            break;

          string_array strs2 = Tokenize( str2, " " );
          for( int i = 0; i < strs2.size(); i++ )
          {
            if( i % 2 == 0 )
              stoichs.emplace_back( std::stod( strs2[i] ));
            else
              speciesNames.emplace_back( strs2[i] );
          }
        }

        GEOSX_ERROR_IF( num1 != speciesNames.size() || num1 != stoichs.size() || speciesName != speciesNames[0], "Internal error when reading database" );

        bool notFound = 0;

        speciesIndices.resize( num1 );
        speciesIndices[0] = count;

        for( localIndex i = 1; i < num1; ++i )
        {
          auto it = primarySpeciesMap.find( speciesNames[i] );
          if( it != primarySpeciesMap.end())
          {
            speciesIndices[i] = it->second;
          }
          else
          {
            notFound = 1;
            break;
          }
        }

        if( notFound )
        {
          speciesIndices.clear();
        }
        else
        {
          logKs.clear();
          while( is.getline( buf, buf_size ))
          {
            std::string str2( buf );
            auto found2 = str2.find( "*" );
            if( found2 != std::string::npos )
              break;

            string_array strs2 = Tokenize( str2, " " );
            for( int i = 0; i < strs2.size(); i++ )
              logKs.emplace_back( std::stod( strs2[i] ));
          }
        }
      }
    }
  }

  /* read gases species */
  speciesIndices.clear();
  stoichs.clear();
  logKs.clear();
  speciesNames.clear();

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    {
      auto found = str.find( "+-------" );
      if( found!=std::string::npos )
      {
        if( speciesIndices.size() > 1 )
        {
          speciesProperties entry;
          entry.name = speciesName;
          entry.type = speciesType::gas;
          entry.MW = MW;
          entry.ionSize = 0;
          entry.charge = 0;
          entry.speciesIndices = speciesIndices;
          entry.stoichs = stoichs;
          entry.logKs = logKs;


          auto iter = secondarySpeciesMap.find( speciesName );
          if( iter != secondarySpeciesMap.end() || inputSecondarySpeciesFlag )
          {
            m_secondarySpecies.emplace_back( entry );
            count++;
          }

          speciesIndices.clear();
          stoichs.clear();
          logKs.clear();
          speciesNames.clear();
        }

        is.getline( buf, buf_size );
        std::string str2( buf );
        string_array strs2 = Tokenize( str2, " " );
        speciesName = strs2[0];

        auto found2 = str2.find( "solid solutions" );
        if( found2 != std::string::npos )
        {
          break;
        }
      }
    }

    {
      auto found = str.find( "mol.wt." );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        MW = std::stod( strs[3] ) * 0.001;
      }
    }

    {
      auto found = str.find( "DHazero" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        ionSize = std::stod( strs[3] ) * 1e-8;
      }
    }

    {
      auto found = str.find( "charge" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        charge = std::stoi( strs[2] );
      }
    }

    {
      auto found = str.find( "aqueous dissociation reaction" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        integer num1 = std::stoi( strs[0] );
        speciesNames.clear();
        stoichs.clear();

        while( is.getline( buf, buf_size ))
        {
          std::string str2( buf );
          auto found2 = str2.find( "* Log K" );
          if( found2 != std::string::npos )
            break;

          string_array strs2 = Tokenize( str2, " " );
          for( int i = 0; i < strs2.size(); i++ )
          {
            if( i % 2 == 0 )
              stoichs.emplace_back( std::stod( strs2[i] ));
            else
              speciesNames.emplace_back( strs2[i] );
          }
        }

        GEOSX_ERROR_IF( num1 != speciesNames.size() || num1 != stoichs.size() || speciesName != speciesNames[0], "Internal error when reading database" );

        bool notFound = 0;

        speciesIndices.resize( num1 );
        speciesIndices[0] = count;

        for( int i = 1; i < num1; i++ )
        {
          auto it = primarySpeciesMap.find( speciesNames[i] );
          if( it != primarySpeciesMap.end())
          {
            speciesIndices[i] = it->second;
          }
          else
          {
            notFound = 1;
            break;
          }
        }

        if( notFound )
        {
          speciesIndices.clear();
        }
        else
        {
          logKs.clear();
          while( is.getline( buf, buf_size ))
          {
            std::string str2( buf );
            auto found2 = str2.find( "*" );
            if( found2 != std::string::npos )
              break;

            string_array strs2 = Tokenize( str2, " " );

            for( int i = 0; i < strs2.size(); i++ )
              logKs.emplace_back( std::stod( strs2[i] ));

          }
        }
      }
    }
  }

  is.close();
}

// Not sure what is this for
REGISTER_CATALOG_ENTRY( ThermoDatabaseBase,
                        EQ36Database,
                        const Path &,
                        const string_array &,
                        const string_array & )

}

}
