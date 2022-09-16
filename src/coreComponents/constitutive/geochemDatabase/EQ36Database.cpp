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
  std::ifstream EQ36File ( m_databaseFileName );
  constexpr std::streamsize buffer_size = 256;
  char fileLine[buffer_size];

  while( EQ36File.getline( fileLine, buffer_size ))
  {
    std::string fileLineString ( fileLine );
    auto found = fileLineString.find( "Temperature grid" );
    if( found!=std::string::npos )
      break;
  }

  // read and store temperature information
  readActivityCoefficientBlock( EQ36File, buffer_size, "Pressure grid", m_actCoeffParameters.temperatures ) 
  // read and ignore pressure and pressure envelope information
  // As of now, we don't make use of pressure and hence I am not storing it
  readAndIgnoreActivityCoefficientBlock( EQ36File, buffer_size, "Pressure envelope") 
  readAndIgnoreActivityCoefficientBlock( EQ36File, buffer_size, "Debye-Huckel A_gamma") 
  // read and store Debye-Huckel A information
  readActivityCoefficientBlock( EQ36File, buffer_size, "Debye-Huckel A_H", m_actCoeffParameters.DebyeHuckelAs ) 
  // read and ignore Debye-Huckel A_H information
  readAndIgnoreActivityCoefficientBlock( EQ36File, buffer_size, "Debye-Huckel B_gamma") 
  // read and store Debye-Huckel B information and convert it in the correct units
  readActivityCoefficientBlock( EQ36File, buffer_size, "Debye-Huckel B_H", m_actCoeffParameters.DebyeHuckelBs ) 
  m_actCoeffParameters.DebyeHuckelBs = m_actCoeffParameters.DebyeHuckelBs*1e8
  // read and ignore Debye-Huckel B_H information
  readAndIgnoreActivityCoefficientBlock( EQ36File, buffer_size, "B-dot") 
  // read and store WATEQB-Dot information
  readActivityCoefficientBlock( EQ36File, buffer_size, "B-dot_H", m_actCoeffParameters.WATEQBDots ) 

  // if the size of the variables to be interpolated are not the same, throw an error
  GEOSX_ERROR_IF(
    m_actCoeffParameters.temperatures.size() != m_actCoeffParameters.DebyeHuckelAs.size() || m_actCoeffParameters.temperatures.size() != m_actCoeffParameters.DebyeHuckelBs.size() || m_actCoeffParameters.temperatures.size() != m_actCoeffParameters.WATEQBDots.size(),
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


void EQ36Database::readActivityCoefficientBlock( std::ifstream const EQ36File, 
                                                 std::streamsize const bufferSize,
                                                 string const nextBlockString,
                                                 arraySlice1d< real64 > const & readVariable ) const

{

  // This function reads a block from the EQ3/6 Database file and returns the read enteries in a vector.
  char fileLine[bufferSize];
  // I am not sure if information on which line was read is preserved when passing information through a function. 
  // That is critical because if the file is going to be read from the first line every single time then this won't work
  while( EQ36File.getline( fileLine, bufferSize ))         //read the line
  {
    std::string fileLineString ( fileLine );
    auto found = fileLineString.find( nextBlockString );   // can you find the required string in the line
    if( found!=std::string::npos )                         // if yes, you have reached the next block and should exit
      break;

    string_array lineEntries = Tokenize( str, " " );      // separate the different numeric enteries and store them
    for( int i = 0; i < lineEntries.size(); i++ )
    {
      readVariable.emplace_back( std::stod( lineEntries[i] ));
    }
  }

}

void EQ36Database::readAndIgnoreActivityCoefficientBlock( std::ifstream const EQ36File, 
                                                          std::streamsize const bufferSize,
                                                          string const nextBlockString ) const

{

  // This function reads a block from the EQ3/6 Database file and ignores it as it isn't useful for us
  char fileLine[bufferSize];
  // I am not sure if information on which line was read is preserved when passing information through a function. 
  // That is critical becuase if the file is going to be read from the first line every single time then this won't work
  while( EQ36File.getline( fileLine, bufferSize ))
  {
    std::string fileLineString ( fileLine );
    auto found = fileLineString.find( nextBlockString );
    if( found!=std::string::npos )
      break;

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
