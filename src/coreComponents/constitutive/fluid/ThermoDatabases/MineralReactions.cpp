
/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

/**
 * @file MineralReactions.cpp
 */

#include "MineralReactions.hpp"

using namespace std;

namespace geosx
{

using namespace stringutilities;

namespace constitutive
{

MineralReactions::MineralReactions( const Path & fileName,
                                    const string_array & basisSpeciesNames ):
  KineticReactionsBase( fileName )
{

  ReadMineralReactions( basisSpeciesNames );

}

void MineralReactions::ReadMineralReactions( const string_array & basisSpeciesNames )
{

  std::ifstream is( m_fileName );

  constexpr std::streamsize buf_size = 256;
  char buf[buf_size];

  unordered_map< string, int > basisSpeciesMap;

  bool isNotFirstLine = 0;

  localIndex numBasisSpecies = basisSpeciesNames.size();

  for( localIndex ic = 0; ic < numBasisSpecies; ++ic )
  {
    basisSpeciesMap[basisSpeciesNames[ic] ] = int(ic);
  }

  string mineralName;
  real64 MW{};
  real64 density{};
  real64 logK{};
  real64 E{};
  real64 rateConst{};
  array1d< real64 > stochs;
  array1d< localIndex > basisSpeciesIndices;
  array1d< string > speciesNames;

  while( is.getline( buf, buf_size ))
  {
    std::string const str( buf );

    {
      auto found = str.find( "+-------" );
      if( found!=std::string::npos && isNotFirstLine )
      {

        KineticReaction entry;
        entry.name = mineralName;
        entry.MW = MW;
        entry.density = density;
        entry.logK = logK;
        entry.E = E;
        entry.rateConst = rateConst;

        basisSpeciesIndices.resize( speciesNames.size()-1 );

        array1d< real64 > stochsNew;

        for( localIndex ic = 1; ic < speciesNames.size(); ++ic )
        {

          auto it = basisSpeciesMap.find( speciesNames[ic] );
          if( it != basisSpeciesMap.end())
          {
            basisSpeciesIndices[ic-1] = it->second;
            stochsNew.emplace_back( stochs[ic] );
          }
          else
          {
            GEOSX_ERROR( "Internal error when reading mineral reactions" );
          }

        }

        entry.basisSpeciesIndices = basisSpeciesIndices;
        entry.stochs = stochsNew;

        m_kineticReactions.emplace_back( entry );

        stochs.clear();
        speciesNames.clear();
        basisSpeciesIndices.clear();

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
      auto found = str.find( "density" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        density = std::stod( strs[3] );
      }
    }

    {

      auto found = str.find( "name" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        mineralName = strs[3];
      }

    }

    {
      auto found = str.find( "Ea" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        E = std::stod( strs[3] );
      }
    }


    {
      auto found = str.find( "logK" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        logK = std::stod( strs[3] );
      }
    }


    {
      auto found = str.find( "rateConst" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        rateConst = std::stod( strs[3] );
      }
    }

    {

      auto found = str.find( "dissociation reaction" );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        localIndex num1 = localIndex( std::stoi( strs[0] ));

        while( is.getline( buf, buf_size ))
        {

          std::string str2( buf );
          auto found2 = str2.find( "* End reaction" );
          if( found2 != std::string::npos )
            break;

          string_array strs2 = Tokenize( str2, " " );
          //localIndex num2 = strs2.size();

          for( localIndex i = 0; i < strs2.size(); ++i )
          {
            if( i % 2 == 0 )
              stochs.emplace_back( std::stod( strs2[i] ));
            else
              speciesNames.emplace_back( strs2[i] );
          }

        }

        GEOSX_ERROR_IF( num1 != speciesNames.size() || num1 != stochs.size(), "Internal error when reading mineral reactions" );

      }
    }

    isNotFirstLine = 1;

  }

  is.close();

}

REGISTER_CATALOG_ENTRY( KineticReactionsBase,
                        MineralReactions,
                        const Path &, const string_array & )

}

}
