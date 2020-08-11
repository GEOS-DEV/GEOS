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
 * @file GraphFromText.cpp
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "common/Path.hpp"
#include "GraphFromText.hpp"
#include "managers/Functions/FunctionManager.hpp"
#include "common/Logger.hpp"

namespace geosx
{

using namespace dataRepository;


GraphFromText::GraphFromText( const std::string & name,
                              Group * const parent ):
  GraphBase( name, parent ),m_file("")
  {
     setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

    // This enables logLevel filtering
    enableLogLevelInput();

    registerWrapper( viewKeyStruct::fileString, &m_file )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the graph file." );
    
  }


GraphFromText::~GraphFromText()
{}

void GraphFromText::GenerateGraph()
{
  std::ifstream infile(m_file);
  std::string line;
  while (std::getline(infile, line))
  {
    std::istringstream iss(line);
    int a, b, c;
    std::cout<<"test";
    if (!(iss >> a >> b >> c)) { break; } // error
    std::cout<<a<<" "<<b<<" "<<c<<"\n" ;   
  }
}

REGISTER_CATALOG_ENTRY( GraphBase, GraphFromText, std::string const &, Group * const )
} /* namespace geosx */
