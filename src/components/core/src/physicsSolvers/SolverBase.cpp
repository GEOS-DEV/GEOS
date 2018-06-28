// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/*
 * SolverBase.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: rrsettgast
 */

#include "SolverBase.hpp"


namespace geosx
{

using namespace dataRepository;

SolverBase::SolverBase( std::string const & name,
                        ManagedGroup * const parent ):
  ManagedGroup( name, parent ),
  m_verboseLevel(0)
{
  this->RegisterViewWrapper( "verbosity", &m_verboseLevel, 0 );
}

SolverBase::~SolverBase()
{}

SolverBase::CatalogInterface::CatalogType& SolverBase::GetCatalog()
{
  static SolverBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void SolverBase::FillDocumentationNode()
{


  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName(this->CatalogName());    // If this method lived in Managed
                                            // groups, this could be done
                                            // automatically
  docNode->setSchemaType("Node");

  docNode->AllocateChildNode( keys::courant,
                              keys::courant,
                              -1,
                              "real64",
                              "real64",
                              "courant Number",
                              "courant Number",
                              "0.7",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::maxDt,
                              keys::maxDt,
                              -1,
                              "real64",
                              "real64",
                              "Maximum Stable Timestep",
                              "Maximum Stable Timestep",
                              "0.0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::verboseLevelString,
                              viewKeyStruct::verboseLevelString,
                              -1,
                              "integer",
                              "integer",
                              "verbosity level",
                              "verbosity level",
                              "0",
                              "",
                              0,
                              1,
                              0 );


}



void SolverBase::TimeStep( real64 const& time_n,
                           real64 const& dt,
                           const int cycleNumber,
                           ManagedGroup * domain )
{
}



void SolverBase::CreateChild( string const & childKey, string const & childName )
{
  if( CatalogInterface::hasKeyName(childKey) )
  {
    std::cout << "Adding Solver of type " << childKey << ", named " << childName << std::endl;
    this->RegisterGroup( childName, CatalogInterface::Factory( childKey, childName, this ) );
  }
}


//void SolverBase::Initialize( dataRepository::ManagedGroup& /*domain*/ )
//{
//  *(this->getData<real64>(keys::courant)) =
// std::numeric_limits<real64>::max();
//}

} /* namespace ANST */
