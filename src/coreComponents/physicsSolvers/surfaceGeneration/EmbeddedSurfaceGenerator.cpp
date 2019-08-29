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
 * @file EmbeddedSurfaceGenerator.cpp
 */

#include "EmbeddedSurfaceGenerator.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "MPI_Communications/SpatialPartition.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"

#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/EmbeddedSurfaceRegion.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"

#include <set>

namespace geosx
{

EmbeddedSurfaceGenerator::EmbeddedSurfaceGenerator(const std::string& name,
                                                   ManagedGroup * const parent ):
  SolverBase( name, parent ),
  m_failCriterion( 1 ),
  m_separableFaceSet()
{
  this->RegisterViewWrapper( viewKeyStruct::failCriterionString,
  &this->m_failCriterion,
  0 );

  this->RegisterViewWrapper( viewKeyStruct::fractureRegionNameString, &m_fractureRegionName, 0 )->
  setInputFlag(dataRepository::InputFlags::OPTIONAL)->
  setApplyDefaultValue("FractureRegion");
}

EmbeddedSurfaceGenerator::~EmbeddedSurfaceGenerator()
{}

void EmbeddedSurfaceGenerator::RegisterDataOnMesh( ManagedGroup * const MeshBody )
{
  // Matteo: to be filled in
  std::cout << "1. Register data on mesh \n";


}


real64 EmbeddedSurfaceGenerator::SolverStep( real64 const& time_n,
                             real64 const& dt,
                             integer const cycleNumber,
                             DomainPartition * domain )
{
  std::cout << "SolverStep \n";
  real64 test = 0;
  return test;
}


void EmbeddedSurfaceGenerator::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager )
{

}


void EmbeddedSurfaceGenerator::postRestartInitialization( ManagedGroup * const domain0 )
{
}


REGISTER_CATALOG_ENTRY( SolverBase,
                        EmbeddedSurfaceGenerator,
                        std::string const &, dataRepository::ManagedGroup * const )

} /* namespace geosx */
