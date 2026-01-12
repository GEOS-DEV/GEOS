/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file WellManager.cpp
 */

#include "WellManager.hpp"

#include "mesh/DomainPartition.hpp"
#include "mesh/PerforationFields.hpp"
#include "mesh/WellElementRegion.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "physicsSolvers/fluidFlow/wells/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBase.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"
#include "fileIO/Outputs/OutputBase.hpp"
#include "functions/FunctionManager.hpp"
namespace geos
{

using namespace dataRepository;
using namespace fields;

WellManager::WellManager( string const & name,
                                Group * const parent )
  : PhysicsSolverBase( name, parent )  

{
  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE ); 
}
Group * WellManager::createChild( string const & childKey, string const & childName )
{
  static std::set< string > const childTypes = {
    keys::compositionalMultiphaseWell,
    keys::singlePhaseWell,
    PhysicsSolverBase::groupKeyStruct::linearSolverParametersString(),
    PhysicsSolverBase::groupKeyStruct::nonlinearSolverParametersString(),
  };
  GEOS_ERROR_IF( childTypes.count( childKey ) == 0,
                 CatalogInterface::unknownTypeError( childKey, getDataContext(), childTypes ),
                 getDataContext() );
  if( childKey == keys::compositionalMultiphaseWell )
  {
    return &registerGroup< CompositionalMultiphaseWell >( childName );
  }
  else if( childKey == keys::singlePhaseWell )
  {
    return &registerGroup< SinglePhaseWell >( childName );
  }
  else
  {
    PhysicsSolverBase::createChild( childKey, childName );
    return nullptr;
  }
}

void WellManager::expandObjectCatalogs()
{
  createChild( keys::compositionalMultiphaseWell, keys::compositionalMultiphaseWell );
  createChild( keys::singlePhaseWell, keys::singlePhaseWell );
}
WellSolverBase & WellManager::getWell( WellElementSubRegion const & subRegion )
{
  return this->getGroup< WellSolverBase >( subRegion.getWellControlsName());
}

void WellManager::implicitStepSetup( real64 const & time_n,
                                                     real64 const & dt,
                                                     DomainPartition & domain )
{

  forDiscretizationOnMeshTargets ( domain.getMeshBodies(), [&] ( string const &,
                                                                 MeshLevel & mesh,
                                                                 string_array const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {

      WellSolverBase & well = getWell( subRegion );
      well.implicitStepSetup( time_n, dt, domain );
    } )
    ;
  } );
}
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, WellManager, string const &, Group * const )
} // namespace geos
