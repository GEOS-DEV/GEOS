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
 * @file SolidMechanicsPostEquilibrationStep.cpp
 */

#include "SolidMechanicsPostEquilibrationStep.hpp"

#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "mainInterface/ProblemManager.hpp"

namespace geosx
{

using namespace constitutive;
using namespace dataRepository;

SolidMechanicsPostEquilibrationStep::SolidMechanicsPostEquilibrationStep( const string & name,
                                                                          Group * const parent ):
  TaskBase( name, parent ),
  m_solidSolverName()
{
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::solidSolverNameString(), &m_solidSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the solid mechanics solver" );
}

SolidMechanicsPostEquilibrationStep::~SolidMechanicsPostEquilibrationStep()
{}

void SolidMechanicsPostEquilibrationStep::postProcessInput()
{
  ProblemManager & problemManager = this->getGroupByPath< ProblemManager >( "/Problem" );
  PhysicsSolverManager & physicsSolverManager = problemManager.getPhysicsSolverManager();

  GEOSX_THROW_IF( !physicsSolverManager.hasGroup( m_solidSolverName ),
                  GEOSX_FMT( "Task {}: physics solver named {} not found",
                             getName(), m_solidSolverName ),
                  InputError );

  m_solidSolver = &physicsSolverManager.getGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );
}

bool SolidMechanicsPostEquilibrationStep::execute( real64 const time_n,
                                                   real64 const GEOSX_UNUSED_PARAM( dt ),
                                                   integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                                   integer const GEOSX_UNUSED_PARAM( eventCounter ),
                                                   real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                                                   DomainPartition & domain )
{
  m_solidSolver->forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "Task `{}`: at time {}s, physics solver `{}` sets total displacement and velocity for all the nodes of the mesh",
                                          getName(), time_n, m_solidSolverName ) );

    NodeManager & nodeManager = mesh.getNodeManager();
    ElementRegionManager & elementRegionManager = mesh.getElemManager();

    // Step 1: zero out velocity, incremental displacement, and displacement after the equilibration step
    arrayView2d< real64, nodes::VELOCITY_USD > const velocity = nodeManager.velocity();
    arrayView2d< real64, nodes::INCR_DISPLACEMENT_USD > const incrementalDisp = nodeManager.incrementalDisplacement();
    arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const disp = nodeManager.totalDisplacement();
    velocity.zero();
    incrementalDisp.zero();
    disp.zero();

    // Step 2: for plasticity, re-activate the yield surface for a plastic calculation with some known consolidation condition
    elementRegionManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                                       [&]( localIndex const,
                                                                            CellElementSubRegion & subRegion )
    {
      string const & solidMaterialName = subRegion.getReference< string >( SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString() );
      Group & constitutiveModels = subRegion.getGroup( ConstitutiveManager::groupKeyStruct::constitutiveModelsString() );

      GEOSX_LOG_LEVEL_RANK_0( 2, GEOSX_FMT( "Task `{}`: at time {}s, solid model `{}` performs a post-equilibration step on subRegion `{}`. "
                                            "This step is only performed by plasticity solid models.",
                                            getName(), time_n, solidMaterialName, subRegion.getName() ) );

      SolidBase & constitutiveRelation = constitutiveModels.getGroup< SolidBase >( solidMaterialName );
      constitutiveRelation.applyPostEquilibrationStep();
    } );
  } );

  return false;
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        SolidMechanicsPostEquilibrationStep,
                        string const &, dataRepository::Group * const )

} /* namespace geosx */
