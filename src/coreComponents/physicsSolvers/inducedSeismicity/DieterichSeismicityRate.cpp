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
 * @file DieterichSeismicityRate.cpp
 */

// Source includes
#include "DieterichSeismicityRate.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"

namespace geos
{

namespace dataRepository
{
namespace keys
{}
}

using namespace dataRepository;
using namespace fields;
using namespace constitutive;

DieterichSeismicityRate::DieterichSeismicityRate( const string & name,
                                                  Group * const parent ):
  SeismicityRateBase( name, parent )
{
  this->registerWrapper( viewKeyStruct::directEffectString(), &m_directEffect ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Rate-and-state friction direct effect parameter" );
  this->registerWrapper( viewKeyStruct::backgroundStressingRateString(), &m_backgroundStressingRate ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Background stressing rate" );
}

DieterichSeismicityRate::~DieterichSeismicityRate()
{
  // TODO Auto-generated destructor stub
}

void DieterichSeismicityRate::registerDataOnMesh( Group & meshBodies )
{
  SeismicityRateBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      subRegion.registerField< inducedSeismicity::logDenom >( getName() );
    } );
  } );
}

real64 DieterichSeismicityRate::solverStep( real64 const & time_n,
                                            real64 const & dt,
                                            const int cycleNumber,
                                            DomainPartition & domain )
{
  // Save initial stress state on pre-defined fault orientations to field variables
  initializeFaultTraction( time_n, cycleNumber, domain );

  // Call member variable stress solver to update the stress state
  real64 dtStress = m_stressSolver->solverStep( time_n, dt, cycleNumber, domain );

  // Loop over subRegions to update stress on faults and solver for seismicity rate
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      // project new stress state to update stress on fault
      if( subRegion.hasWrapper( SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString() ) )
      {
        updateFaultTraction( subRegion );
      }

      // solve for the seismicity rate given new stresses on faults
      integralSolverStep( time_n, dtStress, subRegion );
    } );
  } );

  // return time step size achieved by stress solver
  return dtStress;
}

// Solve integral solution to ODE
void DieterichSeismicityRate::integralSolverStep( real64 const & time_n,
                                                  real64 const & dt,
                                                  ElementSubRegionBase & subRegion )
{
  solverHelper solverHelperStruct( subRegion );

  solverHelperStruct.computeSeismicityRate( subRegion, time_n, dt,
                                            m_directEffect, m_backgroundStressingRate );
}

void DieterichSeismicityRate::initializePreSubGroups()
{
  SeismicityRateBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 > const tempLogDenom = subRegion.getField< inducedSeismicity::logDenom >();
      tempLogDenom.setValues< parallelHostPolicy >( 0.0 );
    } );
  } );
}

void checkExpArgument( real64 arg )
{
  GEOS_UNUSED_VAR( arg );
  // TODO:
  // 1. CHECK IF CLOSE TO LITHOSTATIC PRESSURE
  // 2. CHECK IF STRESSING RATE IS TOO LARGE
  // 3. CHECK IF a IS TOO SMALL
}

REGISTER_CATALOG_ENTRY( SolverBase, DieterichSeismicityRate, string const &, Group * const )
} /* namespace geos */
