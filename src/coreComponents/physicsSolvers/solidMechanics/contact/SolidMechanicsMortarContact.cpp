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

/*
 * SolidMechanicsMortarContact.cpp
 */

#include "mesh/DomainPartition.hpp"
#include "SolidMechanicsMortarContact.hpp"

#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsConformingContactKernelsBase.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsALMKernels.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsALMKernelsBase.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsALMSimultaneousKernels.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsDisplacementJumpUpdateKernels.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsContactFaceBubbleKernels.hpp"
#include "physicsSolvers/solidMechanics/contact/LogLevelsInfo.hpp"
#include "physicsSolvers/solidMechanics/contact/ContactFields.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"

#include "constitutive/contact/FrictionSelector.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

SolidMechanicsMortarContact::SolidMechanicsMortarContact( const string & name,
                                                                                    Group * const parent ):
  ContactSolverBase( name, parent )
{}

SolidMechanicsMortarContact::~SolidMechanicsMortarContact()
{}

void SolidMechanicsMortarContact::registerDataOnMesh( dataRepository::Group & meshBodies )
{
  GEOS_MARK_FUNCTION;
  ContactSolverBase::registerDataOnMesh( meshBodies );
}

void SolidMechanicsMortarContact::setupDofs( DomainPartition const & domain,
                                             DofManager & dofManager ) const
{
  GEOS_MARK_FUNCTION;
  SolidMechanicsLagrangianFEM::setupDofs( domain, dofManager );
}

void SolidMechanicsMortarContact::setupSystem( DomainPartition & domain,
                                                            DofManager & dofManager,
                                                            CRSMatrix< real64, globalIndex > & localMatrix,
                                                            ParallelVector & rhs,
                                                            ParallelVector & solution,
                                                            bool const setSparsity )
{

  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( domain, dofManager, localMatrix, rhs, solution, setSparsity );

}

void SolidMechanicsMortarContact::implicitStepSetup( real64 const & time_n,
                                                     real64 const & dt,
                                                     DomainPartition & domain )
{

  GEOS_MARK_FUNCTION;
  SolidMechanicsLagrangianFEM::implicitStepSetup( time_n, dt, domain );

  //////////////////////////////////////////////////////////////////////////////////
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    std::cout << meshName << std::endl;
    elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion1 )
    {
      forAll< parallelHostPolicy >( subRegion1.size(), [=] ( localIndex const kfe )
      {
        std::cout << kfe << std::endl;
      });
    });
  } );
  //////////////////////////////////////////////////////////////////////////////////

  GEOS_ERROR("Mortar solver is not implemented yet. ");

}

void SolidMechanicsMortarContact::assembleSystem( real64 const time,
                                                               real64 const dt,
                                                               DomainPartition & domain,
                                                               DofManager const & dofManager,
                                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                               arrayView1d< real64 > const & localRhs )
{

  GEOS_MARK_FUNCTION;
  synchronizeFractureState( domain );

  SolidMechanicsLagrangianFEM::assembleSystem( time,
                                               dt,
                                               domain,
                                               dofManager,
                                               localMatrix,
                                               localRhs );

  //ParallelMatrix parallel_matrix;
  //parallel_matrix.create( localMatrix.toViewConst(), dofManager.numLocalDofs(), MPI_COMM_GEOS );
  //parallel_matrix.write("mech.mtx");

}

void SolidMechanicsMortarContact::implicitStepComplete( real64 const & time_n,
                                                                     real64 const & dt,
                                                                     DomainPartition & domain )
{

  SolidMechanicsLagrangianFEM::implicitStepComplete( time_n, dt, domain );

}

real64 SolidMechanicsMortarContact::calculateResidualNorm( real64 const & time,
                                                                        real64 const & dt,
                                                                        DomainPartition const & domain,
                                                                        DofManager const & dofManager,
                                                                        arrayView1d< real64 const > const & localRhs )
{

  GEOS_MARK_FUNCTION;
  real64 const solidResidualNorm = SolidMechanicsLagrangianFEM::calculateResidualNorm( time, dt, domain, dofManager, localRhs );

  return solidResidualNorm;

}

void SolidMechanicsMortarContact::applySystemSolution( DofManager const & dofManager,
                                                                    arrayView1d< real64 const > const & localSolution,
                                                                    real64 const scalingFactor,
                                                                    real64 const dt,
                                                                    DomainPartition & domain )
{

  GEOS_MARK_FUNCTION;
  SolidMechanicsLagrangianFEM::applySystemSolution( dofManager,
                                                    localSolution,
                                                    scalingFactor,
                                                    dt,
                                                    domain );

}

void SolidMechanicsMortarContact::updateState( DomainPartition & domain )
{
  GEOS_UNUSED_VAR( domain );
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SolidMechanicsMortarContact, string const &, Group * const )
} /* namespace geos */
