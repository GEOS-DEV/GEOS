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

/*
 * ContactSolverBase.cpp
 */

#include "ContactSolverBase.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactSelector.hpp"
#include "constitutive/solid/ElasticIsotropic.hpp"
#include "finiteElement/elementFormulations/FiniteElementBase.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
#include "mesh/DomainPartition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/contact/ContactExtrinsicData.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

ContactSolverBase::ContactSolverBase( const string & name,
                                      Group * const parent ):
  SolverBase( name, parent ),
  m_solidSolver( nullptr )
{
  registerWrapper( viewKeyStruct::solidSolverNameString(), &m_solidSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the solid mechanics solver in the rock matrix" );

  registerWrapper( viewKeyStruct::contactRelationNameString(), &m_contactRelationName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of contact relation to enforce constraints on fracture boundary." );

  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );
}


void ContactSolverBase::postProcessInput()
{
  m_solidSolver = &this->getParent().getGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );
  SolverBase::postProcessInput();
}

void ContactSolverBase::registerDataOnMesh( dataRepository::Group & meshBodies )
{
  using namespace extrinsicMeshData::contact;

  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getMeshLevel( 0 );

    ElementRegionManager & elemManager = meshLevel.getElemManager();
    {
      elemManager.forElementRegions< SurfaceElementRegion >( [&] ( SurfaceElementRegion & region )
      {
        region.forElementSubRegions< SurfaceElementSubRegion >( [&]( SurfaceElementSubRegion & subRegion )
        {
          subRegion.registerExtrinsicData< dispJump >( getName() ).
            reference().resizeDimension< 1 >( 3 );

          subRegion.registerExtrinsicData< deltaDispJump >( getName() ).
            reference().resizeDimension< 1 >( 3 );

          subRegion.registerExtrinsicData< oldDispJump >( getName() ).
            reference().resizeDimension< 1 >( 3 );

          subRegion.registerExtrinsicData< traction >( getName() ).
            reference().resizeDimension< 1 >( 3 );

          subRegion.registerWrapper< array1d< integer > >( viewKeyStruct::fractureStateString() ).
            setPlotLevel( PlotLevel::LEVEL_0 ).
            setRegisteringObjects( this->getName()).
            setDescription( "An array that holds the fracture state." );
          initializeFractureState( subRegion, viewKeyStruct::fractureStateString() );

          subRegion.registerWrapper< array1d< integer > >( viewKeyStruct::oldFractureStateString() ).
            setPlotLevel( PlotLevel::NOPLOT ).
            setRegisteringObjects( this->getName()).
            setDescription( "An array that holds the fracture state." );
          initializeFractureState( subRegion, viewKeyStruct::oldFractureStateString() );

        } );
      } );
    }
  } );
}


void ContactSolverBase::initializeFractureState( SurfaceElementSubRegion & subRegion,
                                                 string const & fieldName ) const
{
  GEOSX_MARK_FUNCTION;
  arrayView1d< integer > const & fractureState = subRegion.getReference< array1d< integer > >( fieldName );
  fractureState.setValues< parallelHostPolicy >( FractureState::Stick );
}

real64 ContactSolverBase::solverStep( real64 const & time_n,
                                      real64 const & dt,
                                      int const cycleNumber,
                                      DomainPartition & domain )
{
  real64 dtReturn = dt;

  implicitStepSetup( time_n,
                     dt,
                     domain );

  setupSystem( domain,
               m_dofManager,
               m_localMatrix,
               m_rhs,
               m_solution );

  // currently the only method is implicit time integration
  dtReturn = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  // final step for completion of timestep. Typically secondary variable updates and cleanup.
  implicitStepComplete( time_n, dtReturn, domain );

  return dtReturn;
}

void ContactSolverBase::computeFractureStateStatistics( MeshLevel const & mesh,
                                                        globalIndex & numStick,
                                                        globalIndex & numSlip,
                                                        globalIndex & numOpen ) const
{
  ElementRegionManager const & elemManager = mesh.getElemManager();

  array1d< globalIndex > localCounter( 3 );

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & subRegion )
  {
    arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
    arrayView1d< integer const > const & fractureState = subRegion.getReference< array1d< integer > >( viewKeyStruct::fractureStateString() );

    RAJA::ReduceSum< parallelHostReduce, localIndex > stickCount( 0 ), slipCount( 0 ), openCount( 0 );
    forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
    {
      if( ghostRank[kfe] < 0 )
      {
        switch( fractureState[kfe] )
        {
          case FractureState::Stick:
            {
              stickCount += 1;
              break;
            }
          case FractureState::NewSlip:
          case FractureState::Slip:
            {
              slipCount += 1;
              break;
            }
          case FractureState::Open:
            {
              openCount += 1;
              break;
            }
        }
      }
    } );

    localCounter[0] += stickCount.get();
    localCounter[1] += slipCount.get();
    localCounter[2] += openCount.get();
  } );

  array1d< globalIndex > totalCounter( 3 );

  MpiWrapper::allReduce( localCounter.data(),
                         totalCounter.data(),
                         3,
                         MPI_SUM,
                         MPI_COMM_GEOSX );

  numStick = totalCounter[0];
  numSlip  = totalCounter[1];
  numOpen  = totalCounter[2];
}

void ContactSolverBase::outputConfigurationStatistics( DomainPartition const & domain ) const
{
  if( getLogLevel() >=1 )
  {
    globalIndex numStick = 0;
    globalIndex numSlip  = 0;
    globalIndex numOpen  = 0;

    forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                 MeshLevel const & mesh,
                                                 arrayView1d< string const > const & )
    {
      computeFractureStateStatistics( mesh, numStick, numSlip, numOpen );

      GEOSX_LOG_RANK_0( GEOSX_FMT( "  Number of element for each fracture state:"
                                   " stick: {:12} | slip:  {:12} | open:  {:12}",
                                   numStick, numSlip, numOpen ) );
    } );
  }
}

void ContactSolverBase::applyBoundaryConditions( real64 const time,
                                                 real64 const dt,
                                                 DomainPartition & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->applyBoundaryConditions( time,
                                          dt,
                                          domain,
                                          dofManager,
                                          localMatrix,
                                          localRhs );
}

real64 ContactSolverBase::explicitStep( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                        real64 const & dt,
                                        const int GEOSX_UNUSED_PARAM( cycleNumber ),
                                        DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  GEOSX_MARK_FUNCTION;
  GEOSX_ERROR( "ExplicitStep non available for contact solvers." );
  return dt;
}



void ContactSolverBase::synchronizeFractureState( DomainPartition & domain ) const
{
  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::fractureStateString() ) );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & )
  {
    CommunicationTools::getInstance().synchronizeFields( fieldNames,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );
}

} /* namespace geosx */
