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

#ifndef GEOSX_TESTSINGLEFLOWUTILS_HPP
#define GEOSX_TESTSINGLEFLOWUTILS_HPP

#include "codingUtilities/UnitTestUtilities.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "mesh/MeshManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVM.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"

namespace geosx
{

namespace testing
{

using namespace geosx::constitutive;

void checkDerivative( real64 const valueEps,
                      real64 const value,
                      real64 const deriv,
                      real64 const eps,
                      real64 const relTol,
                      real64 const absTol,
                      string const & name,
                      string const & var )
{
  real64 const numDeriv = (valueEps - value) / eps;
  checkRelativeError( deriv, numDeriv, relTol, absTol, "d(" + name + ")/d(" + var + ")" );
}

void checkDerivative( real64 const valueEps,
                      real64 const value,
                      real64 const deriv,
                      real64 const eps,
                      real64 const relTol,
                      string const & name,
                      string const & var )
{ return checkDerivative( valueEps, value, deriv, eps, relTol, DEFAULT_ABS_TOL, name, var ); }


void fillNumericalJacobian( arrayView1d< real64 const > const & residual,
                            arrayView1d< real64 const > const & residualOrig,
                            globalIndex const dofIndex,
                            real64 const eps,
                            CRSMatrixView< real64, globalIndex const > const & jacobian )
{
  forAll< parallelDevicePolicy<> >( residual.size(), [=] GEOSX_HOST_DEVICE ( localIndex const row )
  {
    real64 const dRdX = ( residual[row] - residualOrig[row] ) / eps;
    if( fabs( dRdX ) > 0.0 )
    {
      jacobian.addToRow< parallelDeviceAtomic >( row, &dofIndex, &dRdX, 1 );
    }
  } );
}

void setupProblemFromXML( ProblemManager & problemManager, char const * const xmlInput )
{
  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult xmlResult = xmlDocument.load_buffer( xmlInput, strlen( xmlInput ) );
  if( !xmlResult )
  {
    GEOSX_LOG_RANK_0( "XML parsed with errors!" );
    GEOSX_LOG_RANK_0( "Error description: " << xmlResult.description());
    GEOSX_LOG_RANK_0( "Error offset: " << xmlResult.offset );
  }

  int mpiSize = MpiWrapper::commSize( MPI_COMM_GEOSX );

  dataRepository::Group & commandLine =
    problemManager.getGroup< dataRepository::Group >( problemManager.groupKeys.commandLine );

  commandLine.registerWrapper< integer >( problemManager.viewKeys.xPartitionsOverride.key() ).
    setApplyDefaultValue( mpiSize );

  xmlWrapper::xmlNode xmlProblemNode = xmlDocument.child( dataRepository::keys::ProblemManager );
  problemManager.processInputFileRecursive( xmlProblemNode );

  DomainPartition & domain = problemManager.getDomainPartition();

  constitutive::ConstitutiveManager & constitutiveManager = domain.getConstitutiveManager();
  xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child( constitutiveManager.getName().c_str());
  constitutiveManager.processInputFileRecursive( topLevelNode );

  MeshManager & meshManager = problemManager.getGroup< MeshManager >( problemManager.groupKeys.meshManager );
  meshManager.generateMeshLevels( domain );

  ElementRegionManager & elementManager = domain.getMeshBody( 0 ).getBaseDiscretization().getElemManager();
  topLevelNode = xmlProblemNode.child( elementManager.getName().c_str());
  elementManager.processInputFileRecursive( topLevelNode );

  problemManager.problemSetup();
  problemManager.applyInitialConditions();
}

void testMobilityNumericalDerivatives( SinglePhaseFVM< SinglePhaseBase > & solver,
                                       DomainPartition & domain,
                                       bool const isThermal,
                                       real64 const perturbParameter,
                                       real64 const relTol )
{
  solver.forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                         [&]( string const,
                                              MeshLevel & mesh,
                                              arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elementRegionManager = mesh.getElemManager();
    elementRegionManager.forElementSubRegions( regionNames,
                                               [&]( localIndex const,
                                                    ElementSubRegionBase & subRegion )
    {
      SCOPED_TRACE( subRegion.getParent().getParent().getName() + "/" + subRegion.getName() );

      arrayView1d< real64 > const pres =
        subRegion.getField< fields::flow::pressure >();
      arrayView1d< real64 const > const pres_n =
        subRegion.getField< fields::flow::pressure_n >();

      arrayView1d< real64 > const mob =
        subRegion.getField< fields::flow::mobility >();

      arrayView1d< real64 > const dMob_dPres =
        subRegion.getField< fields::flow::dMobility_dPressure >();

      // reset the solver state to zero out variable updates
      solver.resetStateToBeginningOfStep( domain );

      // make a copy of unperturbed values of component fractions
      array1d< real64 > mobOrig( subRegion.size() );
      mobOrig.setValues< serialPolicy >( mob );

      // Step 1: update pressure and check derivatives

      {
        // perturb pressure in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          real64 const dP = perturbParameter * ( pres[ei] + perturbParameter );
          pres[ei] += dP;
        } );

        solver.updateFluidState( subRegion );

        // check values in each cell
        forAll< serialPolicy >( subRegion.size(), [=, &mobOrig] ( localIndex const ei )
        {
          SCOPED_TRACE( "Element " + std::to_string( ei ) );

          real64 const delta = pres[ei] - pres_n[ei];
          checkDerivative( mob[ei],
                           mobOrig[ei],
                           dMob_dPres[ei],
                           delta,
                           relTol,
                           "mob",
                           "Pres" );
        } );
      }

      // Step 2: update temperature and check derivatives

      if( isThermal )
      {
        arrayView1d< real64 > const temp =
          subRegion.getField< fields::flow::temperature >();
        arrayView1d< real64 const > const temp_n =
          subRegion.getField< fields::flow::temperature_n >();

        arrayView1d< real64 > const dMob_dTemp =
          subRegion.getField< fields::flow::dMobility_dTemperature >();

        // reset the solver state to zero out variable updates (resetting the whole domain is overkill...)
        solver.resetStateToBeginningOfStep( domain );

        // perturb temperature in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          real64 const dT = perturbParameter * ( temp[ei] + perturbParameter );
          temp[ei] += dT;
        } );

        solver.updateFluidState( subRegion );

        // check values in each cell
        forAll< serialPolicy >( subRegion.size(), [=, &mobOrig] ( localIndex const ei )
        {
          SCOPED_TRACE( "Element " + std::to_string( ei ) );

          real64 const delta = temp[ei] - temp_n[ei];
          checkDerivative( mob[ei],
                           mobOrig[ei],
                           dMob_dTemp[ei],
                           delta,
                           relTol,
                           "mob",
                           "Temp" );
        } );
      }

    } );
  } );
}

template< typename SINGLE_PHASE_SOLVER, typename LAMBDA >
void fillCellCenteredNumericalJacobian( SINGLE_PHASE_SOLVER & solver,
                                        DomainPartition & domain,
                                        bool const isThermal,
                                        real64 const perturbParameter,
                                        arrayView1d< real64 > residual,
                                        arrayView1d< real64 > residualOrig,
                                        CRSMatrixView< real64, globalIndex > jacobian,
                                        CRSMatrixView< real64, globalIndex > jacobianFD,
                                        LAMBDA assembleFunction )
{
  DofManager const & dofManager = solver.getDofManager();
  string const elemDofKey = dofManager.getKey( SINGLE_PHASE_SOLVER::viewKeyStruct::elemDofFieldString() );

  solver.forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                       MeshLevel & mesh,
                                                                       arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
      arrayView1d< globalIndex const > const & elemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( elemDofKey );

      arrayView1d< real64 > const pres =
        subRegion.getField< fields::flow::pressure >();

      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        if( elemGhostRank[ei] >= 0 )
        {
          continue;
        }

        solver.resetStateToBeginningOfStep( domain );

        // Step 1: compute numerical derivatives wrt pressure

        {
          pres.move( hostMemorySpace, true ); // to get the correct pres after reset
          real64 const dP = perturbParameter * ( pres[ei] + perturbParameter );
          pres[ei] += dP;
#if defined(GEOSX_USE_CUDA)
          pres.move( parallelDeviceMemorySpace, false );
#endif

          solver.updateState( domain );

          residual.zero();
          jacobian.zero();
          assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

          fillNumericalJacobian( residual.toViewConst(),
                                 residualOrig.toViewConst(),
                                 elemDofNumber[ei],
                                 dP,
                                 jacobianFD.toViewConstSizes() );
        }

        // Step 2: compute numerical derivatives wrt temperature

        if( isThermal )
        {
          arrayView1d< real64 > const temp =
            subRegion.getField< fields::flow::temperature >();
          temp.move( hostMemorySpace, false );

          solver.resetStateToBeginningOfStep( domain );

          temp.move( hostMemorySpace, true );
          real64 const dT = perturbParameter * ( temp[ei] + perturbParameter );
          temp[ei] += dT;
#if defined(GEOSX_USE_CUDA)
          temp.move( parallelDeviceMemorySpace, false );
#endif

          // here, we make sure that rock internal energy is updated
          // in other words, a call to updateFluidState would not work
          solver.updateState( domain );

          residual.zero();
          jacobian.zero();
          assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

          fillNumericalJacobian( residual.toViewConst(),
                                 residualOrig.toViewConst(),
                                 elemDofNumber[ei] + 1,
                                 dT,
                                 jacobianFD.toViewConstSizes() );
        }

      }
    } );
  } );
}


} // namespace testing

} // namespace geosx

#endif //GEOSX_TESTSINGLEFLOWUTILS_HPP
