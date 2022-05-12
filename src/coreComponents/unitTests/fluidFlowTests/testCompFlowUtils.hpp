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

#ifndef GEOSX_TESTCOMPFLOWUTILS_HPP
#define GEOSX_TESTCOMPFLOWUTILS_HPP

#include "codingUtilities/UnitTestUtilities.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "mesh/MeshManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"

namespace geosx
{

namespace testing
{

using namespace geosx::constitutive;
using namespace geosx::constitutive::multifluid;

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

template< int USD1, int USD2, int USD3 >
void checkDerivative( arraySlice1d< real64 const, USD1 > const & valueEps,
                      arraySlice1d< real64 const, USD2 > const & value,
                      arraySlice1d< real64 const, USD3 > const & deriv,
                      real64 const eps,
                      real64 const relTol,
                      real64 const absTol,
                      string const & name,
                      string const & var,
                      arrayView1d< string const > const & labels )
{
  localIndex const size = labels.size( 0 );

  for( localIndex i = 0; i < size; ++i )
  {
    checkDerivative( valueEps[i], value[i], deriv[i], eps, relTol, absTol,
                     name + "[" + labels[i] + "]", var );
  }
}

template< int DIM, int USD1, int USD2, int USD3, typename ... Args >
void checkDerivative( ArraySlice< real64 const, DIM, USD1 > const & valueEps,
                      ArraySlice< real64 const, DIM, USD2 > const & value,
                      ArraySlice< real64 const, DIM, USD3 > const & deriv,
                      real64 const eps,
                      real64 const relTol,
                      real64 const absTol,
                      string const & name,
                      string const & var,
                      arrayView1d< string const > const & labels,
                      Args ... label_lists )
{
  localIndex const size = labels.size( 0 );

  for( localIndex i = 0; i < size; ++i )
  {
    checkDerivative( valueEps[i], value[i], deriv[i], eps, relTol, absTol,
                     name + "[" + labels[i] + "]", var, label_lists ... );
  }
}

template< int DIM, int USD1, int USD2, int USD3, typename ... Args >
void checkDerivative( ArraySlice< real64 const, DIM, USD1 > const & valueEps,
                      ArraySlice< real64 const, DIM, USD2 > const & value,
                      ArraySlice< real64 const, DIM, USD3 > const & deriv,
                      real64 const eps,
                      real64 const relTol,
                      string const & name,
                      string const & var,
                      arrayView1d< string const > const & labels,
                      Args ... label_lists )
{ return checkDerivative( valueEps, value, deriv, eps, relTol, DEFAULT_ABS_TOL, name, var, labels, label_lists ... ); }

// invert compositional derivative array layout to move innermost slice on the top
// (this is needed so we can use checkDerivative() to check derivative w.r.t. for each compositional var)
template< int USD >
array1d< real64 > invertLayout( arraySlice1d< real64 const, USD > const & input,
                                localIndex N )
{
  array1d< real64 > output( N );
  for( int i = 0; i < N; ++i )
  {
    output[i] = input[i];
  }

  return output;
}

template< int USD >
array2d< real64 > invertLayout( arraySlice2d< real64 const, USD > const & input,
                                localIndex N1,
                                localIndex N2 )
{
  array2d< real64 > output( N2, N1 );

  for( localIndex i = 0; i < N1; ++i )
  {
    for( localIndex j = 0; j < N2; ++j )
    {
      output( j, i ) = input( i, j );
    }
  }

  return output;
}

template< int USD >
array3d< real64 > invertLayout( arraySlice3d< real64 const, USD > const & input,
                                localIndex N1,
                                localIndex N2,
                                localIndex N3 )
{
  array3d< real64 > output( N3, N1, N2 );

  for( localIndex i = 0; i < N1; ++i )
  {
    for( localIndex j = 0; j < N2; ++j )
    {
      for( localIndex k = 0; k < N3; ++k )
      {
        output( k, i, j ) = input( i, j, k );
      }
    }
  }

  return output;
}

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

  ElementRegionManager & elementManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getElemManager();
  topLevelNode = xmlProblemNode.child( elementManager.getName().c_str());
  elementManager.processInputFileRecursive( topLevelNode );

  problemManager.problemSetup();
  problemManager.applyInitialConditions();
}

void testCompositionNumericalDerivatives( CompositionalMultiphaseFVM & solver,
                                          DomainPartition & domain,
                                          real64 const perturbParameter,
                                          real64 const relTol )
{
  integer const numComp = solver.numFluidComponents();

  solver.forMeshTargets( domain.getMeshBodies(),
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

      string const & fluidName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
      arrayView1d< string const > const & components = fluid.componentNames();

      arrayView2d< real64, compflow::USD_COMP > const compDens =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity >();
      arrayView2d< real64 const, compflow::USD_COMP > const compDens_n =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity_n >();

      arrayView2d< real64, compflow::USD_COMP > const compFrac =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompFraction >();

      arrayView3d< real64, compflow::USD_COMP_DC > const dCompFrac_dCompDens =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::dGlobalCompFraction_dGlobalCompDensity >();

      // reset the solver state to zero out variable updates
      solver.resetStateToBeginningOfStep( domain );

      // make a copy of unperturbed values of component fractions
      array2d< real64, compflow::LAYOUT_COMP > compFracOrig( subRegion.size(), numComp );
      compFracOrig.setValues< serialPolicy >( compFrac );

      // update component density and check derivatives
      for( integer jc = 0; jc < numComp; ++jc )
      {
        // reset the solver state to zero out variable updates (resetting the whole domain is overkill...)
        solver.resetStateToBeginningOfStep( domain );

        // perturb a single component density in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          real64 const dRho = perturbParameter * ( compDens[ei][jc] + perturbParameter );
          compDens[ei][jc] += dRho;
        } );

        // recompute component fractions
        solver.updateComponentFraction( subRegion );

        // check values in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          SCOPED_TRACE( "Element " + std::to_string( ei ) );

          auto dZ_dRho = invertLayout( dCompFrac_dCompDens[ei].toSliceConst(), numComp, numComp );
          string var = "compDens[" + components[jc] + "]";

          real64 const delta = compDens[ei][jc] - compDens_n[ei][jc];
          checkDerivative( compFrac[ei].toSliceConst(),
                           compFracOrig[ei].toSliceConst(),
                           dZ_dRho[jc].toSliceConst(),
                           delta,
                           relTol,
                           "compFrac",
                           var,
                           components );
        } );
      }
    } );
  } );
}

void testPhaseVolumeFractionNumericalDerivatives( CompositionalMultiphaseFVM & solver,
                                                  DomainPartition & domain,
                                                  bool const isThermal,
                                                  real64 const perturbParameter,
                                                  real64 const relTol )
{
  integer const numComp = solver.numFluidComponents();
  integer const numPhase = solver.numFluidPhases();

  solver.forMeshTargets( domain.getMeshBodies(),
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

      string const & fluidName = subRegion.getReference< string >( CompositionalMultiphaseFVM::viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
      arrayView1d< string const > const & components = fluid.componentNames();
      arrayView1d< string const > const & phases = fluid.phaseNames();

      arrayView1d< real64 > const pres =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::pressure >();
      arrayView1d< real64 const > const pres_n =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::pressure_n >();

      arrayView2d< real64, compflow::USD_COMP > const compDens =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity >();
      arrayView2d< real64, compflow::USD_COMP > const compDens_n =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity_n >();

      arrayView2d< real64, compflow::USD_PHASE > const phaseVolFrac =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseVolumeFraction >();

      arrayView2d< real64, compflow::USD_PHASE > const dPhaseVolFrac_dPres =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::dPhaseVolumeFraction_dPressure >();

      arrayView3d< real64, compflow::USD_PHASE_DC > const dPhaseVolFrac_dCompDens =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::dPhaseVolumeFraction_dGlobalCompDensity >();

      // reset the solver state to zero out variable updates
      solver.resetStateToBeginningOfStep( domain );

      // make a copy of unperturbed values of component fractions
      array2d< real64, compflow::LAYOUT_PHASE > phaseVolFracOrig( subRegion.size(), numPhase );
      phaseVolFracOrig.setValues< serialPolicy >( phaseVolFrac );

      // Step 1: update pressure and check derivatives

      {
        // perturb pressure in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          real64 const dP = perturbParameter * ( pres[ei] + perturbParameter );
          pres[ei] += dP;
        } );

        // recompute component fractions
        solver.updateFluidState( subRegion );

        // check values in each cell
        forAll< serialPolicy >( subRegion.size(), [=, &phaseVolFracOrig] ( localIndex const ei )
        {
          SCOPED_TRACE( "Element " + std::to_string( ei ) );

          real64 const delta = pres[ei] - pres_n[ei];
          checkDerivative( phaseVolFrac[ei].toSliceConst(),
                           phaseVolFracOrig[ei].toSliceConst(),
                           dPhaseVolFrac_dPres[ei].toSliceConst(),
                           delta,
                           relTol,
                           "phaseVolFrac",
                           "Pres",
                           phases );
        } );
      }

      // Step 2: update component density and check derivatives

      for( integer jc = 0; jc < numComp; ++jc )
      {
        // reset the solver state to zero out variable updates (resetting the whole domain is overkill...)
        solver.resetStateToBeginningOfStep( domain );

        // perturb a single component density in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          real64 const dRho = perturbParameter * ( compDens[ei][jc] + perturbParameter );
          compDens[ei][jc] += dRho;
        } );

        // recompute component fractions
        solver.updateFluidState( subRegion );

        // check values in each cell
        forAll< serialPolicy >( subRegion.size(), [=, &phaseVolFracOrig] ( localIndex const ei )
        {
          SCOPED_TRACE( "Element " + std::to_string( ei ) );

          auto dS_dRho = invertLayout( dPhaseVolFrac_dCompDens[ei].toSliceConst(), numPhase, numComp );
          string var = "compDens[" + components[jc] + "]";

          real64 const delta = compDens[ei][jc] - compDens_n[ei][jc];
          checkDerivative( phaseVolFrac[ei].toSliceConst(),
                           phaseVolFracOrig[ei].toSliceConst(),
                           dS_dRho[jc].toSliceConst(),
                           delta,
                           relTol,
                           "phaseVolFrac",
                           var,
                           phases );
        } );
      }

      // Step 3: update temperature and check derivatives

      if( isThermal )
      {
        arrayView1d< real64 > const temp =
          subRegion.getExtrinsicData< extrinsicMeshData::flow::temperature >();
        arrayView1d< real64 > const temp_n =
          subRegion.getExtrinsicData< extrinsicMeshData::flow::temperature_n >();

        arrayView2d< real64, compflow::USD_PHASE > const dPhaseVolFrac_dTemp =
          subRegion.getExtrinsicData< extrinsicMeshData::flow::dPhaseVolumeFraction_dTemperature >();

        // reset the solver state to zero out variable updates
        solver.resetStateToBeginningOfStep( domain );

        // perturb temperature in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          real64 const dT = perturbParameter * ( temp[ei] + perturbParameter );
          temp[ei] += dT;
        } );

        // recompute all fluid properties
        solver.updateFluidState( subRegion );

        // check values in each cell
        forAll< serialPolicy >( subRegion.size(), [=, &phaseVolFracOrig] ( localIndex const ei )
        {
          SCOPED_TRACE( "Element " + std::to_string( ei ) );

          real64 const delta = temp[ei] - temp_n[ei];
          checkDerivative( phaseVolFrac[ei].toSliceConst(),
                           phaseVolFracOrig[ei].toSliceConst(),
                           dPhaseVolFrac_dTemp[ei].toSliceConst(),
                           delta,
                           relTol,
                           "phaseVolFrac",
                           "Temp",
                           phases );
        } );
      }

    } );
  } );
}

void testPhaseMobilityNumericalDerivatives( CompositionalMultiphaseFVM & solver,
                                            DomainPartition & domain,
                                            bool const isThermal,
                                            real64 const perturbParameter,
                                            real64 const relTol )
{
  integer const numComp = solver.numFluidComponents();
  integer const numPhase = solver.numFluidPhases();

  solver.forMeshTargets( domain.getMeshBodies(),
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

      string const & fluidName = subRegion.getReference< string >( CompositionalMultiphaseFVM::viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
      arrayView1d< string const > const & components = fluid.componentNames();
      arrayView1d< string const > const & phases = fluid.phaseNames();

      arrayView1d< real64 > const pres =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::pressure >();
      arrayView1d< real64 const > const pres_n =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::pressure_n >();

      arrayView2d< real64, compflow::USD_COMP > const compDens =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity >();
      arrayView2d< real64 const, compflow::USD_COMP > const compDens_n =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity_n >();

      arrayView2d< real64, compflow::USD_PHASE > const phaseMob =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseMobility >();

      arrayView2d< real64, compflow::USD_PHASE > const dPhaseMob_dPres =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::dPhaseMobility_dPressure >();

      arrayView3d< real64, compflow::USD_PHASE_DC > const dPhaseMob_dCompDens =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::dPhaseMobility_dGlobalCompDensity >();

      // reset the solver state to zero out variable updates
      solver.resetStateToBeginningOfStep( domain );

      // make a copy of unperturbed values of component fractions
      array2d< real64, compflow::LAYOUT_PHASE > phaseMobOrig( subRegion.size(), numPhase );
      phaseMobOrig.setValues< serialPolicy >( phaseMob );

      // Step 1: update pressure and check derivatives

      {
        // perturb pressure in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          real64 const dP = perturbParameter * ( pres[ei] + perturbParameter );
          pres[ei] += dP;
        } );

        // recompute component fractions
        solver.updateFluidState( subRegion );

        // check values in each cell
        forAll< serialPolicy >( subRegion.size(), [=, &phaseMobOrig] ( localIndex const ei )
        {
          SCOPED_TRACE( "Element " + std::to_string( ei ) );

          real64 const delta = pres[ei] - pres_n[ei];
          checkDerivative( phaseMob[ei].toSliceConst(),
                           phaseMobOrig[ei].toSliceConst(),
                           dPhaseMob_dPres[ei].toSliceConst(),
                           delta,
                           relTol,
                           "phaseMob",
                           "Pres",
                           phases );
        } );
      }

      // Step 2: update component density and check derivatives

      for( integer jc = 0; jc < numComp; ++jc )
      {
        // reset the solver state to zero out variable updates (resetting the whole domain is overkill...)
        solver.resetStateToBeginningOfStep( domain );

        // perturb a single component density in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          real64 const dRho = perturbParameter * ( compDens[ei][jc] + perturbParameter );
          compDens[ei][jc] += dRho;
        } );

        // recompute component fractions
        solver.updateFluidState( subRegion );

        // check values in each cell
        forAll< serialPolicy >( subRegion.size(), [=, &phaseMobOrig] ( localIndex const ei )
        {
          SCOPED_TRACE( "Element " + std::to_string( ei ) );

          auto dS_dRho = invertLayout( dPhaseMob_dCompDens[ei].toSliceConst(), numPhase, numComp );
          string var = "compDens[" + components[jc] + "]";

          real64 const delta = compDens[ei][jc] - compDens_n[ei][jc];
          checkDerivative( phaseMob[ei].toSliceConst(),
                           phaseMobOrig[ei].toSliceConst(),
                           dS_dRho[jc].toSliceConst(),
                           delta,
                           relTol,
                           "phaseMob",
                           var,
                           phases );
        } );
      }

      // Step 3: update temperature and check derivatives

      if( isThermal )
      {
        arrayView1d< real64 > const temp =
          subRegion.getExtrinsicData< extrinsicMeshData::flow::temperature >();
        arrayView1d< real64 const > const temp_n =
          subRegion.getExtrinsicData< extrinsicMeshData::flow::temperature_n >();

        arrayView2d< real64, compflow::USD_PHASE > const dPhaseMob_dTemp =
          subRegion.getExtrinsicData< extrinsicMeshData::flow::dPhaseMobility_dTemperature >();

        // reset the solver state to zero out variable updates (resetting the whole domain is overkill...)
        solver.resetStateToBeginningOfStep( domain );

        // perturb temperature in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          real64 const dT = perturbParameter * ( temp[ei] + perturbParameter );
          temp[ei] += dT;
        } );

        // recompute component fractions
        solver.updateFluidState( subRegion );

        // check values in each cell
        forAll< serialPolicy >( subRegion.size(), [=, &phaseMobOrig] ( localIndex const ei )
        {
          SCOPED_TRACE( "Element " + std::to_string( ei ) );

          real64 const delta = temp[ei] - temp_n[ei];
          checkDerivative( phaseMob[ei].toSliceConst(),
                           phaseMobOrig[ei].toSliceConst(),
                           dPhaseMob_dTemp[ei].toSliceConst(),
                           delta,
                           relTol,
                           "phaseMob",
                           "Temp",
                           phases );
        } );
      }

    } );
  } );
}

template< typename COMPOSITIONAL_SOLVER, typename LAMBDA >
void fillCellCenteredNumericalJacobian( COMPOSITIONAL_SOLVER & solver,
                                        DomainPartition & domain,
                                        bool const isThermal,
                                        real64 const perturbParameter,
                                        arrayView1d< real64 > residual,
                                        arrayView1d< real64 > residualOrig,
                                        CRSMatrixView< real64, globalIndex > jacobian,
                                        CRSMatrixView< real64, globalIndex > jacobianFD,
                                        LAMBDA assembleFunction )
{
  integer const numComp = solver.numFluidComponents();

  DofManager const & dofManager = solver.getDofManager();
  string const elemDofKey = dofManager.getKey( COMPOSITIONAL_SOLVER::viewKeyStruct::elemDofFieldString() );

  solver.forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
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
        subRegion.getExtrinsicData< extrinsicMeshData::flow::pressure >();
      arrayView2d< real64, compflow::USD_COMP > const compDens =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity >();

      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        if( elemGhostRank[ei] >= 0 )
        {
          continue;
        }

        // Warning! To compute a correct totalDensity, we have to reset compDens
        solver.resetStateToBeginningOfStep( domain );

        real64 totalDensity = 0.0;
        compDens.move( LvArray::MemorySpace::host, true ); // to get the correct compDens after reset
        for( integer ic = 0; ic < numComp; ++ic )
        {
          totalDensity += compDens[ei][ic];
        }

        // Step 1: compute numerical derivatives wrt pressure

        {
          pres.move( LvArray::MemorySpace::host, true ); // to get the correct pres after reset
          real64 const dP = perturbParameter * ( pres[ei] + perturbParameter );
          pres[ei] += dP;
#if defined(GEOSX_USE_CUDA)
          pres.move( LvArray::MemorySpace::cuda, false );
#endif


          mesh.getElemManager().forElementSubRegions( regionNames,
                                                      [&]( localIndex const,
                                                           ElementSubRegionBase & subRegion2 )
          {
            solver.updateFluidState( subRegion2 );
          } );

          residual.zero();
          jacobian.zero();
          assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

          fillNumericalJacobian( residual.toViewConst(),
                                 residualOrig.toViewConst(),
                                 elemDofNumber[ei],
                                 dP,
                                 jacobianFD.toViewConstSizes() );
        }

        // Step 2: compute numerical derivatives wrt component densities

        for( integer jc = 0; jc < numComp; ++jc )
        {
          solver.resetStateToBeginningOfStep( domain );

          compDens.move( LvArray::MemorySpace::host, true ); // to get the correct compDens after reset
          real64 const dRho = perturbParameter * totalDensity;
          compDens[ei][jc] += dRho;
#if defined(GEOSX_USE_CUDA)
          compDens.move( LvArray::MemorySpace::cuda, false );
#endif


          mesh.getElemManager().forElementSubRegions( regionNames,
                                                      [&]( localIndex const,
                                                           ElementSubRegionBase & subRegion2 )
          {
            solver.updateFluidState( subRegion2 );
          } );

          residual.zero();
          jacobian.zero();
          assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

          fillNumericalJacobian( residual.toViewConst(),
                                 residualOrig.toViewConst(),
                                 elemDofNumber[ei] + jc + 1,
                                 dRho,
                                 jacobianFD.toViewConstSizes() );
        }

        // Step 3: compute numerical derivatives wrt temperature

        if( isThermal )
        {
          arrayView1d< real64 > const temp =
            subRegion.getExtrinsicData< extrinsicMeshData::flow::temperature >();
          temp.move( LvArray::MemorySpace::host, false );

          solver.resetStateToBeginningOfStep( domain );

          temp.move( LvArray::MemorySpace::host, true );
          real64 const dT = perturbParameter * ( temp[ei] + perturbParameter );
          temp[ei] += dT;
#if defined(GEOSX_USE_CUDA)
          temp.move( LvArray::MemorySpace::cuda, false );
#endif

          // here, we make sure that rock internal energy is updated
          // in other words, a call to updateFluidState would not work
          solver.updateState( domain );

          residual.zero();
          jacobian.zero();
          assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

          fillNumericalJacobian( residual.toViewConst(),
                                 residualOrig.toViewConst(),
                                 elemDofNumber[ei] + numComp + 1,
                                 dT,
                                 jacobianFD.toViewConstSizes() );
        }

      }
    } );
  } );
}


} // namespace testing

} // namespace geosx

#endif //GEOSX_TESTCOMPFLOWUTILS_HPP
