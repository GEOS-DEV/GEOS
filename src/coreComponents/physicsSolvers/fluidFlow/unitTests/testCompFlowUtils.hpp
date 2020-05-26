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

#ifndef GEOSX_TESTCOMPFLOWUTILS_HPP
#define GEOSX_TESTCOMPFLOWUTILS_HPP

#include "codingUtilities/UnitTestUtilities.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "meshUtilities/MeshManager.hpp"
#include "managers/ProblemManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFlow.hpp"

namespace geosx
{

namespace testing
{

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

void checkDerivative( arraySlice1d< real64 const > const & valueEps,
                      arraySlice1d< real64 const > const & value,
                      arraySlice1d< real64 const > const & deriv,
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

template< int DIM, typename ... Args >
void checkDerivative( ArraySlice< real64 const, DIM > const & valueEps,
                      ArraySlice< real64 const, DIM > const & value,
                      ArraySlice< real64 const, DIM > const & deriv,
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

template< int DIM, typename ... Args >
void checkDerivative( ArraySlice< real64 const, DIM > const & valueEps,
                      ArraySlice< real64 const, DIM > const & value,
                      ArraySlice< real64 const, DIM > const & deriv,
                      real64 const eps,
                      real64 const relTol,
                      string const & name,
                      string const & var,
                      arrayView1d< string const > const & labels,
                      Args ... label_lists )
{ return checkDerivative( valueEps, value, deriv, eps, relTol, DEFAULT_ABS_TOL, name, var, labels, label_lists ... ); }

// invert compositional derivative array layout to move innermost slice on the top
// (this is needed so we can use checkDerivative() to check derivative w.r.t. for each compositional var)
array1d< real64 > invertLayout( arraySlice1d< real64 const > const & input,
                                localIndex N )
{
  array1d< real64 > output( N );
  for( int i = 0; i < N; ++i )
  {
    output[i] = input[i];
  }

  return output;
}

array2d< real64 > invertLayout( arraySlice2d< real64 const > const & input,
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

array3d< real64 > invertLayout( arraySlice3d< real64 const > const & input,
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
  forAll< parallelHostPolicy >( residual.size(), [=]( localIndex const row )
  {
    real64 const dRdX = ( residual[row] - residualOrig[row] ) / eps;
    if( std::fabs( dRdX ) > 0.0 )
    {
      jacobian.addToRow< serialAtomic >( row, &dofIndex, &dRdX, 1 );
    }
  } );
}

template< typename LAMBDA >
void testNumericalJacobian( CompositionalMultiphaseFlow & solver,
                            DomainPartition & domain,
                            double perturbParameter,
                            double relTol,
                            LAMBDA assembleFunction )
{
  localIndex const NC = solver.numFluidComponents();

  CRSMatrix< real64, globalIndex > const & jacobian = solver.getLocalMatrix();
  array1d< real64 > const & residual = solver.getLocalRhs();
  DofManager const & dofManager = solver.getDofManager();

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  // assemble the analytical residual
  solver.ResetStateToBeginningOfStep( &domain );

  residual.setValues< parallelHostPolicy >( 0.0 );
  jacobian.setValues< parallelHostPolicy >( 0.0 );

  assembleFunction( jacobian.toViewConstSizes(), residual.toView() );
  residual.move( chai::CPU, false );

  // copy the analytical residual
  array1d< real64 > residualOrig( residual );

  // create the numerical jacobian
  CRSMatrix< real64, globalIndex > jacobianFD( jacobian );
  jacobianFD.setValues< parallelHostPolicy >( 0.0 );

  string const dofKey = dofManager.getKey( CompositionalMultiphaseFlow::viewKeyStruct::dofFieldString );

  solver.forTargetSubRegions( mesh, [&]( localIndex const,
                                         ElementSubRegionBase & subRegion )
  {
    arrayView1d< integer const > const & elemGhostRank =
      subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );

    arrayView1d< globalIndex const > const & dofNumber =
      subRegion.getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< real64 const > const pres =
      subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::pressureString );

    arrayView1d< real64 > const dPres =
      subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::deltaPressureString );

    arrayView2d< real64 const > const compDens =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::globalCompDensityString );
    compDens.move( chai::CPU );

    arrayView2d< real64 > const dCompDens =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::deltaGlobalCompDensityString );

    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      if( elemGhostRank[ei] >= 0 ) continue;

      real64 totalDensity = 0.0;
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        totalDensity += compDens[ei][ic];
      }

      {
        solver.ResetStateToBeginningOfStep( &domain );

        real64 const dP = perturbParameter * ( pres[ei] + perturbParameter );
        dPres.move( chai::CPU, true );
        dPres[ei] = dP;

        solver.forTargetSubRegions( mesh, [&]( localIndex const targetIndex2,
                                               ElementSubRegionBase & subRegion2 )
        {
          solver.UpdateState( subRegion2, targetIndex2 );
        } );

        residual.setValues< parallelHostPolicy >( 0.0 );
        jacobian.setValues< parallelHostPolicy >( 0.0 );
        assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

        fillNumericalJacobian( residual.toViewConst(),
                               residualOrig.toViewConst(),
                               dofNumber[ei],
                               dP,
                               jacobianFD.toViewConstSizes() );
      }

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        solver.ResetStateToBeginningOfStep( &domain );

        real64 const dRho = perturbParameter * totalDensity;
        dCompDens.move( chai::CPU, true );
        dCompDens[ei][jc] = dRho;

        solver.forTargetSubRegions( mesh, [&]( localIndex const targetIndex2,
                                               ElementSubRegionBase & subRegion2 )
        {
          solver.UpdateState( subRegion2, targetIndex2 );
        } );

        residual.setValues< parallelHostPolicy >( 0.0 );
        jacobian.setValues< parallelHostPolicy >( 0.0 );
        assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

        fillNumericalJacobian( residual.toViewConst(),
                               residualOrig.toViewConst(),
                               dofNumber[ei] + jc + 1,
                               dRho,
                               jacobianFD.toViewConstSizes() );
      }
    }
  } );

  // assemble the analytical jacobian
  solver.ResetStateToBeginningOfStep( &domain );

  residual.setValues< parallelHostPolicy >( 0.0 );
  jacobian.setValues< parallelHostPolicy >( 0.0 );
  assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

  compareLocalMatrices( jacobian.toViewConst(), jacobianFD.toViewConst(), relTol );
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

  int mpiSize = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
  dataRepository::Group * commandLine =
    problemManager.GetGroup< dataRepository::Group >( problemManager.groupKeys.commandLine );
  commandLine->registerWrapper< integer >( problemManager.viewKeys.xPartitionsOverride.Key() )->
    setApplyDefaultValue( mpiSize );

  xmlWrapper::xmlNode xmlProblemNode = xmlDocument.child( "Problem" );
  problemManager.InitializePythonInterpreter();
  problemManager.ProcessInputFileRecursive( xmlProblemNode );

  DomainPartition & domain  = *problemManager.getDomainPartition();

  constitutive::ConstitutiveManager & constitutiveManager = *domain.getConstitutiveManager();
  xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child( constitutiveManager.getName().c_str());
  constitutiveManager.ProcessInputFileRecursive( topLevelNode );

  MeshManager & meshManager = *problemManager.GetGroup< MeshManager >( problemManager.groupKeys.meshManager );
  meshManager.GenerateMeshLevels( &domain );

  ElementRegionManager & elementManager = *domain.getMeshBody( 0 )->getMeshLevel( 0 )->getElemManager();
  topLevelNode = xmlProblemNode.child( elementManager.getName().c_str());
  elementManager.ProcessInputFileRecursive( topLevelNode );

  problemManager.ProblemSetup();
}

} // namespace testing

} // namespace geosx

#endif //GEOSX_TESTCOMPFLOWUTILS_HPP
