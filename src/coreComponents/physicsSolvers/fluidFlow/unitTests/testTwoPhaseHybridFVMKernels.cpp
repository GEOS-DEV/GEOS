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

// Source includes
#include "managers/initialization.hpp"
#include "physicsSolvers/fluidFlow/TwoPhaseHybridFVMKernels.hpp"
#include "codingUtilities/UnitTestUtilities.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::TwoPhaseHybridFVMKernels;
using namespace geosx::HybridFVMInnerProduct;
using namespace geosx::testing;


static real64 constexpr relTol = 1e-5;

// helpers

void updatePhaseDensity( real64 const & refPres,
                         real64 const & elemPres,
                         real64 const & dElemPres,
                         array1d< real64 > & elemPhaseDens,
                         array1d< real64 > & dElemPhaseDens_dp )
{
  localIndex constexpr numPhases = TwoPhaseBase::NUM_PHASES;

  // we assume (very) compressible flow to catch wrong derivatives
  stackArray1d< real64, numPhases > phaseCompressibility( numPhases );
  phaseCompressibility( 0 ) = 1e-3;
  phaseCompressibility( 1 ) = 1e-4;
  stackArray1d< real64, numPhases > refPhaseDens( numPhases );
  refPhaseDens( 0 ) = 800;
  refPhaseDens( 0 ) = 1000;

  // hard-coded relationship between pressure and density
  for( localIndex ip = 0; ip < numPhases; ++ip )
  {
    elemPhaseDens( ip ) = refPhaseDens( ip ) * exp( phaseCompressibility( ip ) * ( elemPres + dElemPres - refPres ) );
    dElemPhaseDens_dp( ip ) = phaseCompressibility( ip ) * elemPhaseDens( ip );
  }
}

void updateUpwindedPhaseMobilities( globalIndex const elemDofNumber,
                                    array1d< real64 > const & elemPhaseSat,
                                    array1d< real64 > const & dElemPhaseSat,
                                    array1d< real64 > const & elemPhaseDens,
                                    array1d< real64 > const & dElemPhaseDens_dp,
                                    stackArray2d< real64, MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > & upwPhaseMobility,
                                    stackArray2d< real64, MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > & dUpwPhaseMobility_dp,
                                    stackArray2d< real64, MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > & dUpwPhaseMobility_dS,
                                    stackArray2d< globalIndex, MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > & upwDofNumber,
                                    stackArray1d< globalIndex, MAX_NUM_FACES > & neighborDofNumber )
{
  localIndex constexpr numPhases = TwoPhaseBase::NUM_PHASES;

  // we assume that viscosity is independent of pressure
  stackArray1d< real64, numPhases > elemPhaseVisc( numPhases );
  elemPhaseVisc( 0 ) = 0.001;
  elemPhaseVisc( 1 ) = 0.005;

  // we assume quadratic relperms
  stackArray1d< real64, numPhases > elemPhaseRelPerm( numPhases );
  stackArray1d< real64, numPhases > dElemPhaseRelPerm_dS( numPhases );
  for( localIndex ip = 0; ip < numPhases; ++ip )
  {
    real64 const newSat = elemPhaseSat( ip ) + dElemPhaseSat( ip );
    elemPhaseRelPerm( ip ) = newSat * newSat;
    dElemPhaseRelPerm_dS( ip ) = 2 * newSat;
  }
  dElemPhaseRelPerm_dS( 1 ) *= -1;

  stackArray1d< real64, numPhases > elemPhaseMobility( numPhases );
  stackArray1d< real64, numPhases > dElemPhaseMobility_dp( numPhases );
  stackArray1d< real64, numPhases > dElemPhaseMobility_dS( numPhases );
  for( localIndex ip = 0; ip < numPhases; ++ip )
  {
    real64 const relPermOverVisc = elemPhaseRelPerm( ip ) / elemPhaseVisc( ip );
    real64 const dRelPermOverVisc_dS = dElemPhaseRelPerm_dS( ip ) / elemPhaseVisc( ip );
    elemPhaseMobility( ip ) = elemPhaseDens( ip ) * relPermOverVisc;
    dElemPhaseMobility_dp( ip ) = dElemPhaseDens_dp( ip ) * relPermOverVisc;
    dElemPhaseMobility_dS( ip ) = elemPhaseDens( ip ) * dRelPermOverVisc_dS;
  }

  for( localIndex ifaceLoc = 0; ifaceLoc < upwPhaseMobility.size( 0 ); ++ifaceLoc )
  {
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      upwPhaseMobility( ifaceLoc, ip ) = elemPhaseMobility( ip );
      dUpwPhaseMobility_dp( ifaceLoc, ip ) = dElemPhaseMobility_dp( ip );
      dUpwPhaseMobility_dS( ifaceLoc, ip ) = dElemPhaseMobility_dS( ip );
      upwDofNumber( ifaceLoc, ip ) = elemDofNumber;
    }
    neighborDofNumber( ifaceLoc ) = elemDofNumber;
  }
}


// in this function, we set up a 1-element problem in a tetrahedron
// the QTPFA transmissibility matrix comes from testHybridFVMInnerProducts
void setupProblemForTetra( localIndex & numFacesInElem,
                           array1d< localIndex > & elemToFaces,
                           array2d< real64 > & facePotential,
                           array2d< real64 > & dFacePotential,
                           array1d< real64 > & faceGravCoef,
                           real64 & refPres,
                           real64 & elemPres,
                           real64 & dElemPres,
                           real64 & elemGravCoef,
                           array1d< real64 > & elemPhaseSat,
                           array1d< real64 > & dElemPhaseSat,
                           array1d< real64 > & elemPhaseDens,
                           array1d< real64 > & dElemPhaseDens_dp,
                           stackArray2d< real64, MAX_NUM_FACES *MAX_NUM_FACES > & transMatrix )
{
  localIndex constexpr numPhases = TwoPhaseBase::NUM_PHASES;

  numFacesInElem = 4;

  // we assume that we are in a tetrahedron

  facePotential.resize( numFacesInElem, numPhases );
  dFacePotential.resize( numFacesInElem, numPhases );
  facePotential( 0, 0 ) = 1e5; dFacePotential( 0, 0 ) = 0;
  facePotential( 0, 1 ) = 2e5; dFacePotential( 0, 1 ) = 0;
  facePotential( 1, 0 ) = 3e5; dFacePotential( 1, 0 ) = 0;
  facePotential( 1, 1 ) = 1e5; dFacePotential( 1, 1 ) = 0;
  facePotential( 2, 0 ) = 6e5; dFacePotential( 2, 0 ) = 0;
  facePotential( 2, 1 ) = 1e4; dFacePotential( 2, 1 ) = 0;
  facePotential( 3, 0 ) = 3e4; dFacePotential( 3, 0 ) = 0;
  facePotential( 3, 1 ) = 4e3; dFacePotential( 3, 1 ) = 0;

  elemGravCoef = 5e1;
  faceGravCoef.resize( numFacesInElem );
  faceGravCoef( 0 ) = 1e1;
  faceGravCoef( 1 ) = 1e2;
  faceGravCoef( 2 ) = 3e1;
  faceGravCoef( 3 ) = 7e1;

  elemToFaces.resize( numFacesInElem );
  elemToFaces( 0 ) = 0;
  elemToFaces( 1 ) = 1;
  elemToFaces( 2 ) = 2;
  elemToFaces( 3 ) = 3;

  elemPres  = 1.2e5;
  refPres   = elemPres;
  dElemPres = 0;
  elemPhaseSat.resize( numPhases );
  dElemPhaseSat.resize( numPhases );
  elemPhaseSat( 0 ) = 0.4;
  elemPhaseSat( 1 ) = 0.6;
  dElemPhaseSat( 0 ) = 0;
  dElemPhaseSat( 1 ) = 0;
  elemPhaseDens.resize( numPhases );
  dElemPhaseDens_dp.resize( numPhases );
  updatePhaseDensity( refPres, elemPres, dElemPres, elemPhaseDens, dElemPhaseDens_dp );

  // the transmissibility matrix comes from the HybridFVMInnerProduct unit tests
  transMatrix.resize( numFacesInElem, numFacesInElem );
  transMatrix( 0, 0 ) =  5.25e-12;
  transMatrix( 0, 1 ) =  3.75e-12;
  transMatrix( 0, 2 ) =  2.25e-12;
  transMatrix( 0, 3 ) =  3.75e-12;

  transMatrix( 1, 0 ) =  3.75e-12;
  transMatrix( 1, 1 ) = 12.75e-12;
  transMatrix( 1, 2 ) = -5.25e-12;
  transMatrix( 1, 3 ) =  3.75e-12;

  transMatrix( 2, 0 ) =  2.25e-12;
  transMatrix( 2, 1 ) = -5.25e-12;
  transMatrix( 2, 2 ) = 18.75e-12;
  transMatrix( 2, 3 ) = -0.75e-12;

  transMatrix( 3, 0 ) =  3.75e-12;
  transMatrix( 3, 1 ) =  3.75e-12;
  transMatrix( 3, 2 ) = -0.75e-12;
  transMatrix( 3, 3 ) =  8.25e-12;

}


void setupMatrixAndRhsForTetra( localIndex const numFacesInElem,
                                globalIndex & elemDofNumber,
                                array1d< globalIndex > & faceDofNumber,
                                ParallelMatrix & matrix,
                                ParallelMatrix & matrixPerturb,
                                ParallelMatrix & matrixFD,
                                ParallelVector & rhs,
                                ParallelVector & rhsPerturb )
{
  localIndex constexpr numPhases = TwoPhaseBase::NUM_PHASES;

  elemDofNumber = 0;
  faceDofNumber.resize( numFacesInElem );
  faceDofNumber( 0 ) = 2;
  faceDofNumber( 1 ) = 4;
  faceDofNumber( 2 ) = 6;
  faceDofNumber( 3 ) = 8;

  matrix.createWithGlobalSize( numPhases*(numFacesInElem+1), numPhases*(numFacesInElem+1), MPI_COMM_GEOSX );
  matrixPerturb.createWithGlobalSize( numPhases*(numFacesInElem+1), numPhases*(numFacesInElem+1), MPI_COMM_GEOSX );
  matrixFD.createWithGlobalSize( numPhases*(numFacesInElem+1), numPhases*(numFacesInElem+1), MPI_COMM_GEOSX );
  rhs.createWithGlobalSize( numPhases*(numFacesInElem+1), MPI_COMM_GEOSX );
  rhsPerturb.createWithGlobalSize( numPhases*(numFacesInElem+1), MPI_COMM_GEOSX );

  matrix.open();
  matrixPerturb.open();
  matrixFD.open();

  // set up a full matrix
  for( globalIndex i = 0; i < numPhases*(numFacesInElem+1); ++i )
  {
    for( globalIndex j = 0; j < numPhases*(numFacesInElem+1); ++j )
    {
      matrix.insert( i, j, 1 );
      matrixPerturb.insert( i, j, 1 );
      matrixFD.insert( i, j, 1 );
    }
  }

  matrix.close();
  matrixPerturb.close();
  matrixFD.close();

}

TEST( TwoPhaseHybridFVMKernels, assembleConstraints )
{
  localIndex constexpr numPhases = TwoPhaseBase::NUM_PHASES;

  real64 const perturbParameter = sqrt( std::numeric_limits< real64 >::epsilon() );

  //////////////////////
  // 1) problem setup //
  //////////////////////

  localIndex numFacesInElem;
  array1d< localIndex > elemToFaces;
  array2d< real64 > facePotential;
  array2d< real64 > dFacePotential;
  array1d< real64 > faceGravCoef;
  real64 refPres;
  real64 elemPres;
  real64 dElemPres;
  real64 elemGravCoef;
  array1d< real64 > elemPhaseSat;
  array1d< real64 > dElemPhaseSat;
  array1d< real64 > elemPhaseDens;
  array1d< real64 > dElemPhaseDens_dp;
  stackArray2d< real64, MAX_NUM_FACES *MAX_NUM_FACES > transMatrix;

  setupProblemForTetra( numFacesInElem,
                        elemToFaces,
                        facePotential,
                        dFacePotential,
                        faceGravCoef,
                        refPres,
                        elemPres,
                        dElemPres,
                        elemGravCoef,
                        elemPhaseSat,
                        dElemPhaseSat,
                        elemPhaseDens,
                        dElemPhaseDens_dp,
                        transMatrix );

  globalIndex elemDofNumber;
  array1d< globalIndex > faceDofNumber;
  ParallelMatrix jacobian;
  ParallelMatrix jacobianPerturb;
  ParallelMatrix jacobianFD;
  ParallelVector rhs;
  ParallelVector rhsPerturb;

  setupMatrixAndRhsForTetra( numFacesInElem,
                             elemDofNumber,
                             faceDofNumber,
                             jacobian,
                             jacobianPerturb,
                             jacobianFD,
                             rhs,
                             rhsPerturb );

  real64 const * localRhs = rhs.extractLocalVector();
  real64 const * localRhsPerturb = rhsPerturb.extractLocalVector();

  stackArray2d< real64, MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > oneSidedVolFlux( numFacesInElem, numPhases );
  stackArray2d< real64, MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > dOneSidedVolFlux_dp( numFacesInElem, numPhases );
  stackArray2d< real64, MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > dOneSidedVolFlux_dS( numFacesInElem, numPhases );
  stackArray3d< real64, MAX_NUM_FACES *MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > dOneSidedVolFlux_dfp( numFacesInElem, numFacesInElem, numPhases );

  ///////////////////////////////////////
  // 2) Compute analytical derivatives //
  ///////////////////////////////////////

  FluxKernelHelper::ComputeOneSidedVolFluxes( facePotential,
                                              dFacePotential,
                                              faceGravCoef,
                                              elemToFaces,
                                              elemPres,
                                              dElemPres,
                                              elemGravCoef,
                                              elemPhaseDens,
                                              dElemPhaseDens_dp,
                                              transMatrix,
                                              oneSidedVolFlux,
                                              dOneSidedVolFlux_dp,
                                              dOneSidedVolFlux_dS,
                                              dOneSidedVolFlux_dfp );

  jacobianFD.zero();
  jacobianFD.open();
  jacobian.zero();
  jacobian.open();
  rhs.zero();
  rhs.open();

  FluxKernelHelper::AssembleConstraints( faceDofNumber,
                                         elemToFaces,
                                         elemDofNumber,
                                         oneSidedVolFlux,
                                         dOneSidedVolFlux_dp,
                                         dOneSidedVolFlux_dS,
                                         dOneSidedVolFlux_dfp,
                                         &jacobian,
                                         &rhs );

  jacobian.close();
  rhs.close();

  ///////////////////////////////////////////////////////////////////////////////////
  // 3) Compute finite-difference derivatives with respect to the element pressure //
  ///////////////////////////////////////////////////////////////////////////////////

  real64 const dElemPresPerturb = perturbParameter * (elemPres + perturbParameter);
  // we need to update density to account for the perturbation
  updatePhaseDensity( refPres, elemPres, dElemPresPerturb, elemPhaseDens, dElemPhaseDens_dp );

  FluxKernelHelper::ComputeOneSidedVolFluxes( facePotential,
                                              dFacePotential,
                                              faceGravCoef,
                                              elemToFaces,
                                              elemPres,
                                              dElemPresPerturb,
                                              elemGravCoef,
                                              elemPhaseDens,
                                              dElemPhaseDens_dp,
                                              transMatrix,
                                              oneSidedVolFlux,
                                              dOneSidedVolFlux_dp,
                                              dOneSidedVolFlux_dS,
                                              dOneSidedVolFlux_dfp );

  jacobianPerturb.zero();
  jacobianPerturb.open();
  rhsPerturb.zero();
  rhsPerturb.open();

  FluxKernelHelper::AssembleConstraints( faceDofNumber,
                                         elemToFaces,
                                         elemDofNumber,
                                         oneSidedVolFlux,
                                         dOneSidedVolFlux_dp,
                                         dOneSidedVolFlux_dS,
                                         dOneSidedVolFlux_dfp,
                                         &jacobianPerturb,
                                         &rhsPerturb );

  jacobianPerturb.close();
  rhsPerturb.close();

  for( localIndex lid = 0; lid < rhs.localSize(); ++lid )
  {
    real64 const dR_dp  = ( localRhsPerturb[lid] - localRhs[lid] ) / dElemPresPerturb;
    if( std::fabs( dR_dp ) > 0.0 )
    {
      globalIndex gid = rhs.getGlobalRowID( lid );
      jacobianFD.set( gid, elemDofNumber, dR_dp );
    }
  }

  //////////////////////////////////////////////////////////////////////////////////
  // 4) Compute finite-difference derivatives with respect to the face potentials //
  //////////////////////////////////////////////////////////////////////////////////

  // we need to revert the density to its initial value (it is independent of face potentials)
  updatePhaseDensity( refPres, elemPres, dElemPres, elemPhaseDens, dElemPhaseDens_dp );
  for( localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc )
  {
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      dFacePotential = 0;
      dFacePotential[ifaceLoc][ip] = perturbParameter * (facePotential[ifaceLoc][ip] + perturbParameter);

      FluxKernelHelper::ComputeOneSidedVolFluxes( facePotential,
                                                  dFacePotential,
                                                  faceGravCoef,
                                                  elemToFaces,
                                                  elemPres,
                                                  dElemPres,
                                                  elemGravCoef,
                                                  elemPhaseDens,
                                                  dElemPhaseDens_dp,
                                                  transMatrix,
                                                  oneSidedVolFlux,
                                                  dOneSidedVolFlux_dp,
                                                  dOneSidedVolFlux_dS,
                                                  dOneSidedVolFlux_dfp );

      jacobianPerturb.zero();
      jacobianPerturb.open();
      rhsPerturb.zero();
      rhsPerturb.open();

      FluxKernelHelper::AssembleConstraints( faceDofNumber,
                                             elemToFaces,
                                             elemDofNumber,
                                             oneSidedVolFlux,
                                             dOneSidedVolFlux_dp,
                                             dOneSidedVolFlux_dS,
                                             dOneSidedVolFlux_dfp,
                                             &jacobianPerturb,
                                             &rhsPerturb );

      jacobianPerturb.close();
      rhsPerturb.close();

      for( localIndex lid = 0; lid < rhs.localSize(); ++lid )
      {
        real64 const dR_dfp  = ( localRhsPerturb[lid] - localRhs[lid] ) / dFacePotential[ifaceLoc][ip];
        if( std::fabs( dR_dfp ) > 0.0 )
        {
          globalIndex gid = rhs.getGlobalRowID( lid );
          jacobianFD.set( gid, faceDofNumber[ifaceLoc] + ip, dR_dfp );
        }
      }
    }
  }

  jacobianFD.close();

  compareMatrices( jacobian, jacobianFD, relTol );
}


TEST( TwoPhaseHybridFVMKernels, assembleOneSidedMassFluxes )
{
  localIndex constexpr numPhases = TwoPhaseBase::NUM_PHASES;

  real64 const perturbParameter = sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const dt = 1e3;

  //////////////////////
  // 1) problem setup //
  //////////////////////

  localIndex numFacesInElem;
  array1d< localIndex > elemToFaces;
  array2d< real64 > facePotential;
  array2d< real64 > dFacePotential;
  array1d< real64 > faceGravCoef;
  real64 refPres;
  real64 elemPres;
  real64 dElemPres;
  real64 elemGravCoef;
  array1d< real64 > elemPhaseSat;
  array1d< real64 > dElemPhaseSat;
  array1d< real64 > elemPhaseDens;
  array1d< real64 > dElemPhaseDens_dp;
  array1d< localIndex > phaseToRow;
  stackArray2d< real64, MAX_NUM_FACES *MAX_NUM_FACES > transMatrix;

  phaseToRow.resize( numPhases );
  phaseToRow( 0 ) = 0;
  phaseToRow( 1 ) = 1;

  setupProblemForTetra( numFacesInElem,
                        elemToFaces,
                        facePotential,
                        dFacePotential,
                        faceGravCoef,
                        refPres,
                        elemPres,
                        dElemPres,
                        elemGravCoef,
                        elemPhaseSat,
                        dElemPhaseSat,
                        elemPhaseDens,
                        dElemPhaseDens_dp,
                        transMatrix );

  globalIndex elemDofNumber;
  array1d< globalIndex > faceDofNumber;
  ParallelMatrix jacobian;
  ParallelMatrix jacobianPerturb;
  ParallelMatrix jacobianFD;
  ParallelVector rhs;
  ParallelVector rhsPerturb;

  setupMatrixAndRhsForTetra( numFacesInElem,
                             elemDofNumber,
                             faceDofNumber,
                             jacobian,
                             jacobianPerturb,
                             jacobianFD,
                             rhs,
                             rhsPerturb );

  real64 const * localRhs = rhs.extractLocalVector();
  real64 const * localRhsPerturb = rhsPerturb.extractLocalVector();

  stackArray2d< real64, MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > upwPhaseMobility( numFacesInElem, numPhases );
  stackArray2d< real64, MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > dUpwPhaseMobility_dp( numFacesInElem, numPhases );
  stackArray2d< real64, MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > dUpwPhaseMobility_dS( numFacesInElem, numPhases );
  stackArray2d< globalIndex, MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > upwDofNumber( numFacesInElem, numPhases );
  stackArray1d< globalIndex, MAX_NUM_FACES > neighborDofNumber( numFacesInElem );

  updateUpwindedPhaseMobilities( elemDofNumber,
                                 elemPhaseSat,
                                 dElemPhaseSat,
                                 elemPhaseDens,
                                 dElemPhaseDens_dp,
                                 upwPhaseMobility,
                                 dUpwPhaseMobility_dp,
                                 dUpwPhaseMobility_dS,
                                 upwDofNumber,
                                 neighborDofNumber );

  stackArray2d< real64, MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > oneSidedVolFlux( numFacesInElem, numPhases );
  stackArray2d< real64, MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > dOneSidedVolFlux_dp( numFacesInElem, numPhases );
  stackArray2d< real64, MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > dOneSidedVolFlux_dS( numFacesInElem, numPhases );
  stackArray3d< real64, MAX_NUM_FACES *MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > dOneSidedVolFlux_dfp( numFacesInElem, numFacesInElem, numPhases );

  ///////////////////////////////////////
  // 2) Compute analytical derivatives //
  ///////////////////////////////////////

  FluxKernelHelper::ComputeOneSidedVolFluxes( facePotential,
                                              dFacePotential,
                                              faceGravCoef,
                                              elemToFaces,
                                              elemPres,
                                              dElemPres,
                                              elemGravCoef,
                                              elemPhaseDens,
                                              dElemPhaseDens_dp,
                                              transMatrix,
                                              oneSidedVolFlux,
                                              dOneSidedVolFlux_dp,
                                              dOneSidedVolFlux_dS,
                                              dOneSidedVolFlux_dfp );

  jacobianFD.zero();
  jacobianFD.open();
  jacobian.zero();
  jacobian.open();
  rhs.zero();
  rhs.open();

  FluxKernelHelper::AssembleOneSidedMassFluxes( dt,
                                                faceDofNumber,
                                                elemToFaces,
                                                elemDofNumber,
                                                phaseToRow,
                                                oneSidedVolFlux,
                                                dOneSidedVolFlux_dp,
                                                dOneSidedVolFlux_dS,
                                                dOneSidedVolFlux_dfp,
                                                upwPhaseMobility,
                                                dUpwPhaseMobility_dp,
                                                dUpwPhaseMobility_dS,
                                                upwDofNumber,
                                                neighborDofNumber,
                                                &jacobian,
                                                &rhs );

  jacobian.close();
  rhs.close();

  ///////////////////////////////////////////////////////////////////////////////////
  // 3) Compute finite-difference derivatives with respect to the element pressure //
  ///////////////////////////////////////////////////////////////////////////////////

  real64 const dElemPresPerturb = perturbParameter * (elemPres + perturbParameter);
  // we need to update density and mobility to account for the perturbation
  updatePhaseDensity( refPres, elemPres, dElemPresPerturb, elemPhaseDens, dElemPhaseDens_dp );
  updateUpwindedPhaseMobilities( elemDofNumber, elemPhaseSat, dElemPhaseSat, elemPhaseDens, dElemPhaseDens_dp,
                                 upwPhaseMobility, dUpwPhaseMobility_dp, dUpwPhaseMobility_dS, upwDofNumber, neighborDofNumber );

  FluxKernelHelper::ComputeOneSidedVolFluxes( facePotential,
                                              dFacePotential,
                                              faceGravCoef,
                                              elemToFaces,
                                              elemPres,
                                              dElemPresPerturb,
                                              elemGravCoef,
                                              elemPhaseDens,
                                              dElemPhaseDens_dp,
                                              transMatrix,
                                              oneSidedVolFlux,
                                              dOneSidedVolFlux_dp,
                                              dOneSidedVolFlux_dS,
                                              dOneSidedVolFlux_dfp );

  jacobianPerturb.zero();
  jacobianPerturb.open();
  rhsPerturb.zero();
  rhsPerturb.open();

  FluxKernelHelper::AssembleOneSidedMassFluxes( dt,
                                                faceDofNumber,
                                                elemToFaces,
                                                elemDofNumber,
                                                phaseToRow,
                                                oneSidedVolFlux,
                                                dOneSidedVolFlux_dp,
                                                dOneSidedVolFlux_dS,
                                                dOneSidedVolFlux_dfp,
                                                upwPhaseMobility,
                                                dUpwPhaseMobility_dp,
                                                dUpwPhaseMobility_dS,
                                                upwDofNumber,
                                                neighborDofNumber,
                                                &jacobianPerturb,
                                                &rhsPerturb );


  jacobianPerturb.close();
  rhsPerturb.close();

  for( localIndex lid = 0; lid < rhs.localSize(); ++lid )
  {
    real64 const dR_dp  = ( localRhsPerturb[lid] - localRhs[lid] ) / dElemPresPerturb;
    if( std::fabs( dR_dp ) > 0.0 )
    {
      globalIndex gid = rhs.getGlobalRowID( lid );
      jacobianFD.set( gid, elemDofNumber, dR_dp );
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////
  // 4) Compute finite-difference derivatives with respect to the element primary saturation //
  /////////////////////////////////////////////////////////////////////////////////////////////

  real64 const dElemSatPerturb = perturbParameter * (elemPhaseSat( 0 ) + perturbParameter);
  dElemPhaseSat( 0 ) = dElemSatPerturb;
  dElemPhaseSat( 1 ) = 1 - ( elemPhaseSat( 0 ) + dElemPhaseSat( 0 ) ) - elemPhaseSat( 1 );
  // we need to update density and mobility to account for the perturbation
  updatePhaseDensity( refPres, elemPres, dElemPres, elemPhaseDens, dElemPhaseDens_dp );
  updateUpwindedPhaseMobilities( elemDofNumber, elemPhaseSat, dElemPhaseSat, elemPhaseDens, dElemPhaseDens_dp,
                                 upwPhaseMobility, dUpwPhaseMobility_dp, dUpwPhaseMobility_dS, upwDofNumber, neighborDofNumber );

  FluxKernelHelper::ComputeOneSidedVolFluxes( facePotential,
                                              dFacePotential,
                                              faceGravCoef,
                                              elemToFaces,
                                              elemPres,
                                              dElemPres,
                                              elemGravCoef,
                                              elemPhaseDens,
                                              dElemPhaseDens_dp,
                                              transMatrix,
                                              oneSidedVolFlux,
                                              dOneSidedVolFlux_dp,
                                              dOneSidedVolFlux_dS,
                                              dOneSidedVolFlux_dfp );

  jacobianPerturb.zero();
  jacobianPerturb.open();
  rhsPerturb.zero();
  rhsPerturb.open();

  FluxKernelHelper::AssembleOneSidedMassFluxes( dt,
                                                faceDofNumber,
                                                elemToFaces,
                                                elemDofNumber,
                                                phaseToRow,
                                                oneSidedVolFlux,
                                                dOneSidedVolFlux_dp,
                                                dOneSidedVolFlux_dS,
                                                dOneSidedVolFlux_dfp,
                                                upwPhaseMobility,
                                                dUpwPhaseMobility_dp,
                                                dUpwPhaseMobility_dS,
                                                upwDofNumber,
                                                neighborDofNumber,
                                                &jacobianPerturb,
                                                &rhsPerturb );


  jacobianPerturb.close();
  rhsPerturb.close();

  for( localIndex lid = 0; lid < rhs.localSize(); ++lid )
  {
    real64 const dR_dS  = ( localRhsPerturb[lid] - localRhs[lid] ) / dElemPhaseSat( 0 );
    if( std::fabs( dR_dS ) > 0.0 )
    {
      globalIndex gid = rhs.getGlobalRowID( lid );
      jacobianFD.set( gid, elemDofNumber + 1, dR_dS );
    }
  }

  /////////////////////////////////////////////////////////////////////////////////
  // 5) Compute finite-difference derivatives with respect to the face pressures //
  /////////////////////////////////////////////////////////////////////////////////

  dElemPhaseSat( 0 ) = 0;
  dElemPhaseSat( 1 ) = 0;
  updatePhaseDensity( refPres, elemPres, dElemPres, elemPhaseDens, dElemPhaseDens_dp );
  updateUpwindedPhaseMobilities( elemDofNumber, elemPhaseSat, dElemPhaseSat, elemPhaseDens, dElemPhaseDens_dp,
                                 upwPhaseMobility, dUpwPhaseMobility_dp, dUpwPhaseMobility_dS, upwDofNumber, neighborDofNumber );
  for( localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc )
  {
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      dFacePotential = 0;
      dFacePotential[ifaceLoc][ip] = perturbParameter * (facePotential[ifaceLoc][ip] + perturbParameter);

      FluxKernelHelper::ComputeOneSidedVolFluxes( facePotential,
                                                  dFacePotential,
                                                  faceGravCoef,
                                                  elemToFaces,
                                                  elemPres,
                                                  dElemPres,
                                                  elemGravCoef,
                                                  elemPhaseDens,
                                                  dElemPhaseDens_dp,
                                                  transMatrix,
                                                  oneSidedVolFlux,
                                                  dOneSidedVolFlux_dp,
                                                  dOneSidedVolFlux_dS,
                                                  dOneSidedVolFlux_dfp );

      jacobianPerturb.zero();
      jacobianPerturb.open();
      rhsPerturb.zero();
      rhsPerturb.open();

      FluxKernelHelper::AssembleOneSidedMassFluxes( dt,
                                                    faceDofNumber,
                                                    elemToFaces,
                                                    elemDofNumber,
                                                    phaseToRow,
                                                    oneSidedVolFlux,
                                                    dOneSidedVolFlux_dp,
                                                    dOneSidedVolFlux_dS,
                                                    dOneSidedVolFlux_dfp,
                                                    upwPhaseMobility,
                                                    dUpwPhaseMobility_dp,
                                                    dUpwPhaseMobility_dS,
                                                    upwDofNumber,
                                                    neighborDofNumber,
                                                    &jacobianPerturb,
                                                    &rhsPerturb );


      jacobianPerturb.close();
      rhsPerturb.close();

      for( localIndex lid = 0; lid < rhs.localSize(); ++lid )
      {
        real64 const dR_dfp  = ( localRhsPerturb[lid] - localRhs[lid] ) / dFacePotential[ifaceLoc][ip];
        if( std::fabs( dR_dfp ) > 0.0 )
        {
          globalIndex gid = rhs.getGlobalRowID( lid );
          jacobianFD.set( gid, faceDofNumber[ifaceLoc] + ip, dR_dfp );
        }
      }
    }
  }

  jacobianFD.close();

  compareMatrices( jacobian, jacobianFD, relTol );

}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
