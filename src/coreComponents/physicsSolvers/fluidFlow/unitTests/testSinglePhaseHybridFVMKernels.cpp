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
#include "physicsSolvers/fluidFlow/SinglePhaseHybridFVMKernels.hpp"
#include "codingUtilities/UnitTestUtilities.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::SinglePhaseHybridFVMKernels;
using namespace geosx::HybridFVMInnerProduct;
using namespace geosx::testing;


static real64 constexpr relTol = 1e-5;
static localIndex constexpr NF = 4; // we consider a tetrahedron in this file

// helpers

void updateDensity( real64 const & refPres,
                    real64 const & elemPres,
                    real64 const & dElemPres,
                    real64 & elemDens,
                    real64 & dElemDens_dp )
{
  // we assume (very) compressible flow to catch wrong derivatives
  real64 const compressibility = 0;//1e-3;
  real64 const refDens = 1000;

  // hard-coded relationship between pressure and density
  elemDens = refDens * exp( compressibility * ( elemPres + dElemPres - refPres ) );
  dElemDens_dp = compressibility * elemDens;
}

void updateUpwindedMobilities( globalIndex const elemDofNumber,
                               real64 const & elemDens,
                               real64 const & dElemDens_dp,
                               arraySlice1d< real64 > const & upwMobility,
                               arraySlice1d< real64 > const & dUpwMobility_dp,
                               arraySlice1d< globalIndex > const & upwDofNumber )
{
  // we assume that viscosity is independent of pressure
  real64 const elemVisc = 0.001;
  real64 elemMobility = elemDens / elemVisc;
  real64 dElemMobility_dp = dElemDens_dp / elemVisc;

  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    upwMobility( ifaceLoc ) = elemMobility;
    dUpwMobility_dp( ifaceLoc ) = dElemMobility_dp;
    upwDofNumber( ifaceLoc ) = elemDofNumber;
  }
}


// in this function, we set up a 1-element problem in a tetrahedron
// the QTPFA transmissibility matrix comes from testHybridFVMInnerProducts
void setupProblemForTetra( array1d< localIndex > & elemToFaces,
                           array1d< real64 > & facePres,
                           array1d< real64 > & dFacePres,
                           array1d< real64 > & faceGravCoef,
                           real64 & refPres,
                           real64 & elemPres,
                           real64 & dElemPres,
                           real64 & elemGravCoef,
                           real64 & elemDens,
                           real64 & dElemDens_dp,
                           arraySlice2d< real64 > const & transMatrix )
{
  facePres.resize( NF );
  dFacePres.resize( NF );
  facePres( 0 ) = 1e5; dFacePres( 0 ) = 0;
  facePres( 1 ) = 2e5; dFacePres( 1 ) = 0;
  facePres( 2 ) = 1e4; dFacePres( 2 ) = 0;
  facePres( 3 ) = 5e4; dFacePres( 3 ) = 0;

  elemGravCoef = 5e1;
  faceGravCoef.resize( NF );
  faceGravCoef( 0 ) = 1e1;
  faceGravCoef( 1 ) = 1e2;
  faceGravCoef( 2 ) = 3e1;
  faceGravCoef( 3 ) = 7e1;

  elemToFaces.resize( NF );
  elemToFaces( 0 ) = 0;
  elemToFaces( 1 ) = 1;
  elemToFaces( 2 ) = 2;
  elemToFaces( 3 ) = 3;

  elemPres  = 1.2e5;
  refPres   = elemPres;
  dElemPres = 0;
  updateDensity( refPres, elemPres, dElemPres, elemDens, dElemDens_dp );

  // the transmissibility matrix comes from the HybridFVMInnerProduct unit tests
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


void setupMatrixAndRhsForTetra( globalIndex & elemDofNumber,
                                array1d< globalIndex > & faceDofNumber,
                                ParallelMatrix & matrix,
                                ParallelMatrix & matrixPerturb,
                                ParallelMatrix & matrixFD,
                                ParallelVector & rhs,
                                ParallelVector & rhsPerturb )
{
  elemDofNumber = 0;
  faceDofNumber.resize( NF );
  faceDofNumber( 0 ) = 1;
  faceDofNumber( 1 ) = 2;
  faceDofNumber( 2 ) = 3;
  faceDofNumber( 3 ) = 4;

  matrix.createWithGlobalSize( NF+1, NF+1, MPI_COMM_GEOSX );
  matrixPerturb.createWithGlobalSize( NF+1, NF+1, MPI_COMM_GEOSX );
  matrixFD.createWithGlobalSize( NF+1, NF+1, MPI_COMM_GEOSX );
  rhs.createWithGlobalSize( NF+1, MPI_COMM_GEOSX );
  rhsPerturb.createWithGlobalSize( NF+1, MPI_COMM_GEOSX );

  matrix.open();
  matrixPerturb.open();
  matrixFD.open();

  // the matrices are full for the one-cell problem
  for( globalIndex i = 0; i < NF+1; ++i )
  {
    for( globalIndex j = 0; j < NF+1; ++j )
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

TEST( SinglePhaseHybridFVMKernels, assembleConstraints )
{
  real64 const perturbParameter = sqrt( std::numeric_limits< real64 >::epsilon() );

  //////////////////////
  // 1) problem setup //
  //////////////////////

  array1d< localIndex > elemToFaces;
  array1d< real64 > facePres;
  array1d< real64 > dFacePres;
  array1d< real64 > faceGravCoef;
  real64 refPres;
  real64 elemPres;
  real64 dElemPres;
  real64 elemGravCoef;
  real64 elemDens;
  real64 dElemDens_dp;
  stackArray2d< real64, NF *NF > transMatrix( NF, NF );

  setupProblemForTetra( elemToFaces,
                        facePres,
                        dFacePres,
                        faceGravCoef,
                        refPres,
                        elemPres,
                        dElemPres,
                        elemGravCoef,
                        elemDens,
                        dElemDens_dp,
                        transMatrix );

  globalIndex elemDofNumber;
  array1d< globalIndex > faceDofNumber;
  ParallelMatrix jacobian;
  ParallelMatrix jacobianPerturb;
  ParallelMatrix jacobianFD;
  ParallelVector rhs;
  ParallelVector rhsPerturb;

  setupMatrixAndRhsForTetra( elemDofNumber,
                             faceDofNumber,
                             jacobian,
                             jacobianPerturb,
                             jacobianFD,
                             rhs,
                             rhsPerturb );

  real64 const * localRhs = rhs.extractLocalVector();
  real64 const * localRhsPerturb = rhsPerturb.extractLocalVector();

  stackArray1d< real64, NF > oneSidedVolFlux( NF );
  stackArray1d< real64, NF > dOneSidedVolFlux_dp( NF );
  stackArray2d< real64, NF *NF > dOneSidedVolFlux_dfp( NF, NF );

  ///////////////////////////////////////
  // 2) Compute analytical derivatives //
  ///////////////////////////////////////

  AssemblerKernelHelper::ComputeOneSidedVolFluxes< NF >( facePres,
                                                         dFacePres,
                                                         faceGravCoef,
                                                         elemToFaces,
                                                         elemPres,
                                                         dElemPres,
                                                         elemGravCoef,
                                                         elemDens,
                                                         dElemDens_dp,
                                                         transMatrix,
                                                         oneSidedVolFlux,
                                                         dOneSidedVolFlux_dp,
                                                         dOneSidedVolFlux_dfp );

  jacobianFD.zero();
  jacobianFD.open();
  jacobian.zero();
  jacobian.open();
  rhs.zero();
  rhs.open();

  AssemblerKernelHelper::AssembleConstraints< NF >( faceDofNumber,
                                                    elemToFaces,
                                                    elemDofNumber,
                                                    oneSidedVolFlux,
                                                    dOneSidedVolFlux_dp,
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
  updateDensity( refPres, elemPres, dElemPresPerturb, elemDens, dElemDens_dp );

  oneSidedVolFlux = 0;
  dOneSidedVolFlux_dp = 0;
  dOneSidedVolFlux_dfp = 0;

  AssemblerKernelHelper::ComputeOneSidedVolFluxes< NF >( facePres,
                                                         dFacePres,
                                                         faceGravCoef,
                                                         elemToFaces,
                                                         elemPres,
                                                         dElemPresPerturb,
                                                         elemGravCoef,
                                                         elemDens,
                                                         dElemDens_dp,
                                                         transMatrix,
                                                         oneSidedVolFlux,
                                                         dOneSidedVolFlux_dp,
                                                         dOneSidedVolFlux_dfp );

  jacobianPerturb.zero();
  jacobianPerturb.open();
  rhsPerturb.zero();
  rhsPerturb.open();

  AssemblerKernelHelper::AssembleConstraints< NF >( faceDofNumber,
                                                    elemToFaces,
                                                    elemDofNumber,
                                                    oneSidedVolFlux,
                                                    dOneSidedVolFlux_dp,
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

  /////////////////////////////////////////////////////////////////////////////////
  // 4) Compute finite-difference derivatives with respect to the face pressures //
  /////////////////////////////////////////////////////////////////////////////////

  // we need to revert the density to its initial value (it is independent of face pressure)
  updateDensity( refPres, elemPres, dElemPres, elemDens, dElemDens_dp );
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    dFacePres = 0;
    dFacePres[ifaceLoc] = perturbParameter * (facePres[ifaceLoc] + perturbParameter);

    oneSidedVolFlux = 0;
    dOneSidedVolFlux_dp = 0;
    dOneSidedVolFlux_dfp = 0;

    AssemblerKernelHelper::ComputeOneSidedVolFluxes< NF >( facePres,
                                                           dFacePres,
                                                           faceGravCoef,
                                                           elemToFaces,
                                                           elemPres,
                                                           dElemPres,
                                                           elemGravCoef,
                                                           elemDens,
                                                           dElemDens_dp,
                                                           transMatrix,
                                                           oneSidedVolFlux,
                                                           dOneSidedVolFlux_dp,
                                                           dOneSidedVolFlux_dfp );

    jacobianPerturb.zero();
    jacobianPerturb.open();
    rhsPerturb.zero();
    rhsPerturb.open();

    AssemblerKernelHelper::AssembleConstraints< NF >( faceDofNumber,
                                                      elemToFaces,
                                                      elemDofNumber,
                                                      oneSidedVolFlux,
                                                      dOneSidedVolFlux_dp,
                                                      dOneSidedVolFlux_dfp,
                                                      &jacobianPerturb,
                                                      &rhsPerturb );

    jacobianPerturb.close();
    rhsPerturb.close();

    for( localIndex lid = 0; lid < rhs.localSize(); ++lid )
    {
      real64 const dR_dfp  = ( localRhsPerturb[lid] - localRhs[lid] ) / dFacePres[ifaceLoc];
      if( std::fabs( dR_dfp ) > 0.0 )
      {
        globalIndex gid = rhs.getGlobalRowID( lid );
        jacobianFD.set( gid, faceDofNumber[ifaceLoc], dR_dfp );
      }
    }
  }

  jacobianFD.close();

  compareMatrices( jacobian, jacobianFD, relTol );
}


TEST( SinglePhaseHybridFVMKernels, assembleOneSidedMassFluxes )
{
  real64 const perturbParameter = sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const dt = 1e3;

  //////////////////////
  // 1) problem setup //
  //////////////////////

  array1d< localIndex > elemToFaces;
  array1d< real64 > facePres;
  array1d< real64 > dFacePres;
  array1d< real64 > faceGravCoef;
  real64 refPres;
  real64 elemPres;
  real64 dElemPres;
  real64 elemGravCoef;
  real64 elemDens;
  real64 dElemDens_dp;
  stackArray2d< real64, NF *NF > transMatrix( NF, NF );

  setupProblemForTetra( elemToFaces,
                        facePres,
                        dFacePres,
                        faceGravCoef,
                        refPres,
                        elemPres,
                        dElemPres,
                        elemGravCoef,
                        elemDens,
                        dElemDens_dp,
                        transMatrix );

  globalIndex elemDofNumber;
  array1d< globalIndex > faceDofNumber;
  ParallelMatrix jacobian;
  ParallelMatrix jacobianPerturb;
  ParallelMatrix jacobianFD;
  ParallelVector rhs;
  ParallelVector rhsPerturb;

  setupMatrixAndRhsForTetra( elemDofNumber,
                             faceDofNumber,
                             jacobian,
                             jacobianPerturb,
                             jacobianFD,
                             rhs,
                             rhsPerturb );

  real64 const * localRhs = rhs.extractLocalVector();
  real64 const * localRhsPerturb = rhsPerturb.extractLocalVector();

  stackArray1d< real64, NF > upwMobility( NF );
  stackArray1d< real64, NF > dUpwMobility_dp( NF );
  stackArray1d< globalIndex, NF > upwDofNumber( NF );

  // no upwinding to do since we are in a one-cell problem
  updateUpwindedMobilities( elemDofNumber,
                            elemDens,
                            dElemDens_dp,
                            upwMobility,
                            dUpwMobility_dp,
                            upwDofNumber );

  stackArray1d< real64, NF > oneSidedVolFlux( NF );
  stackArray1d< real64, NF > dOneSidedVolFlux_dp( NF );
  stackArray2d< real64, NF *NF > dOneSidedVolFlux_dfp( NF, NF );

  ///////////////////////////////////////
  // 2) Compute analytical derivatives //
  ///////////////////////////////////////

  AssemblerKernelHelper::ComputeOneSidedVolFluxes< NF >( facePres,
                                                         dFacePres,
                                                         faceGravCoef,
                                                         elemToFaces,
                                                         elemPres,
                                                         dElemPres,
                                                         elemGravCoef,
                                                         elemDens,
                                                         dElemDens_dp,
                                                         transMatrix,
                                                         oneSidedVolFlux,
                                                         dOneSidedVolFlux_dp,
                                                         dOneSidedVolFlux_dfp );

  jacobianFD.zero();
  jacobianFD.open();
  jacobian.zero();
  jacobian.open();
  rhs.zero();
  rhs.open();

  AssemblerKernelHelper::AssembleOneSidedMassFluxes< NF >( dt,
                                                           faceDofNumber,
                                                           elemToFaces,
                                                           elemDofNumber,
                                                           oneSidedVolFlux,
                                                           dOneSidedVolFlux_dp,
                                                           dOneSidedVolFlux_dfp,
                                                           upwMobility,
                                                           dUpwMobility_dp,
                                                           upwDofNumber,
                                                           &jacobian,
                                                           &rhs );


  jacobian.close();
  rhs.close();

  ///////////////////////////////////////////////////////////////////////////////////
  // 3) Compute finite-difference derivatives with respect to the element pressure //
  ///////////////////////////////////////////////////////////////////////////////////

  real64 const dElemPresPerturb = perturbParameter * (elemPres + perturbParameter);
  // we need to update density and mobility to account for the perturbation
  updateDensity( refPres, elemPres, dElemPresPerturb, elemDens, dElemDens_dp );
  updateUpwindedMobilities( elemDofNumber, elemDens, dElemDens_dp, upwMobility, dUpwMobility_dp, upwDofNumber );

  oneSidedVolFlux = 0;
  dOneSidedVolFlux_dp = 0;
  dOneSidedVolFlux_dfp = 0;

  AssemblerKernelHelper::ComputeOneSidedVolFluxes< NF >( facePres,
                                                         dFacePres,
                                                         faceGravCoef,
                                                         elemToFaces,
                                                         elemPres,
                                                         dElemPresPerturb,
                                                         elemGravCoef,
                                                         elemDens,
                                                         dElemDens_dp,
                                                         transMatrix,
                                                         oneSidedVolFlux,
                                                         dOneSidedVolFlux_dp,
                                                         dOneSidedVolFlux_dfp );

  jacobianPerturb.zero();
  jacobianPerturb.open();
  rhsPerturb.zero();
  rhsPerturb.open();

  AssemblerKernelHelper::AssembleOneSidedMassFluxes< NF >( dt,
                                                           faceDofNumber,
                                                           elemToFaces,
                                                           elemDofNumber,
                                                           oneSidedVolFlux,
                                                           dOneSidedVolFlux_dp,
                                                           dOneSidedVolFlux_dfp,
                                                           upwMobility,
                                                           dUpwMobility_dp,
                                                           upwDofNumber,
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

  /////////////////////////////////////////////////////////////////////////////////
  // 4) Compute finite-difference derivatives with respect to the face pressures //
  /////////////////////////////////////////////////////////////////////////////////

  updateDensity( refPres, elemPres, dElemPres, elemDens, dElemDens_dp );
  updateUpwindedMobilities( elemDofNumber, elemDens, dElemDens_dp, upwMobility, dUpwMobility_dp, upwDofNumber );
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    dFacePres = 0;
    dFacePres[ifaceLoc] = perturbParameter * (facePres[ifaceLoc] + perturbParameter);

    oneSidedVolFlux = 0;
    dOneSidedVolFlux_dp = 0;
    dOneSidedVolFlux_dfp = 0;

    AssemblerKernelHelper::ComputeOneSidedVolFluxes< NF >( facePres,
                                                           dFacePres,
                                                           faceGravCoef,
                                                           elemToFaces,
                                                           elemPres,
                                                           dElemPres,
                                                           elemGravCoef,
                                                           elemDens,
                                                           dElemDens_dp,
                                                           transMatrix,
                                                           oneSidedVolFlux,
                                                           dOneSidedVolFlux_dp,
                                                           dOneSidedVolFlux_dfp );

    jacobianPerturb.zero();
    jacobianPerturb.open();
    rhsPerturb.zero();
    rhsPerturb.open();

    AssemblerKernelHelper::AssembleOneSidedMassFluxes< NF >( dt,
                                                             faceDofNumber,
                                                             elemToFaces,
                                                             elemDofNumber,
                                                             oneSidedVolFlux,
                                                             dOneSidedVolFlux_dp,
                                                             dOneSidedVolFlux_dfp,
                                                             upwMobility,
                                                             dUpwMobility_dp,
                                                             upwDofNumber,
                                                             &jacobianPerturb,
                                                             &rhsPerturb );


    jacobianPerturb.close();
    rhsPerturb.close();

    for( localIndex lid = 0; lid < rhs.localSize(); ++lid )
    {
      real64 const dR_dfp  = ( localRhsPerturb[lid] - localRhs[lid] ) / dFacePres[ifaceLoc];
      if( std::fabs( dR_dfp ) > 0.0 )
      {
        globalIndex gid = rhs.getGlobalRowID( lid );
        jacobianFD.set( gid, faceDofNumber[ifaceLoc], dR_dfp );
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
