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
 * @file PreconditionerTwoLevel.cpp
 */

#include "PreconditionerTwoLevel.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

namespace geosx
{

// --------------------------------------------
struct twoLevelStructuredMesh
{
  constexpr static integer nDim = 3;
  integer numDofPerNode;
  localIndex numCoarsePoints[nDim]{};
  integer refinementFactors[nDim]{};

  //twoLevelStructuredMesh()
  //{
  //  readFromFile( "twoLevelMesh.txt" );
  //}	  
  
  //void printInfo()
  //{
  //  std::cout << "twoLevelMesh Info" << std::endl;
  //  std::cout << "nDim: " << nDim << std::endl;
  //  std::cout << "numDofPerNode: " << numDofPerNode << std::endl;
  //  std::cout << "numCoarsePoints:";
  //  for( int iDim = 0; iDim < nDim; ++iDim )
  //  {
  //    std::cout << " " << numCoarsePoints[iDim];
  //  }
  //  std::cout << std::endl;
  //  std::cout << "refinementFactors:";
  //  for( int iDim = 0; iDim < nDim; ++iDim )
  //  {
  //    std::cout << " " << refinementFactors[iDim];
  //  }
  //  std::cout << std::endl;
  //  return;
  //}

  void readFromFile( string const & filename )
  {	  
    std::ifstream TLMeshFile( filename );
    GEOSX_ERROR_IF( !TLMeshFile, GEOSX_FMT( "Unable to open file for writing: {}", filename ) );
    std::string str;
    TLMeshFile >> str;
    TLMeshFile >> numDofPerNode;
    TLMeshFile >> str;
    TLMeshFile >> numCoarsePoints[0];
    TLMeshFile >> numCoarsePoints[1];
    TLMeshFile >> numCoarsePoints[2];
    TLMeshFile >> str;
    TLMeshFile >> refinementFactors[0];
    TLMeshFile >> refinementFactors[1];
    TLMeshFile >> refinementFactors[2];
    return;
  }	  

};   

namespace
{

GEOSX_HOST_DEVICE
void getIJKFromGlobalID( localIndex const & IND,
                         localIndex const & NI,
                         localIndex const & NJ,
                         localIndex const & NK,
                         localIndex & I,
                         localIndex & J,
                         localIndex & K )
{
  localIndex N2D = NI * NJ;
  I = IND % NI;
  J = ( IND % N2D ) / NI;
  K = IND / N2D;
  return;
}

GEOSX_HOST_DEVICE
localIndex getGlobalIDFromIJK( localIndex const & I,
                               localIndex const & J,
                               localIndex const & K,
                               localIndex const & NI,
                               localIndex const & NJ,
                               localIndex const & NK )
{
  return K * NI * NJ + J * NI + I;
}

void prolongVector( arrayView1d< real64 const > const & coarseValues,
		    arrayView1d< real64 > const & fineValues,
                    twoLevelStructuredMesh const & mesh )
{
  integer numDofPerNode = mesh.numDofPerNode;

  // Coarse mesh number of points
  integer const NI = mesh.numCoarsePoints[0];
  integer const NJ = mesh.numCoarsePoints[1];
  integer const NK = mesh.numCoarsePoints[2];
  localIndex const NNODES = NI*NK*NK;
  
  // Refinement factors in I, J and K
  integer const RI = mesh.refinementFactors[0];
  integer const RJ = mesh.refinementFactors[1];
  integer const RK = mesh.refinementFactors[2];

  // Compute fine mesh number of points
  localIndex const ni = ( NI - 1 ) * RI + 1;
  localIndex const nj = ( NJ - 1 ) * RJ + 1;
  localIndex const nk = ( NK - 1 ) * RK + 1;
  localIndex const nnodes = ni*nj*nk;

  // Check vector size correctness
  // ...
  

  // Compute prolonged vector entries
  forAll< parallelDevicePolicy<> >( nnodes, [=] GEOSX_HOST_DEVICE ( localIndex const iNode )
  {
    // Get fine node multi-index
    integer i;
    integer j;
    integer k;
    getIJKFromGlobalID( iNode, ni, nj, nk, i, j, k );

    // Compute I, J, K indeces of the bottom, south, west coarse vertex closest to the
    // processed fine scale node along with fine scale offsets
    integer const I = i / RI;
    integer const J = j / RJ;
    integer const K = k / RK;
      
    integer const OI = i % RI;
    integer const OJ = j % RJ;
    integer const OK = k % RK;

    // Determine whether fine scale node
    // - is a coarse VERTEX
    // - belongs to a coarse scale EDGE or FACE
    // - is INTERNAL
    integer DI = 0;
    integer DJ = 0;
    integer DK = 0;

    if( OI == 0 )
    {
      if( OJ == 0 )
      {
        if( OK == 0 )
        {
          // Fine node is a coarse vertex
          // DI, DJ, DK are correct 
        }
        else
        {
          // Fine scale node belongs to an edge in the K-direction
          DK = 1;
        }
      }
      else
      {
        if( OK == 0 )
        {
          // Fine scale node belongs to an edge in the J-direction
          DJ = 1;
        }
        else
        {
          // Fine scale node belongs to a face with fixed I index
          DJ = 1;
          DK = 1;
        }
      }
    }
    else
    {
      if( OJ == 0 )
      {
        if( OK == 0 )
        {
          //  Fine scale node belongs to an edge in the I-direction
          DI = 1;
        }
        else
        {
          // Fine scale node belongs to a face with fixed J index
          DI = 1;
          DK = 1;
        }
      }
      else
      {
        if( OK == 0 )
        {
          // Fine scale node belongs to a face with fixed K index
          DI = 1;
          DJ = 1;
        }
        else
        {
          // Fine scale node is INTERNAL
          DI = 1;
          DJ = 1;
          DK = 1;
        }
      }
    }
    // Interpolate
    for( integer iDof = 0; iDof < numDofPerNode; ++iDof )
    {
      real64 value = 0.0;
      for( integer lk = 0; lk <= DK ; ++lk )
      {
        for( integer lj = 0; lj <= DJ ; ++lj )
        {
          for( integer li = 0; li <= DI ; ++li )
          {
#if GEOSX_ENABLE_CUDA
            localIndex const ind = getGlobalIDFromIJK( I+li, J+lj, K+lk, NI, NJ, NK ) + iDof * NNODES;
#else
            localIndex const ind = getGlobalIDFromIJK( I+li, J+lj, K+lk, NI, NJ, NK ) * numDofPerNode + iDof;
#endif
            real64 const wI = DI > 0
                            ? 1.0 - (real64)li + pow(-1.0, 1+li) * (real64)OI / RI
                            : 1;
            real64 const wJ = DJ > 0
                            ? 1.0 - (real64)lj + pow(-1.0, 1+lj) * (real64)OJ / RJ
                            : 1;
            real64 const wK = DK > 0
                            ? 1.0 - (real64)lk + pow(-1.0, 1+lk) * (real64)OK / RK
                            : 1;
            value += wI * wJ * wK * coarseValues[ ind ];
	    ////////////////
	    //GEOSX_LOG_RANK_VAR( iNode );
	    //GEOSX_LOG_RANK_VAR( wI );
	    //GEOSX_LOG_RANK_VAR( wJ );
	    //GEOSX_LOG_RANK_VAR( wK );
	    //GEOSX_LOG_RANK_VAR( ind );
	    //GEOSX_LOG_RANK_VAR( coarseValues[ ind ] );
	    //GEOSX_LOG_RANK_VAR( wI * wJ * wK * coarseValues[ ind ] );
	    //GEOSX_LOG_RANK_VAR( value );
	    ////////////////
          }
        }
      }
#if GEOSX_ENABLE_CUDA
      fineValues[ iNode + iDof * nnodes ] = value;
#else
      fineValues[ iNode*numDofPerNode + iDof ] = value;
#endif
    }
  } );

  
  return;
}

void restrictVector( arrayView1d< real64 const > const & fineValues,
		     arrayView1d< real64 > const & coarseValues,
                     twoLevelStructuredMesh const & mesh )
{
  integer numDofPerNode = mesh.numDofPerNode;

  // Coarse mesh number of points
  integer const NI = mesh.numCoarsePoints[0];
  integer const NJ = mesh.numCoarsePoints[1];
  integer const NK = mesh.numCoarsePoints[2];
  localIndex const NNODES = NI*NK*NK;
  
  // Refinement factors in I, J and K
  integer const RI = mesh.refinementFactors[0];
  integer const RJ = mesh.refinementFactors[1];
  integer const RK = mesh.refinementFactors[2];

  // Compute fine mesh number of points
  localIndex const ni = ( NI - 1 ) * RI + 1;
  localIndex const nj = ( NJ - 1 ) * RJ + 1;
  localIndex const nk = ( NK - 1 ) * RK + 1;
  localIndex const nnodes = ni*nj*nk;

  // Check vector size correctness
  // ...

  // Compute prolonged vector entries
  //for( localIndex iNode = 0; iNode < NNODES; ++iNode )
  forAll< parallelDevicePolicy<> >( NNODES, [=] GEOSX_HOST_DEVICE ( localIndex const iNode )
  {
    // Get coarse node multi-index and corresponding fine node multi-index
    localIndex I;
    localIndex J;
    localIndex K;
    getIJKFromGlobalID( iNode, NI, NJ, NK, I, J, K );
    // std::cout << iNode << ": " << I << ", " << J << ", " << K << std::endl;
    localIndex const i = I * RI;
    localIndex const j = J * RJ;
    localIndex const k = K * RK;

    // Compute support region for basis function associated to the coarse node 
    localIndex const diMin = I > 0 ? -RI + 1 : 0;
    localIndex const diMax = I < NI - 1 ? RI - 1 : 0;
    localIndex const djMin = J > 0 ? -RJ + 1 : 0;
    localIndex const djMax = J < NJ - 1 ? RJ - 1 : 0;
    localIndex const dkMin = K > 0 ? -RK + 1 : 0;
    localIndex const dkMax = K < NK - 1 ? RK - 1 : 0;

    // Populate restricted vector
    for( integer iDof = 0; iDof < numDofPerNode; ++iDof )
    {
      real64 value = 0.0;
      for( integer dk = dkMin; dk <= dkMax; ++dk )
      {
        for( integer dj = djMin; dj <= djMax; ++dj )
        {
          for( integer di = diMin; di <= diMax; ++di )
          {
#if GEOSX_ENABLE_CUDA
            localIndex const ind = getGlobalIDFromIJK( i+di, j+dj, k+dk, ni, nj, nk ) + iDof * nnodes;
#else
            localIndex const ind = getGlobalIDFromIJK( i+di, j+dj, k+dk, ni, nj, nk ) * numDofPerNode + iDof;
#endif
            real64 const wI = di > 0
                            ? 1.0 - (real64)di / RI
                            : 1.0 + (real64)di / RI;
            real64 const wJ = dj > 0
                            ? 1.0 - (real64)dj / RJ
                            : 1.0 + (real64)dj / RJ;
            real64 const wK = dk > 0
                            ? 1.0 - (real64)dk / RK
                            : 1.0 + (real64)dk / RK;
            value += wI * wJ * wK * fineValues[ ind ];
          }
        }
      }
#if GEOSX_ENABLE_CUDA
      coarseValues[ iNode + iDof * NNODES ] = value;
#else
      coarseValues[ iNode*numDofPerNode + iDof ] = value;
#endif
    }
  } );

  return;
}

}
// --------------------------------------------

template< typename LAI >
PreconditionerTwoLevel< LAI >::
PreconditionerTwoLevel ( Matrix const & mat,
		         Matrix const & matCoarse )
  : Base()
{

  GEOSX_LAI_ASSERT( mat.ready() );
  GEOSX_LAI_ASSERT_MSG( mat.numLocalRows() == mat.numLocalCols(), "Fine matrix must be square" );
  this->setup( mat );

  GEOSX_LAI_ASSERT( matCoarse.ready() );
  GEOSX_LAI_ASSERT_MSG( matCoarse.numLocalRows() == matCoarse.numLocalCols(), "Coarse matrix must be square" );
  m_coarseMat = &matCoarse;

  m_diagInv.create( mat.numLocalRows(), mat.comm() );
  mat.extractDiagonal( m_diagInv );
  m_diagInv.reciprocal();

  LinearSolverParameters params; 
  params.isSymmetric = true;

  // params.solverType = LinearSolverParameters::SolverType::direct;
  // params.direct.parallel = 0;

  params.preconditionerType = LinearSolverParameters::PreconditionerType::amg;
  params.amg.smootherType = LinearSolverParameters::AMG::SmootherType::l1jacobi;
  // params.logLevel = 1;

  //if( params.solverType == LinearSolverParameters::SolverType::direct )
  //{
  //  GEOSX_LOG_RANK("Direct coarse solver");
  //  m_coarseSolver = LAI::createSolver( params );
  //}
  //else
  //{
    GEOSX_LOG_RANK("AMG coarse solver");
    m_coarseSolver = LAI::createPreconditioner( params );
  //}

  m_coarseSolver->setup( matCoarse );
  
  m_mesh = std::make_unique< twoLevelStructuredMesh >();
  m_mesh->readFromFile( "twoLevelMesh.txt" );

  // Create temporary vectors
  fineTmp.create( mat.numLocalCols(), mat.comm() );
  fineRhs.create( mat.numLocalCols(), mat.comm() );
  coarseTmp.create( m_coarseMat->numLocalCols(), m_coarseMat->comm() );
  coarseRhs.create( m_coarseMat->numLocalCols(), m_coarseMat->comm() );
}

template< typename LAI >
PreconditionerTwoLevel< LAI >::~PreconditionerTwoLevel() = default;

template< typename LAI >
void PreconditionerTwoLevel< LAI >::apply( Vector const & src,
                                           Vector & dst ) const
{
  // V-cycle with single smoother (Jacobi) smoother application
  
  // ... Presmoothing with dst initialized to zero 
  dst.zero();
  fineRhs.copy( src );
  m_diagInv.pointwiseProduct( fineRhs, fineTmp );
  dst.axpy( 1.0, fineTmp );
  this->matrix().residual( fineTmp, fineRhs, fineRhs );

  // ... Coarse-grid correction
  arrayView1d< real64 > const values_coarseRhs = coarseRhs.open();
  restrictVector( fineRhs.values( ), values_coarseRhs, *m_mesh );
  coarseRhs.close();
  m_coarseSolver->apply( coarseRhs, coarseTmp );

  arrayView1d< real64 > const values_fineTmp = fineTmp.open();
  prolongVector( coarseTmp.values( ), values_fineTmp, *m_mesh );
  fineTmp.close();

  dst.axpy( 1.0, fineTmp );
  this->matrix().residual( fineTmp, fineRhs, fineRhs );

  // ... Postsmoothing
  m_diagInv.pointwiseProduct( fineRhs, fineTmp );
  dst.axpy( 1.0, fineTmp );
}

//template< typename LAI >
//void PreconditionerTwoLevel< LAI >::clear()
//{
//  Base::clear();
//  m_precond->clear();
//  m_matSC.reset();
//}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class PreconditionerTwoLevel< TrilinosInterface >;
#endif

#ifdef GEOSX_USE_HYPRE
template class PreconditionerTwoLevel< HypreInterface >;
#endif

#ifdef GEOSX_USE_PETSC
template class PreconditionerTwoLevel< PetscInterface >;
#endif

}
