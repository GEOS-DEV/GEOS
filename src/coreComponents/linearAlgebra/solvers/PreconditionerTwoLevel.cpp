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
namespace
{

struct twoLevelStructuredMesh
{
  integer nDim = 3;
  localIndex numCoarsePoints[3] = { 4, 3, 3 };
  integer refinementFactor[3] = { 4, 2, 3 };
  
  void printInfo()
  {
    std::cout << "twoLevelMesh Info" << std::endl;
    std::cout << "nDim: " << nDim << std::endl;
    std::cout << "numCoarsePoints:";
    for( int iDim = 0; iDim < nDim; ++iDim )
    {
      std::cout << " " << numCoarsePoints[iDim];
    }
    std::cout << std::endl;
    std::cout << "refinementFactor:";
    for( int iDim = 0; iDim < nDim; ++iDim )
    {
      std::cout << " " << refinementFactor[iDim];
    }
    std::cout << std::endl;
    return;
  }

};   

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
                    twoLevelStructuredMesh const & mesh,
                    integer const & numDofPerNode )
{
  // Coarse mesh number of points
  integer const NI = mesh.numCoarsePoints[0];
  integer const NJ = mesh.numCoarsePoints[1];
  integer const NK = mesh.numCoarsePoints[2];
  localIndex const NNODES = NI*NK*NK;
  
  // Refinement factors in I, J and K
  integer const RI = mesh.refinementFactor[0];
  integer const RJ = mesh.refinementFactor[1];
  integer const RK = mesh.refinementFactor[2];

  // Compute fine mesh number of points
  localIndex const ni = ( NI - 1 ) * RI + 1;
  localIndex const nj = ( NJ - 1 ) * RJ + 1;
  localIndex const nk = ( NK - 1 ) * RK + 1;
  localIndex const nnodes = ni*nj*nk;

  // Check vector size correctness
  // ...
  

  // Compute prolonged vector entries
  integer i;
  integer j;
  integer k;
  for( int iNode = 0; iNode < nnodes; ++iNode )
  {
    // Get fine node multi-index
    getIJKFromGlobalID( iNode, ni, nj, nk, i, j, k );
    // std::cout << iNode << ": " << i << ", " << j << ", " << k << std::endl;

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
            localIndex const ind = getGlobalIDFromIJK( I+li, J+lj, K+lk, NI, NJ, NK ) * numDofPerNode + iDof;
            real64 const wI = DI > 0
                            ? 1.0 - (real64)li + std::pow(-1.0, 1+li) * (real64)OI / RI
                            : 1;
            real64 const wJ = DJ > 0
                            ? 1.0 - (real64)lj + std::pow(-1.0, 1+lj) * (real64)OJ / RJ
                            : 1;
            real64 const wK = DK > 0
                            ? 1.0 - (real64)lk + std::pow(-1.0, 1+lk) * (real64)OK / RK
                            : 1;
            value += wI * wJ * wK * coarseValues[ ind ];
          }
        }
      }
      fineValues[ iNode*numDofPerNode + iDof ] = value;
    }
  }

  
  return;
}

void restrictVector( arrayView1d< real64 const > const & fineValues,
		     arrayView1d< real64 > const & coarseValues,
                     twoLevelStructuredMesh const & mesh,
                     integer const & numDofPerNode )
{
  // Coarse mesh number of points
  integer const NI = mesh.numCoarsePoints[0];
  integer const NJ = mesh.numCoarsePoints[1];
  integer const NK = mesh.numCoarsePoints[2];
  localIndex const NNODES = NI*NK*NK;
  
  // Refinement factors in I, J and K
  integer const RI = mesh.refinementFactor[0];
  integer const RJ = mesh.refinementFactor[1];
  integer const RK = mesh.refinementFactor[2];

  // Compute fine mesh number of points
  localIndex const ni = ( NI - 1 ) * RI + 1;
  localIndex const nj = ( NJ - 1 ) * RJ + 1;
  localIndex const nk = ( NK - 1 ) * RK + 1;
  localIndex const nnodes = ni*nj*nk;

  // Check vector size correctness
  // ...

  // Compute prolonged vector entries
  for( localIndex iNode = 0; iNode < NNODES; ++iNode )
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
            localIndex const ind = getGlobalIDFromIJK( i+di, j+dj, k+dk, ni, nj, nk ) * numDofPerNode + iDof;
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
      coarseValues[ iNode*numDofPerNode + iDof ] = value;
    }
  }

  return;
}

}
// --------------------------------------------

template< typename LAI >
PreconditionerTwoLevel< LAI >::
PreconditionerTwoLevel ()
  : Base()
{
  twoLevelStructuredMesh mesh;
  integer const numDofPerNode = 3;
  localIndex const localCoarseSize = mesh.numCoarsePoints[0] * mesh.numCoarsePoints[1] * mesh.numCoarsePoints[2];
  localIndex const coarseSpaceDim = localCoarseSize * numDofPerNode;
  localIndex const localFineSize = ( ( mesh.numCoarsePoints[0] - 1 ) * mesh.refinementFactor[0] + 1 ) 
		                 * ( ( mesh.numCoarsePoints[1] - 1 ) * mesh.refinementFactor[1] + 1 )
				 * ( ( mesh.numCoarsePoints[2] - 1 ) * mesh.refinementFactor[2] + 1 );
  localIndex const fineSpaceDim = localFineSize * numDofPerNode;

  Vector vH;
  Vector wH;
  Vector vh;
  MPI_Comm const comm = MPI_COMM_GEOSX;
  int const rank = MpiWrapper::commRank( comm );
  int const nproc = MpiWrapper::commSize( comm );
  
  vH.create( coarseSpaceDim, comm );
  wH.create( coarseSpaceDim, comm );
  vh.create( fineSpaceDim, comm );
  arrayView1d< real64 > const values = vH.open();
  values[   0 ] = 8.984169707017889e-01;
  values[   1 ] = 1.730556422508585e-01;
  values[   2 ] = 4.984730689116268e-01;
  values[   3 ] = 2.262703081123629e-02;
  values[   4 ] = 9.414463161042915e-01;
  values[   5 ] = 6.676273433857585e-01;
  values[   6 ] = 7.476695656612741e-01;
  values[   7 ] = 7.065466594578906e-01;
  values[   8 ] = 6.705204299312063e-01;
  values[   9 ] = 5.298601187793164e-01;
  values[  10 ] = 6.883938266853351e-01;
  values[  11 ] = 8.464315203708332e-01;
  values[  12 ] = 6.785590194065089e-01;
  values[  13 ] = 4.773825703942811e-01;
  values[  14 ] = 8.405881411943199e-01;
  values[  15 ] = 8.237660055988379e-01;
  values[  16 ] = 6.639589894287550e-01;
  values[  17 ] = 9.546640985209540e-01;
  values[  18 ] = 8.348865374630587e-01;
  values[  19 ] = 3.138753818623718e-01;
  values[  20 ] = 9.080192489606731e-02;
  values[  21 ] = 8.744201525171559e-01;
  values[  22 ] = 1.866423577745185e-01;
  values[  23 ] = 4.855050980254463e-01;
  values[  24 ] = 3.674526305826995e-01;
  values[  25 ] = 1.173122833641906e-01;
  values[  26 ] = 7.589345951251416e-01;
  values[  27 ] = 4.895759772029249e-01;
  values[  28 ] = 3.863219491556974e-01;
  values[  29 ] = 6.498751386303387e-01;
  values[  30 ] = 4.357528394476876e-01;
  values[  31 ] = 5.046460598360366e-01;
  values[  32 ] = 3.007770337177867e-01;
  values[  33 ] = 2.146779382915700e-01;
  values[  34 ] = 1.473342711458868e-01;
  values[  35 ] = 8.233130359881409e-01;
  values[  36 ] = 4.151743724253540e-01;
  values[  37 ] = 9.133220439323685e-01;
  values[  38 ] = 4.123744225844106e-01;
  values[  39 ] = 5.350365064480540e-01;
  values[  40 ] = 8.706237991792931e-03;
  values[  41 ] = 4.990207991399526e-01;
  values[  42 ] = 9.818468585627238e-02;
  values[  43 ] = 6.961308820215216e-01;
  values[  44 ] = 5.692237745925788e-02;
  values[  45 ] = 2.877581005668916e-01;
  values[  46 ] = 3.147101113953451e-01;
  values[  47 ] = 5.280159866782138e-01;
  values[  48 ] = 9.624630982308021e-01;
  values[  49 ] = 3.253391743050699e-01;
  values[  50 ] = 7.467690503925523e-01;
  values[  51 ] = 4.390037803917550e-01;
  values[  52 ] = 2.754445076951573e-02;
  values[  53 ] = 7.642763681461628e-02;
  values[  54 ] = 1.629611534086595e-01;
  values[  55 ] = 8.221273654961192e-01;
  values[  56 ] = 1.487300662275823e-01;
  values[  57 ] = 5.891046555546018e-01;
  values[  58 ] = 8.738749524224014e-01;
  values[  59 ] = 4.724670802736430e-02;
  values[  60 ] = 6.035232300905453e-01;
  values[  61 ] = 7.800933310093185e-01;
  values[  62 ] = 6.285970419375732e-01;
  values[  63 ] = 9.795804536483746e-01;
  values[  64 ] = 5.577203103561332e-01;
  values[  65 ] = 9.921687439673234e-01;
  values[  66 ] = 1.297723532738355e-01;
  values[  67 ] = 8.546265810377777e-01;
  values[  68 ] = 5.755380415817175e-01;
  values[  69 ] = 2.261509893455913e-01;
  values[  70 ] = 6.761157212459386e-01;
  values[  71 ] = 5.846388626177723e-01;
  values[  72 ] = 4.293384459064786e-01;
  values[  73 ] = 3.512763711339811e-01;
  values[  74 ] = 3.101741591591093e-01;
  values[  75 ] = 3.901417866500019e-01;
  values[  76 ] = 2.710763527614077e-01;
  values[  77 ] = 2.463981969637962e-01;
  values[  78 ] = 6.578893457789561e-01;
  values[  79 ] = 5.191700191790319e-01;
  values[  80 ] = 4.450980676269678e-02;
  values[  81 ] = 7.985932520668604e-01;
  values[  82 ] = 8.688837429561447e-01;
  values[  83 ] = 1.387128427344385e-01;
  values[  84 ] = 7.009368844277101e-01;
  values[  85 ] = 2.558834553019614e-01;
  values[  86 ] = 1.776722341122916e-02;
  values[  87 ] = 5.616519368822925e-01;
  values[  88 ] = 1.873775539704424e-01;
  values[  89 ] = 5.018605249428387e-02;
  values[  90 ] = 3.323915742354102e-01;
  values[  91 ] = 8.936256612673193e-01;
  values[  92 ] = 3.331755152079796e-01;
  values[  93 ] = 6.182523732470672e-01;
  values[  94 ] = 8.981638052968632e-01;
  values[  95 ] = 2.635205945271855e-01;
  values[  96 ] = 5.394916409379239e-01;
  values[  97 ] = 8.432502544620608e-01;
  values[  98 ] = 1.536310267204350e-01;
  values[  99 ] = 3.238441101768030e-01;
  values[ 100 ] = 1.558311737470012e-01;
  values[ 101 ] = 7.219896896661782e-02;
  values[ 102 ] = 3.502688218288149e-01;
  values[ 103 ] = 9.680155816687932e-01;
  values[ 104 ] = 5.266514038048760e-01;
  values[ 105 ] = 2.206796923324806e-01;
  values[ 106 ] = 8.348798827094354e-01;
  values[ 107 ] = 5.176335258800429e-01;
  vH.close();

  // Prolong vector
  arrayView1d< real64 const > const vH_input = vH.values();
  arrayView1d< real64 > const vh_output = vh.open();
  prolongVector( vH_input, vh_output, mesh, numDofPerNode );
  vh.close();

  vh.write("vh");

  arrayView1d< real64 const > const vh_input = vh.values();
  arrayView1d< real64 > const wH_output = wH.open();
  restrictVector( vh_input, wH_output, mesh, numDofPerNode );
  wH.close();

  wH.write("wH");


}

//template< typename LAI >
//void PreconditionerTwoLevel< LAI >::setup( Matrix const & mat )
//{
//  Base::setup( mat );
//  // TODO: if matrix structure hasn't changed, can just copy entries into existing m_matSC
//  mat.separateComponentFilter( m_matSC, m_numComp );
//  m_precond->setup( m_matSC );
//}
//
//template< typename LAI >
//void PreconditionerTwoLevel< LAI >::apply( Vector const & src,
//                                                    Vector & dst ) const
//{
//  m_precond->apply( src, dst );
//}
//
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
