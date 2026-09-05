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

/**
 * @file HypreRieszMFD.cpp
 */

#include "HypreRieszMFD.hpp"

#include "linearAlgebra/common/common.hpp"

#include "common/MpiWrapper.hpp"

#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_ls.h>
#include <_hypre_utilities.h>
#include <_hypre_parcsr_mv.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <iterator>
#include <numeric>
#include <string>
#include <utility>

namespace geos
{

namespace hypre
{

namespace
{

/// dof markers set by the physics solver
constexpr HYPRE_Int condensedFaceMarker = 0;
constexpr HYPRE_Int liveFaceMarker = 1;
constexpr HYPRE_Int pressureMarker = 2;

using RowMap = std::map< HYPRE_Int, HYPRE_Real >;

/// Owning container for the Riesz-map preconditioner. The hypre_Solver base is first so the
/// handle is also usable where hypre dispatches through the base struct (e.g. as an MGR solver).
struct RieszMFDData
{
  hypre_Solver base{};

  stdVector< HYPRE_Int > blockRow;      ///< row of each local dof within its block
  stdVector< HYPRE_Int > kindOf;        ///< marker of each local dof
  stdVector< HYPRE_Int > dofOfFlux0;    ///< block row -> local dof, condensed faces
  stdVector< HYPRE_Int > dofOfFlux1;    ///< block row -> local dof, live faces
  stdVector< HYPRE_Int > dofOfPres;     ///< block row -> local dof, pressures

  HYPRE_IJMatrix curlIJ{};              ///< discrete curl of the live-face sub-complex
  HYPRE_IJMatrix gradIJ{};              ///< discrete gradient of the live-face sub-complex
  HYPRE_IJVector coordIJ[3]{};          ///< active vertex coordinates

  stdVector< HYPRE_Int > isMfd;         ///< 1 for the pressure of an MFD cell, 0 for a TPFA cell
  stdVector< HYPRE_Real > normScale;    ///< (l_e/D)^2 of each pressure row: geometric factor of the L2 weight

  HYPRE_IJMatrix nIJ{};                 ///< M + B_M^T W_M^{-1} B_M on the live faces
  HYPRE_IJMatrix pIJ{};                 ///< pressure block: two-point couplings, L2 mass of the MFD cells, interface links
  HYPRE_IJMatrix b0tIJ{};               ///< F0, condensed-face rows x pressures
  stdVector< HYPRE_Real > invD0;        ///< inverse diagonal of the condensed-face closure rows
  stdVector< HYPRE_Int > pinned;        ///< 1 for a live flux dof out of the space
  stdVector< HYPRE_Real > invDiagPinned;  ///< exact diagonal solve of the pinned dofs
  stdVector< HYPRE_Real > rowScale;     ///< gamma / s_e: takes a conservation residual to the symmetrized system
  stdVector< HYPRE_Real > invWm;        ///< 1 / (W_M + C_M) on the MFD pressures, 0 elsewhere
  bool adsActive = false;               ///< false when no live dof is left for the flux block
  HYPRE_Solver ads{};
  HYPRE_Solver pcg{};                   ///< CG on the flux block preconditioned by ADS
  HYPRE_Solver amg{};

  HYPRE_IJVector r1{}, x1{}, rt{}, xt{}, xp{}, t0{};   ///< block-sized work vectors
};

/// PCG re-setup hook: ADS is set up once, before PCG
HYPRE_Int noSetup( HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector )
{
  return 0;
}

HYPRE_IJMatrix createIJFromRowMaps( stdVector< RowMap > const & rows,
                                    HYPRE_BigInt const numCols )
{
  HYPRE_BigInt const numRows = LvArray::integerConversion< HYPRE_BigInt >( rows.size() );
  HYPRE_IJMatrix ij{};
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixCreate( MPI_COMM_WORLD, 0, numRows - 1, 0, numCols - 1, &ij ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixSetObjectType( ij, HYPRE_PARCSR ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixInitialize( ij ) );
  stdVector< HYPRE_BigInt > cols;
  stdVector< HYPRE_Real > vals;
  for( HYPRE_BigInt i = 0; i < numRows; ++i )
  {
    RowMap const & row = rows[i];
    HYPRE_Int nnz = LvArray::integerConversion< HYPRE_Int >( row.size() );
    if( nnz == 0 )
    {
      continue;
    }
    cols.clear();
    vals.clear();
    for( auto const & [c, v] : row )
    {
      cols.push_back( c );
      vals.push_back( v );
    }
    GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixSetValues( ij, 1, &nnz, &i, cols.data(), vals.data() ) );
  }
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixAssemble( ij ) );
  return ij;
}

HYPRE_IJMatrix createIJFromCSR( arrayView1d< globalIndex const > const & rowPtr,
                                arrayView1d< globalIndex const > const & cols,
                                arrayView1d< real64 const > const & vals,
                                globalIndex const numCols )
{
  HYPRE_BigInt const numRows = LvArray::integerConversion< HYPRE_BigInt >( rowPtr.size() - 1 );
  HYPRE_IJMatrix ij{};
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixCreate( MPI_COMM_WORLD, 0, numRows - 1, 0, numCols - 1, &ij ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixSetObjectType( ij, HYPRE_PARCSR ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixInitialize( ij ) );
  stdVector< HYPRE_BigInt > colIdx;
  for( HYPRE_BigInt i = 0; i < numRows; ++i )
  {
    HYPRE_Int nnz = LvArray::integerConversion< HYPRE_Int >( rowPtr[i + 1] - rowPtr[i] );
    if( nnz == 0 )
    {
      continue;
    }
    colIdx.assign( nnz, 0 );
    for( HYPRE_Int k = 0; k < nnz; ++k )
    {
      colIdx[k] = LvArray::integerConversion< HYPRE_BigInt >( cols[rowPtr[i] + k] );
    }
    GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixSetValues( ij, 1, &nnz, &i, colIdx.data(), &vals[rowPtr[i]] ) );
  }
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixAssemble( ij ) );
  return ij;
}

HYPRE_IJVector createIJVector( HYPRE_BigInt const n, HYPRE_Real const * const vals )
{
  HYPRE_IJVector ij{};
  GEOS_LAI_CHECK_ERROR( HYPRE_IJVectorCreate( MPI_COMM_WORLD, 0, n - 1, &ij ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_IJVectorSetObjectType( ij, HYPRE_PARCSR ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_IJVectorInitialize( ij ) );
  if( n > 0 )
  {
    stdVector< HYPRE_BigInt > idx( n );
    std::iota( idx.begin(), idx.end(), HYPRE_BigInt( 0 ) );
    if( vals != nullptr )
    {
      GEOS_LAI_CHECK_ERROR( HYPRE_IJVectorSetValues( ij, n, idx.data(), vals ) );
    }
  }
  GEOS_LAI_CHECK_ERROR( HYPRE_IJVectorAssemble( ij ) );
  return ij;
}

/// ADS builds its vector interpolation from the linear part of the coordinate functions; a
/// projected coordinate system carries an offset that dwarfs the domain and turns that linear
/// part into a rounding error, so the coordinates are shifted to the origin
HYPRE_IJVector createShiftedCoordinates( arrayView1d< real64 const > const & vals )
{
  HYPRE_BigInt const n = LvArray::integerConversion< HYPRE_BigInt >( vals.size() );
  real64 const lo = n > 0 ? *std::min_element( vals.begin(), vals.end() ) : 0.0;
  stdVector< HYPRE_Real > shifted( n );
  for( HYPRE_BigInt i = 0; i < n; ++i )
  {
    shifted[i] = vals[i] - lo;
  }
  return createIJVector( n, shifted.data() );
}

hypre_ParCSRMatrix * parCSROf( HYPRE_IJMatrix const ij )
{
  hypre_ParCSRMatrix * mat{};
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixGetObject( ij, reinterpret_cast< void * * >( &mat ) ) );
  return mat;
}

hypre_ParVector * parVectorOf( HYPRE_IJVector const ij )
{
  hypre_ParVector * vec{};
  GEOS_LAI_CHECK_ERROR( HYPRE_IJVectorGetObject( ij, reinterpret_cast< void * * >( &vec ) ) );
  return vec;
}

HYPRE_Real * dataOf( HYPRE_IJVector const ij )
{
  return hypre_VectorData( hypre_ParVectorLocalVector( parVectorOf( ij ) ) );
}

void destroyIJ( HYPRE_IJMatrix & ij )
{
  if( ij )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixDestroy( ij ) );
    ij = nullptr;
  }
}

void destroyIJ( HYPRE_IJVector & ij )
{
  if( ij )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_IJVectorDestroy( ij ) );
    ij = nullptr;
  }
}

void destroySetupObjects( RieszMFDData & data )
{
  if( data.pcg )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRPCGDestroy( data.pcg ) );
    data.pcg = nullptr;
  }
  if( data.ads )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_ADSDestroy( data.ads ) );
    data.ads = nullptr;
  }
  if( data.amg )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGDestroy( data.amg ) );
    data.amg = nullptr;
  }
  destroyIJ( data.nIJ );
  destroyIJ( data.pIJ );
  destroyIJ( data.b0tIJ );
  destroyIJ( data.r1 );
  destroyIJ( data.x1 );
  destroyIJ( data.rt );
  destroyIJ( data.xt );
  destroyIJ( data.xp );
  destroyIJ( data.t0 );
}

/**
 * The assembled system, in the unknowns (q1 live fluxes, q0 condensed fluxes, p pressures), is
 *
 *   [ M      0    F  ] [ q1 ]      F  = -gamma B^T,  B the +-1 face-to-cell incidence of the live faces,
 *   [ 0      D0   F0 ] [ q0 ]      s_e B_e the conservation row of cell e (its own scale s_e),
 *   [ s B    0    C  ] [ p  ]      C  = s (C_acc + S_TT), S_TT the two-point Laplacian of the TPFA cells.
 *
 * Multiplying the conservation rows by gamma / s_e gives the symmetrized saddle point with
 * gamma B and its transpose, on which the Riesz map of the product norm
 *
 *   H(div; MFD region) x L2(MFD cells) x H1(TPFA cells)
 *
 * is P = diag( M + gamma^2 B_M^T W_M^{-1} B_M ,  W_M + C_M ,  C_T + S_TT ):
 * the flux block by CG under ADS, the MFD pressures by their diagonal, the TPFA pressures by one
 * AMG cycle. The interface faces are not weighted: they are bounded by the normal-trace inequality.
 * The condensed fluxes are recovered exactly from their closure rows.
 */
HYPRE_Int rieszSetup( HYPRE_Solver solver,
                      HYPRE_ParCSRMatrix A,
                      HYPRE_ParVector,
                      HYPRE_ParVector )
{
  RieszMFDData & data = *reinterpret_cast< RieszMFDData * >( solver );
  destroySetupObjects( data );

  hypre_ParCSRMatrix * const mat = reinterpret_cast< hypre_ParCSRMatrix * >( A );
  hypre_CSRMatrix * const diag = hypre_ParCSRMatrixDiag( mat );
  HYPRE_Int const n = hypre_CSRMatrixNumRows( diag );
  HYPRE_Int const * const ia = hypre_CSRMatrixI( diag );
  HYPRE_Int const * const ja = hypre_CSRMatrixJ( diag );
  HYPRE_Real const * const va = hypre_CSRMatrixData( diag );
  GEOS_ERROR_IF_NE_MSG( LvArray::integerConversion< std::size_t >( n ), data.kindOf.size(),
                        "Riesz MFD preconditioner: the dof markers do not match the system size" );

  HYPRE_Int const nf0 = LvArray::integerConversion< HYPRE_Int >( data.dofOfFlux0.size() );
  HYPRE_Int const nf1 = LvArray::integerConversion< HYPRE_Int >( data.dofOfFlux1.size() );
  HYPRE_Int const np = LvArray::integerConversion< HYPRE_Int >( data.dofOfPres.size() );

  // 1) read the blocks off the assembled system
  stdVector< RowMap > m1( nf1 ), c( np ), b0t( nf0 );
  stdVector< RowMap > b1( np );                         // conservation rows x live faces
  stdVector< stdVector< std::pair< HYPRE_Int, HYPRE_Real > > > b1Cols( nf1 );  // live face -> (cell, value)
  HYPRE_Real gamma = 0.0;                                // magnitude of the (1,2) entries: |F| = gamma
  data.invD0.assign( nf0, 0.0 );
  for( HYPRE_Int i = 0; i < n; ++i )
  {
    HYPRE_Int const ki = data.kindOf[i];
    HYPRE_Int const ri = data.blockRow[i];
    for( HYPRE_Int k = ia[i]; k < ia[i + 1]; ++k )
    {
      HYPRE_Int const j = ja[k];
      HYPRE_Int const kj = data.kindOf[j];
      HYPRE_Int const rj = data.blockRow[j];
      HYPRE_Real const v = va[k];
      if( ki == liveFaceMarker && kj == liveFaceMarker )
      {
        m1[ri][rj] += v;
      }
      else if( ki == pressureMarker && kj == liveFaceMarker )
      {
        b1[ri][rj] += v;
        b1Cols[rj].emplace_back( ri, v );
      }
      else if( ki == liveFaceMarker && kj == pressureMarker )
      {
        gamma = std::max( gamma, std::abs( v ) );
      }
      else if( ki == pressureMarker && kj == pressureMarker )
      {
        c[ri][rj] += v;
      }
      else if( ki == condensedFaceMarker && kj == pressureMarker )
      {
        b0t[ri][rj] += v;
      }
      else if( ki == condensedFaceMarker && kj == condensedFaceMarker && i == j )
      {
        data.invD0[ri] = v != 0.0 ? 1.0 / v : 0.0;
      }
    }
  }

  // 2) a live flux dof the divergence never touches (no-flow face) or whose row couples
  // negligibly against its own diagonal is out of the space: solved by its diagonal, outside the block
  HYPRE_Real bMaxAll = 0.0;
  for( HYPRE_Int e = 0; e < np; ++e )
  {
    for( auto const & [f, bf] : b1[e] )
    {
      bMaxAll = std::max( bMaxAll, std::abs( bf ) );
    }
  }
  data.pinned.assign( nf1, 0 );
  data.invDiagPinned.assign( nf1, 0.0 );
  HYPRE_Int numFree = 0;
  for( HYPRE_Int f = 0; f < nf1; ++f )
  {
    HYPRE_Real colMax = 0.0;
    for( auto const & [e, be] : b1Cols[f] )
    {
      colMax = std::max( colMax, std::abs( be ) );
    }
    HYPRE_Real offMax = colMax;
    for( auto const & [g, v] : m1[f] )
    {
      offMax = g != f ? std::max( offMax, std::abs( v ) ) : offMax;
    }
    auto const dIt = m1[f].find( f );
    HYPRE_Real const diagAbs = dIt != m1[f].end() ? std::abs( dIt->second ) : 0.0;
    if( colMax <= 1e-12 * bMaxAll || offMax <= 1e-10 * diagAbs )
    {
      data.pinned[f] = 1;
      data.invDiagPinned[f] = diagAbs > 0.0 ? 1.0 / dIt->second : 1.0;
      for( auto const & [g, v] : m1[f] )
      {
        if( g != f )
        {
          m1[g].erase( f );
        }
      }
      HYPRE_Real const diagVal = diagAbs > 0.0 ? dIt->second : 1.0;
      m1[f].clear();
      m1[f][f] = diagVal;
      b1Cols[f].clear();
    }
    else
    {
      ++numFree;
    }
  }
  for( HYPRE_Int e = 0; e < np; ++e )
  {
    for( auto it = b1[e].begin(); it != b1[e].end(); )
    {
      it = data.pinned[it->first] ? b1[e].erase( it ) : std::next( it );
    }
  }

  // 3) the scale s_e of each conservation row is the magnitude of its incidence entries; rows
  // without a live face take the median of the others (the scaling is per field)
  stdVector< HYPRE_Real > s( np, 0.0 ), known;
  for( HYPRE_Int e = 0; e < np; ++e )
  {
    for( auto const & [f, bf] : b1[e] )
    {
      s[e] = std::max( s[e], std::abs( bf ) );
    }
    if( s[e] > 0.0 )
    {
      known.push_back( s[e] );
    }
  }
  HYPRE_Real sMedian = gamma > 0.0 ? gamma : 1.0;
  if( !known.empty() )
  {
    std::nth_element( known.begin(), known.begin() + known.size() / 2, known.end() );
    sMedian = known[known.size() / 2];
  }
  if( gamma == 0.0 )
  {
    gamma = sMedian;
  }
  data.rowScale.assign( np, 0.0 );
  for( HYPRE_Int e = 0; e < np; ++e )
  {
    data.rowScale[e] = gamma / ( s[e] > 0.0 ? s[e] : sMedian );
  }

  // 4) MFD pressures: W_M = (l_e/D)^2 gamma^2 sum_f 1/M_ff, the L2 mass of the cell in the units of
  // the symmetrized operator, plus its accumulation
  stdVector< HYPRE_Real > invDiagM1( nf1, 0.0 );
  for( HYPRE_Int f = 0; f < nf1; ++f )
  {
    auto const it = m1[f].find( f );
    invDiagM1[f] = ( it != m1[f].end() && it->second != 0.0 ) ? 1.0 / it->second : 0.0;
  }
  // A cell belongs to the H1 block iff its two-point couplings connect it to a Dirichlet
  // contribution: a TPFA-flagged cell whose faces are all live has no coupling, and an island of
  // TPFA cells enclosed by MFD cells has a singular (pure Neumann) sub-block whose level is set by
  // the interface fluxes; both live in the L2 block with the saddle point
  stdVector< HYPRE_Int > root( np );
  std::iota( root.begin(), root.end(), 0 );
  auto const findRoot = [&]( HYPRE_Int e )
  {
    while( root[e] != e )
    {
      root[e] = root[root[e]];
      e = root[e];
    }
    return e;
  };
  for( HYPRE_Int e = 0; e < np; ++e )
  {
    for( auto const & [g, v] : c[e] )
    {
      if( g != e && v != 0.0 )
      {
        root[findRoot( e )] = findRoot( g );
      }
    }
  }
  std::map< HYPRE_Int, bool > componentGrounded;
  for( HYPRE_Int e = 0; e < np; ++e )
  {
    HYPRE_Real diagVal = 0.0, offSum = 0.0;
    bool coupled = false;
    for( auto const & [g, v] : c[e] )
    {
      diagVal += g == e ? v : 0.0;
      offSum += g != e ? std::abs( v ) : 0.0;
      coupled = coupled || ( g != e && v != 0.0 );
    }
    if( coupled )
    {
      bool & grounded = componentGrounded[findRoot( e )];
      grounded = grounded || diagVal > offSum * ( 1.0 + 1e-10 );
    }
  }
  data.invWm.assign( np, 0.0 );
  stdVector< HYPRE_Real > massOf( np, 0.0 );   // W_M of the L2 cells
  for( HYPRE_Int e = 0; e < np; ++e )
  {
    auto const cIt = componentGrounded.find( findRoot( e ) );
    bool const h1 = cIt != componentGrounded.end() && cIt->second;
    if( data.isMfd[e] || !h1 )
    {
      data.isMfd[e] = 1;
      HYPRE_Real w = 0.0;
      for( auto const & [f, bf] : b1[e] )
      {
        w += invDiagM1[f];
      }
      w *= data.normScale[e] * gamma * gamma;
      massOf[e] = w;
      auto const it = c[e].find( e );
      HYPRE_Real const acc = it != c[e].end() ? data.rowScale[e] * it->second : 0.0;
      data.invWm[e] = ( w + acc ) > 0.0 ? 1.0 / ( w + acc ) : 0.0;
    }
    else
    {
      data.isMfd[e] = 0;
    }
  }

  // 5) flux block N = M + gamma^2 B_M^T W_M^{-1} B_M, the divergence taken into MFD cells only
  for( HYPRE_Int e = 0; e < np; ++e )
  {
    if( !data.isMfd[e] || s[e] <= 0.0 )
    {
      continue;
    }
    HYPRE_Real const wInv = gamma * gamma * data.invWm[e] / ( s[e] * s[e] );
    for( auto const & [f, bf] : b1[e] )
    {
      for( auto const & [g, bg] : b1[e] )
      {
        m1[f][g] += bf * bg * wInv;
      }
    }
  }

  // 6) pressure block, one sparse SPD matrix over all cells: the two-point couplings wherever a
  // face was condensed (R^{1/2} C R^{1/2}, accumulation included), the L2 mass W_M of the L2 cells,
  // and the two-point energy gamma^2 / M_ff of each interface face, which weighs the mean pressure
  // of a patch like its neighbours instead of by its (h^2 smaller) L2 mass
  stdVector< RowMap > pp( np );
  for( HYPRE_Int e = 0; e < np; ++e )
  {
    for( auto const & [g, v] : c[e] )
    {
      pp[e][g] += std::sqrt( data.rowScale[e] * data.rowScale[g] ) * v;
    }
    pp[e][e] += massOf[e];
  }
  for( HYPRE_Int f = 0; f < nf1; ++f )
  {
    if( b1Cols[f].size() != 2 )
    {
      continue;
    }
    HYPRE_Int const e = b1Cols[f][0].first;
    HYPRE_Int const g = b1Cols[f][1].first;
    if( data.isMfd[e] == data.isMfd[g] )
    {
      continue;
    }
    HYPRE_Real const tau = gamma * gamma * invDiagM1[f];
    pp[e][e] += tau;
    pp[g][g] += tau;
    pp[e][g] -= tau;
    pp[g][e] -= tau;
  }

  data.nIJ = createIJFromRowMaps( m1, nf1 );
  data.pIJ = createIJFromRowMaps( pp, np );
  data.b0tIJ = createIJFromRowMaps( b0t, np );
  data.r1 = createIJVector( nf1, nullptr );
  data.x1 = createIJVector( nf1, nullptr );
  data.rt = createIJVector( np, nullptr );
  data.xt = createIJVector( np, nullptr );
  data.xp = createIJVector( np, nullptr );
  data.t0 = createIJVector( nf0, nullptr );

  // 7) CG under one ADS cycle on the flux block (the reference's ads-cg)
  data.adsActive = numFree > 0;
  if( data.adsActive )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_ADSCreate( &data.ads ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_ADSSetDiscreteCurl( data.ads, reinterpret_cast< HYPRE_ParCSRMatrix >( parCSROf( data.curlIJ ) ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_ADSSetDiscreteGradient( data.ads, reinterpret_cast< HYPRE_ParCSRMatrix >( parCSROf( data.gradIJ ) ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_ADSSetCoordinateVectors( data.ads,
                                                         reinterpret_cast< HYPRE_ParVector >( parVectorOf( data.coordIJ[0] ) ),
                                                         reinterpret_cast< HYPRE_ParVector >( parVectorOf( data.coordIJ[1] ) ),
                                                         reinterpret_cast< HYPRE_ParVector >( parVectorOf( data.coordIJ[2] ) ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_ADSSetCycleType( data.ads, 13 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_ADSSetMaxIter( data.ads, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_ADSSetTol( data.ads, 0.0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_ADSSetPrintLevel( data.ads, 0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_ADSSetAMGOptions( data.ads, 10, 1, 3, 0.25, 0, 0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_ADSSetAMSOptions( data.ads, 11, 10, 1, 3, 0.25, 0, 0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_ADSSetup( data.ads,
                                          reinterpret_cast< HYPRE_ParCSRMatrix >( parCSROf( data.nIJ ) ),
                                          reinterpret_cast< HYPRE_ParVector >( parVectorOf( data.r1 ) ),
                                          reinterpret_cast< HYPRE_ParVector >( parVectorOf( data.x1 ) ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRPCGCreate( MPI_COMM_WORLD, &data.pcg ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_PCGSetTol( data.pcg, 1e-2 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_PCGSetMaxIter( data.pcg, 50 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_PCGSetTwoNorm( data.pcg, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_PCGSetPrintLevel( data.pcg, 0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_PCGSetLogging( data.pcg, 0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRPCGSetPrecond( data.pcg,
                                                     reinterpret_cast< HYPRE_PtrToParSolverFcn >( HYPRE_ADSSolve ),
                                                     reinterpret_cast< HYPRE_PtrToParSolverFcn >( noSetup ),
                                                     data.ads ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRPCGSetup( data.pcg,
                                                reinterpret_cast< HYPRE_ParCSRMatrix >( parCSROf( data.nIJ ) ),
                                                reinterpret_cast< HYPRE_ParVector >( parVectorOf( data.r1 ) ),
                                                reinterpret_cast< HYPRE_ParVector >( parVectorOf( data.x1 ) ) ) );
  }

  // 8) one BoomerAMG cycle on the pressure block
  if( np > 0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &data.amg ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( data.amg, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( data.amg, 0.0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( data.amg, 0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCoarsenType( data.amg, 8 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpType( data.amg, 6 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetAggNumLevels( data.amg, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( data.amg, 8 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( data.amg, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetStrongThreshold( data.amg, 0.25 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetup( data.amg,
                                                reinterpret_cast< HYPRE_ParCSRMatrix >( parCSROf( data.pIJ ) ),
                                                reinterpret_cast< HYPRE_ParVector >( parVectorOf( data.rt ) ),
                                                reinterpret_cast< HYPRE_ParVector >( parVectorOf( data.xt ) ) ) );
  }

  hypre_SolverSetIsSetup( &data.base );
  return 0;
}

HYPRE_Int rieszSolve( HYPRE_Solver solver,
                      HYPRE_ParCSRMatrix,
                      HYPRE_ParVector f,
                      HYPRE_ParVector u )
{
  RieszMFDData & data = *reinterpret_cast< RieszMFDData * >( solver );
  HYPRE_Real const * const r = hypre_VectorData( hypre_ParVectorLocalVector( reinterpret_cast< hypre_ParVector * >( f ) ) );
  HYPRE_Real * const x = hypre_VectorData( hypre_ParVectorLocalVector( reinterpret_cast< hypre_ParVector * >( u ) ) );

  HYPRE_Int const nf0 = LvArray::integerConversion< HYPRE_Int >( data.dofOfFlux0.size() );
  HYPRE_Int const nf1 = LvArray::integerConversion< HYPRE_Int >( data.dofOfFlux1.size() );
  HYPRE_Int const np = LvArray::integerConversion< HYPRE_Int >( data.dofOfPres.size() );

  // pressures: the conservation residual taken to the symmetrized system, one AMG cycle on the block
  HYPRE_Real * const xp = dataOf( data.xp );
  HYPRE_Real * const rt = dataOf( data.rt );
  HYPRE_Real * const xt = dataOf( data.xt );
  for( HYPRE_Int e = 0; e < np; ++e )
  {
    rt[e] = data.rowScale[e] * r[data.dofOfPres[e]];
    xt[e] = 0.0;
  }
  if( np > 0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSolve( data.amg,
                                                reinterpret_cast< HYPRE_ParCSRMatrix >( parCSROf( data.pIJ ) ),
                                                reinterpret_cast< HYPRE_ParVector >( parVectorOf( data.rt ) ),
                                                reinterpret_cast< HYPRE_ParVector >( parVectorOf( data.xt ) ) ) );
  }
  for( HYPRE_Int e = 0; e < np; ++e )
  {
    xp[e] = xt[e];
    x[data.dofOfPres[e]] = xp[e];
  }

  // live fluxes: the flux block alone (block-diagonal map), pinned dofs by their diagonal
  if( nf1 > 0 )
  {
    HYPRE_Real * const r1 = dataOf( data.r1 );
    HYPRE_Real * const x1 = dataOf( data.x1 );
    for( HYPRE_Int i = 0; i < nf1; ++i )
    {
      r1[i] = r[data.dofOfFlux1[i]];
      x1[i] = 0.0;
    }
    if( data.adsActive )
    {
      // a CG that stops on its iteration cap reports a convergence error: not fatal here
      HYPRE_ParCSRPCGSolve( data.pcg,
                            reinterpret_cast< HYPRE_ParCSRMatrix >( parCSROf( data.nIJ ) ),
                            reinterpret_cast< HYPRE_ParVector >( parVectorOf( data.r1 ) ),
                            reinterpret_cast< HYPRE_ParVector >( parVectorOf( data.x1 ) ) );
      HYPRE_ClearAllErrors();
    }
    for( HYPRE_Int i = 0; i < nf1; ++i )
    {
      x[data.dofOfFlux1[i]] = data.pinned[i] ? data.invDiagPinned[i] * r1[i] : x1[i];
    }
  }

  // condensed fluxes: exact closure q0 = D0^{-1} ( r0 - F0 p )
  if( nf0 > 0 )
  {
    GEOS_LAI_CHECK_ERROR( hypre_ParCSRMatrixMatvec( 1.0, parCSROf( data.b0tIJ ), parVectorOf( data.xp ),
                                                    0.0, parVectorOf( data.t0 ) ) );
    HYPRE_Real const * const t0 = dataOf( data.t0 );
    for( HYPRE_Int k = 0; k < nf0; ++k )
    {
      HYPRE_Int const dof = data.dofOfFlux0[k];
      x[dof] = data.invD0[k] * ( r[dof] - t0[k] );
    }
  }
  return 0;
}

HYPRE_Int rieszDestroy( HYPRE_Solver solver )
{
  RieszMFDData * const data = reinterpret_cast< RieszMFDData * >( solver );
  destroySetupObjects( *data );
  destroyIJ( data->curlIJ );
  destroyIJ( data->gradIJ );
  for( int d = 0; d < 3; ++d )
  {
    destroyIJ( data->coordIJ[d] );
  }
  delete data;
  return 0;
}

} // namespace

void createRieszMFD( LinearSolverParameters const & params,
                     HyprePrecWrapper & precond )
{
  GEOS_ERROR_IF( MpiWrapper::commSize( MPI_COMM_WORLD ) > 1,
                 "Riesz MFD preconditioner: only serial runs are supported at the moment" );
  GEOS_ERROR_IF( params.mgr.customPointMarkers.empty(),
                 "Riesz MFD preconditioner: the solver must provide the dof markers (0 = condensed face, 1 = live face, 2 = pressure)" );

  auto * const data = new RieszMFDData{};
  localIndex const n = params.mgr.customPointMarkers.size();
  data->kindOf.resize( n );
  data->blockRow.resize( n );
  for( localIndex i = 0; i < n; ++i )
  {
    HYPRE_Int const kind = LvArray::integerConversion< HYPRE_Int >( params.mgr.customPointMarkers[i] );
    data->kindOf[i] = kind;
    stdVector< HYPRE_Int > & block = kind == condensedFaceMarker ? data->dofOfFlux0 :
                                     kind == liveFaceMarker ? data->dofOfFlux1 : data->dofOfPres;
    GEOS_ERROR_IF( kind != condensedFaceMarker && kind != liveFaceMarker && kind != pressureMarker,
                   GEOS_FMT( "Riesz MFD preconditioner: unknown dof marker {}", kind ) );
    data->blockRow[i] = LvArray::integerConversion< HYPRE_Int >( block.size() );
    block.push_back( LvArray::integerConversion< HYPRE_Int >( i ) );
  }

  LinearSolverParameters::ADSAuxData const & aux = params.adsAuxData;
  GEOS_ERROR_IF( aux.mfdCell.size() != n || aux.pressureNormScale.size() != n,
                 "Riesz MFD preconditioner: the solver must provide the MFD flag and the L2 scale of every pressure dof" );
  data->isMfd.resize( data->dofOfPres.size() );
  data->normScale.resize( data->dofOfPres.size() );
  for( std::size_t e = 0; e < data->dofOfPres.size(); ++e )
  {
    data->isMfd[e] = aux.mfdCell[data->dofOfPres[e]];
    data->normScale[e] = aux.pressureNormScale[data->dofOfPres[e]];
  }
  globalIndex const numFluxRows = aux.cRowPtr.size() - 1;
  GEOS_ERROR_IF_NE_MSG( numFluxRows, LvArray::integerConversion< globalIndex >( data->dofOfFlux1.size() ),
                        "Riesz MFD preconditioner: the discrete curl's rows are not the live flux dofs" );
  if( numFluxRows > 0 )
  {
    globalIndex const numEdges = aux.gRowPtr.size() - 1;
    globalIndex const numVertices = aux.xCoords.size();
    data->curlIJ = createIJFromCSR( aux.cRowPtr, aux.cCols, aux.cVals, numEdges );
    data->gradIJ = createIJFromCSR( aux.gRowPtr, aux.gCols, aux.gVals, numVertices );
    data->coordIJ[0] = createShiftedCoordinates( aux.xCoords );
    data->coordIJ[1] = createShiftedCoordinates( aux.yCoords );
    data->coordIJ[2] = createShiftedCoordinates( aux.zCoords );
  }

  data->base.setup = reinterpret_cast< HYPRE_PtrToSolverFcn >( rieszSetup );
  data->base.solve = reinterpret_cast< HYPRE_PtrToSolverFcn >( rieszSolve );
  data->base.destroy = reinterpret_cast< HYPRE_PtrToDestroyFcn >( rieszDestroy );
  data->base.is_setup = 0;

  precond.ptr = reinterpret_cast< HYPRE_Solver >( data );
  precond.setup = rieszSetup;
  precond.solve = rieszSolve;
  precond.destroy = rieszDestroy;
}

} // namespace hypre

} // namespace geos
