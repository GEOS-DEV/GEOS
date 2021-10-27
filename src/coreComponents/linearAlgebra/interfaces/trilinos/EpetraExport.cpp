/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file EpetraExport.cpp
 */

#include "EpetraExport.hpp"

#include "common/MpiWrapper.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraMatrix.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraUtils.hpp"

#include <Epetra_Map.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_FEVector.h>
#include <Epetra_Import.h>

namespace geosx
{

EpetraExport::EpetraExport() = default;

EpetraExport::EpetraExport( EpetraMatrix const & mat,
                            integer targetRank )
  : m_targetRank( targetRank )
{
  globalIndex const numGlobalRows = mat.numGlobalRows();
  localIndex const numLocalRows = ( m_targetRank == MpiWrapper::commRank( mat.getComm() ) ) ? numGlobalRows : 0;
  m_serialMap = std::make_unique< Epetra_Map >( numGlobalRows, numLocalRows, 0, mat.unwrapped().Comm() );
  m_serialImport = std::make_unique< Epetra_Import >( *m_serialMap, mat.unwrapped().RowMap() );
}

EpetraExport::~EpetraExport() = default;

template< typename OFFSET_TYPE, typename COLUMN_TYPE >
void EpetraExport::exportCRS( EpetraMatrix const & mat,
                              arrayView1d< OFFSET_TYPE > const & rowOffsets,
                              arrayView1d< COLUMN_TYPE > const & colIndices,
                              arrayView1d< real64 > const & values ) const
{

  int const rank = MpiWrapper::commRank( mat.getComm() );
  Epetra_CrsMatrix const * localMatrix = &mat.unwrapped();

  // import on target rank if needed
  std::unique_ptr< Epetra_CrsMatrix > serialMatrix;
  if( m_targetRank >= 0 )
  {
//    GEOSX_ERROR_IF( rank == m_targetRank && !( rowOffsets && colIndices && values ),
//                    "EpetraExport: must pass non-null pointers on the target rank" );
    serialMatrix = std::make_unique< Epetra_CrsMatrix >( mat.unwrapped(), *m_serialImport, m_serialMap.get() );
    localMatrix = serialMatrix.get();
  }

  // export the raw CRS data
  if( m_targetRank < 0 || m_targetRank == rank )
  {
    int * ia;
    int * ja;
    real64 * va;
    GEOSX_LAI_CHECK_ERROR( localMatrix->ExtractCrsDataPointers( ia, ja, va ) );

    // contains the global ID of local columns
    globalIndex const * const globalColumns = localMatrix->ColMap().MyGlobalElements64();

    ///////////////
    rowOffsets.move( LvArray::MemorySpace::host, false );
    colIndices.move( LvArray::MemorySpace::host, false );
    values.move( LvArray::MemorySpace::host, false );
    /////////

    std::transform( ia, ia + localMatrix->NumMyRows() + 1, rowOffsets.data(), LvArray::integerConversion< OFFSET_TYPE, int > );
    std::transform( ja, ja + localMatrix->NumMyNonzeros(), colIndices.data(),
                    [globalColumns]( int const i ){ return LvArray::integerConversion< COLUMN_TYPE >( globalColumns[i] ); } );
    std::copy( va, va + localMatrix->NumMyNonzeros(), values.data() );
  }
}

void EpetraExport::exportVector( EpetraVector const & vec,
                                 real64 * values ) const
{
  if( m_targetRank >= 0 )
  {
    int const rank = MpiWrapper::commRank( vec.getComm() );
    GEOSX_ERROR_IF( rank == m_targetRank && !values, "EpetraExport: must pass non-null pointers on the target rank" );

    // Create a local vector that directly wraps the user-provided buffer and gather
    Epetra_MultiVector localVec( View, *m_serialMap, values, 0, 1 );
    localVec.Import( vec.unwrapped(), *m_serialImport, Insert );
  }
  else
  {
    real64 const * const data = vec.extractLocalVector();
    std::copy( data, data + vec.localSize(), values );
  }
}

void EpetraExport::importVector( real64 const * values,
                                 EpetraVector & vec ) const
{
  if( m_targetRank >= 0 )
  {
    int const rank = MpiWrapper::commRank( vec.getComm() );
    GEOSX_ERROR_IF( rank == m_targetRank && !values, "EpetraExport: must pass non-null pointers on the target rank" );

    // HACK: const_cast required in order to create an Epetra vector that wraps user data;
    //       we promise the values won't be changed since the only thing we do is scatter.
    Epetra_MultiVector localVector( View, *m_serialMap, const_cast< real64 * >( values ), LvArray::integerConversion< int >( vec.globalSize() ), 1 );
    vec.unwrapped().Export( localVector, *m_serialImport, Insert );
  }
  else
  {
    real64 * const data = vec.extractLocalVector();
    std::copy( values, values + vec.localSize(), data );
  }
}

/**
 * Explicit template instantiation macro for EpetraExport::exportCRS.
 * We need to explicitly instantiate this function template because:
 * - we want CRS consumers to specify their own destination buffer types;
 * - we're "hiding" Epetra headers from being included by consumer code.
 */
#define INST_EPETRA_EXPORT_CRS( OFFSET_TYPE, COLUMN_TYPE ) \
  template void \
  EpetraExport::exportCRS< OFFSET_TYPE, COLUMN_TYPE >( EpetraMatrix const &, \
                                                       arrayView1d< OFFSET_TYPE > const &, \
                                                       arrayView1d< COLUMN_TYPE > const &, \
                                                       arrayView1d< real64 > const & ) const

// Add other instantiations as needed (only use built-in types)
INST_EPETRA_EXPORT_CRS( int, int );
INST_EPETRA_EXPORT_CRS( long, long );
INST_EPETRA_EXPORT_CRS( long long, long long );

#undef INST_EPETRA_EXPORT_CRS

} // namespace geosx
