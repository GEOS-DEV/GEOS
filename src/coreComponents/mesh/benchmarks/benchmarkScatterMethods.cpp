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
 * @file benchmarkScatterMethods.cpp
 * @brief Benchmark comparing VTK mesh scattering methods with load balance analysis
 */

#include "common/DataTypes.hpp"
#include "common/format/Format.hpp"
#include "common/logger/Logger.hpp"
#include "common/MpiWrapper.hpp"
#include "mesh/generators/VTKMeshScattering.hpp"

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkMultiProcessController.h>
#include <vtkMPIController.h>
#include <vtkMPICommunicator.h>
#include <vtkMPI.h>
#include <vtkRedistributeDataSetFilter.h>
#include <vtkIntArray.h>
#include <vtkAbstractArray.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>

using namespace geos;
using namespace geos::vtk;

namespace
{

// ============================================================================
// Load balance and compactness metrics
// ============================================================================

struct LoadBalanceStats
{
  long long totalCells;
  long long minCells;
  long long maxCells;
  double avgCells;
  double imbalanceRatio;   // max / avg (1.0 = perfect)
  double bboxOverlap;      // sum-of-local-bbox-volumes / total-bbox-volume (1.0 = ideal)
};

/**
 * @brief Compute bounding-box overlap ratio as a proxy for spatial compactness.
 *
 * Each rank's local bounding box volume is summed across all ranks and divided
 * by the total mesh bounding box volume.  A ratio near 1.0 means partitions
 * tile the domain with little overlap; a ratio near @p nranks means every
 * partition spans the entire mesh (no spatial locality).
 */
double computeBboxOverlap( vtkDataSet * localMesh,
                           double const totalBounds[6],
                           MPI_Comm comm )
{
  bool const hasLocal = localMesh->GetNumberOfCells() > 0;

  double localBounds[6] = { 0, 0, 0, 0, 0, 0 };
  if( hasLocal )
  {
    localMesh->GetBounds( localBounds );
  }

  // Only use non-degenerate dimensions (extent > epsilon).
  double totalVol = 1.0;
  double localVol = hasLocal ? 1.0 : 0.0;
  bool hasDim = false;
  for( int d = 0; d < 3; ++d )
  {
    double const totalExtent = totalBounds[2 * d + 1] - totalBounds[2 * d];
    if( totalExtent > 1.0e-12 )
    {
      hasDim = true;
      totalVol *= totalExtent;
      if( hasLocal )
      {
        localVol *= std::max( localBounds[2 * d + 1] - localBounds[2 * d], 0.0 );
      }
    }
  }

  if( !hasDim || totalVol <= 0.0 )
  {
    return 0.0;
  }

  double const sumLocalVol = MpiWrapper::allReduce( localVol, MpiWrapper::Reduction::Sum, comm );
  return sumLocalVol / totalVol;
}

LoadBalanceStats computeLoadBalance( long long localCells,
                                     vtkDataSet * localMesh,
                                     double const totalBounds[6],
                                     MPI_Comm comm )
{
  LoadBalanceStats stats{};
  int const size = MpiWrapper::commSize( comm );

  stats.totalCells = MpiWrapper::allReduce( localCells, MpiWrapper::Reduction::Sum, comm );
  stats.minCells = MpiWrapper::allReduce( localCells, MpiWrapper::Reduction::Min, comm );
  stats.maxCells = MpiWrapper::allReduce( localCells, MpiWrapper::Reduction::Max, comm );
  stats.avgCells = static_cast< double >( stats.totalCells ) / size;
  stats.imbalanceRatio = ( stats.avgCells > 0 )
    ? stats.maxCells / stats.avgCells : 0.0;
  stats.bboxOverlap = computeBboxOverlap( localMesh, totalBounds, comm );

  return stats;
}

// ============================================================================
// Benchmark result and output
// ============================================================================

struct BenchmarkResult
{
  std::string name;
  double time;
  LoadBalanceStats balance;
  bool cellsPreserved;
};

/**
 * @brief Save a partitioned mesh to a per-method subdirectory.
 *
 * Each rank writes its own .vtu piece.  Only rank 0 writes the .pvtu header.
 * Directory layout:  <dir>/<methodName>/<methodName>_<rank>.vtu
 *                    <dir>/<methodName>/<methodName>.pvtu
 */
void savePartitionedMesh( vtkUnstructuredGrid * mesh,
                          std::string const & dir,
                          std::string const & methodName,
                          MPI_Comm comm )
{
  int const rank = MpiWrapper::commRank( comm );

  // Add RankID cell data
  vtkNew< vtkIntArray > rankArray;
  rankArray->SetName( "RankID" );
  rankArray->SetNumberOfTuples( mesh->GetNumberOfCells() );
  rankArray->FillValue( rank );
  mesh->GetCellData()->AddArray( rankArray );

  std::string const methodDir = GEOS_FMT( "{}/{}", dir, methodName );
  if( rank == 0 )
  {
    std::filesystem::create_directories( methodDir );
  }
  MpiWrapper::barrier( comm );

  // Write using the parallel XML writer
  std::string const pvtuFile = GEOS_FMT( "{}/{}.pvtu", methodDir, methodName );

  vtkNew< vtkMPIController > controller;
  vtkMPICommunicatorOpaqueComm vtkComm( &comm );
  vtkNew< vtkMPICommunicator > communicator;
  communicator->InitializeExternal( &vtkComm );
  controller->SetCommunicator( communicator );

  vtkNew< vtkXMLPUnstructuredGridWriter > writer;
  writer->SetController( controller );
  writer->SetFileName( pvtuFile.c_str() );
  writer->SetInputData( mesh );
  writer->SetNumberOfPieces( MpiWrapper::commSize( comm ) );
  writer->SetStartPiece( rank );
  writer->SetEndPiece( rank );
  writer->Write();

  controller->SetCommunicator( nullptr );

  MpiWrapper::barrier( comm );

  if( rank == 0 )
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "  Saved: {}", pvtuFile ) );
  }
}

BenchmarkResult runBenchmark( std::string const & name,
                              ScatterMethod method,
                              vtkUnstructuredGrid & mesh,
                              arrayView1d< int const > const & partitions,
                              long long totalCells,
                              double const totalBounds[6],
                              std::string const & outputDir,
                              MPI_Comm comm )
{
  BenchmarkResult result;
  result.name = name;

  MpiWrapper::barrier( comm );

  auto start = std::chrono::high_resolution_clock::now();

  vtkSmartPointer< vtkDataSet > scattered;
  if( method == ScatterMethod::kdtree )
  {
    // Call VTK redistribution directly to get raw kdtree metrics,
    // bypassing the contiguous fallback in scatterMesh.
    int const size = MpiWrapper::commSize( comm );
    vtkNew< vtkMPIController > controller;
    vtkMPICommunicatorOpaqueComm vtkComm( &comm );
    vtkNew< vtkMPICommunicator > communicator;
    communicator->InitializeExternal( &vtkComm );
    controller->SetCommunicator( communicator );
    vtkMultiProcessController::SetGlobalController( controller );

    vtkNew< vtkRedistributeDataSetFilter > rdsf;
    rdsf->SetInputDataObject( &mesh );
    rdsf->SetNumberOfPartitions( size );
    rdsf->Update();
    scattered = vtkDataSet::SafeDownCast( rdsf->GetOutputDataObject( 0 ) );
  }
  else
  {
    scattered = scatterMesh( method, mesh, partitions, comm );
  }
  vtkSmartPointer< vtkUnstructuredGrid > output =
    vtkUnstructuredGrid::SafeDownCast( scattered );

  auto end = std::chrono::high_resolution_clock::now();
  double localTime = std::chrono::duration< double >( end - start ).count();

  long long localCells = output->GetNumberOfCells();

  result.time = MpiWrapper::reduce( localTime, MpiWrapper::Reduction::Max, 0, comm );
  result.balance = computeLoadBalance( localCells, output, totalBounds, comm );
  result.cellsPreserved = ( result.balance.totalCells == totalCells );

  GEOS_LOG_RANK_0( GEOS_FMT( "\n=== {} ===\n"
                             "  Time: {:.3f}s\n"
                             "  Cells: {}{}\n"
                             "  Balance: min={} max={} avg={:.1f}"
                             " imbalance={:.6f}x bbox_overlap={:.3f}",
                             name,
                             result.time,
                             result.balance.totalCells,
                             result.cellsPreserved ? ""
                                : GEOS_FMT( "  WARNING: lost {} cells!",
                                            totalCells - result.balance.totalCells ),
                             result.balance.minCells,
                             result.balance.maxCells,
                             result.balance.avgCells,
                             result.balance.imbalanceRatio,
                             result.balance.bboxOverlap ) );

  // Save partitioned mesh immediately if an output directory was provided.
  if( !outputDir.empty() )
  {
    savePartitionedMesh( output, outputDir, name, comm );
  }

  MpiWrapper::barrier( comm );
  return result;
}

void printSummary( stdVector< BenchmarkResult > const & results,
                   long long totalCells,
                   MPI_Comm comm )
{
  int const size = MpiWrapper::commSize( comm );

  std::string const sep( 110, '=' );
  std::string const thin( 110, '-' );

  GEOS_LOG_RANK_0( GEOS_FMT( "\n{}\nSUMMARY ({} ranks, {} cells)\n{}",
                             sep, size, totalCells, sep ) );

  // Table header
  GEOS_LOG_RANK_0( GEOS_FMT( "{:<14}| {:<10}| {:<12}| {:<8}| {:<10}| {:<10}| {:<12}| {:<12}",
                             "Method", "Time (s)", "Cells", "Status",
                             "Min cells", "Max cells", "Imbalance", "BBox Overlap" ) );
  GEOS_LOG_RANK_0( thin );

  // Table rows
  for( auto const & r : results )
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "{:<14}| {:<10.3f}| {:<12}| {:<8}| {:<10}| {:<10}| {:<12.6f}| {:<12.3f}",
                               r.name, r.time,
                               r.balance.totalCells,
                               r.cellsPreserved ? "OK" : "LOSS",
                               r.balance.minCells,
                               r.balance.maxCells,
                               r.balance.imbalanceRatio,
                               r.balance.bboxOverlap ) );
  }
  GEOS_LOG_RANK_0( thin );

  // Find fastest method
  double fastestTime = results[0].time;
  for( auto const & r : results )
  {
    fastestTime = std::min( fastestTime, r.time );
  }

  // Speedup table
  GEOS_LOG_RANK_0( "\nSpeedup relative to each method:" );
  for( auto const & r : results )
  {
    double const speedup = r.time / fastestTime;
    GEOS_LOG_RANK_0( GEOS_FMT( "  {:<14}{:.2f}x slower than fastest{}",
                               r.name, speedup,
                               ( std::abs( r.time - fastestTime ) < 1e-6 ) ? " <-- fastest" : "" ) );
  }

  // Best balance and compactness
  double bestImbalance = results[0].balance.imbalanceRatio;
  double bestOverlap = results[0].balance.bboxOverlap;
  std::string bestBalanceName = results[0].name;
  std::string bestCompactName = results[0].name;
  for( auto const & r : results )
  {
    if( r.cellsPreserved )
    {
      if( r.balance.imbalanceRatio < bestImbalance )
      {
        bestImbalance = r.balance.imbalanceRatio; bestBalanceName = r.name;
      }
      if( r.balance.bboxOverlap < bestOverlap )
      {
        bestOverlap = r.balance.bboxOverlap; bestCompactName = r.name;
      }
    }
  }
  GEOS_LOG_RANK_0( GEOS_FMT( "\nBest load balance: {} (imbalance ratio: {:.6f}x)",
                             bestBalanceName, bestImbalance ) );
  GEOS_LOG_RANK_0( GEOS_FMT( "Best compactness:  {} (bbox overlap: {:.3f})", bestCompactName, bestOverlap ) );
  GEOS_LOG_RANK_0( sep );
}

} // anonymous namespace

int main( int argc, char * argv[] )
{
  MpiWrapper::init( &argc, &argv );

#ifdef GEOS_USE_CHAI
  chai::ArrayManager::getInstance()->disableCallbacks();
#endif

  MPI_Comm const comm = MPI_COMM_WORLD;
  int const size = MpiWrapper::commSize( comm );

  logger::InitializeLogger( comm );

  if( argc < 5 )
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "Usage: {} <input.vtu> <nx> <ny> <nz> [output_dir]\n"
                               "  nx*ny*nz must equal number of MPI ranks (for Cartesian method)\n"
                               "  output_dir (optional): save partitioned meshes (one subfolder per method)",
                               argv[0] ) );
    logger::FinalizeLogger();
    MpiWrapper::finalize();
    return 1;
  }

  int const nx = std::atoi( argv[2] );
  int const ny = std::atoi( argv[3] );
  int const nz = std::atoi( argv[4] );

  std::string outputDir;
  if( argc > 5 )
  {
    outputDir = argv[5];
  }

  GEOS_LOG_RANK_0( "=== Redistribution Performance Comparison ===" );
  GEOS_LOG_RANK_0( GEOS_FMT( "MPI Ranks: {}", size ) );
  GEOS_LOG_RANK_0( GEOS_FMT( "Grid: {} x {} x {}", nx, ny, nz ) );
  GEOS_LOG_RANK_0( GEOS_FMT( "Input file: {}", argv[1] ) );

  // Load mesh on rank 0
  vtkSmartPointer< vtkUnstructuredGrid > mesh;
  long long totalCells = 0;
  double totalBounds[6] = { 0, 0, 0, 0, 0, 0 };

  int const rank = MpiWrapper::commRank( comm );
  if( rank == 0 )
  {
    vtkNew< vtkXMLUnstructuredGridReader > reader;
    reader->SetFileName( argv[1] );
    reader->Update();
    mesh = reader->GetOutput();

    if( !mesh || mesh->GetNumberOfCells() == 0 )
    {
      std::cerr << "Failed to read mesh or mesh is empty\n";
      MPI_Abort( comm, 1 );
    }

    totalCells = mesh->GetNumberOfCells();
    mesh->GetBounds( totalBounds );
  }
  else
  {
    mesh = vtkSmartPointer< vtkUnstructuredGrid >::New();
  }

  MpiWrapper::bcast( &totalCells, 1, 0, comm );
  MpiWrapper::bcast( totalBounds, 6, 0, comm );

  GEOS_LOG_RANK_0( GEOS_FMT( "Loaded {} cells", totalCells ) );
  GEOS_LOG_RANK_0( GEOS_FMT( "Bounds: [{}, {}] x [{}, {}] x [{}, {}]",
                             totalBounds[0], totalBounds[1],
                             totalBounds[2], totalBounds[3],
                             totalBounds[4], totalBounds[5] ) );

  stdVector< BenchmarkResult > results;

  array1d< int > parts( 3 );
  parts[0] = nx; parts[1] = ny; parts[2] = nz;

  auto run = [&]( std::string const & name, ScatterMethod method )
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "\nRunning {} scatter ({} cells, {} ranks)...", name, totalCells, size ) );
    results.push_back( runBenchmark( name, method, *mesh, parts.toViewConst(),
                                     totalCells, totalBounds, outputDir, comm ) );
  };

  run( "KdTree", ScatterMethod::kdtree );
  run( "Contiguous", ScatterMethod::contiguous );

  if( nx * ny * nz == size )
  {
    run( "Cartesian", ScatterMethod::cartesian );
  }
  else
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "\n=== Cartesian ===\n  SKIPPED: nx*ny*nz ({}) != MPI size ({})",
                               nx * ny * nz, size ) );
  }

  run( "RCB", ScatterMethod::rcb );

  // Print summary
  printSummary( results, totalCells, comm );

  logger::FinalizeLogger();
  MpiWrapper::finalize();
  return 0;
}
