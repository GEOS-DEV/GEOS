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

/** @file benchmarkUniformRefinement.cpp
 * Component benchmark: local volume templates, incidence, point IDs and fields.
 * It does not measure the coupled AllMeshes importer or GEOS ghosts/DOFs.
 */

#include "common/TimingMacros.hpp"
#include "mesh/generators/VTKMeshScattering.hpp"
#include "mesh/generators/VTKRefinementCommunication.hpp"
#include "mesh/generators/VTKRefinementFields.hpp"
#include "mesh/generators/VTKRefinementTemplates.hpp"
#include "mesh/generators/VTKRefinementSharing.hpp"

#include <vtkCellData.h>
#include <vtkArrayDispatch.h>
#include <vtkCellType.h>
#include <vtkDataArray.h>
#include <vtkDataArrayAccessor.h>
#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <algorithm>
#include <chrono>
#include <climits>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <sys/resource.h>
#ifdef __linux__
#include <unistd.h>
#endif

using namespace geos;
using namespace geos::vtk::refinement;

namespace
{
using Clock = std::chrono::steady_clock;

struct Options
{
  std::array< int, 3 > cells{ 8, 8, 8 }, grid{ 0, 0, 0 };
  int levels = 1, pointComponents = 0;
  bool weak = false;
  std::string kind = "hex", input, label = "unspecified";
  std::string discovery = "boundary", sharing = "interfaces";
  std::uint64_t chunkBytes = UINT64_C( 1 ) << 30;
};

int parseInteger( char const * value, bool zero = false )
{
  std::size_t consumed = 0;
  long long const parsed = std::stoll( value, &consumed );
  if( value[consumed] != '\0' || parsed < ( zero ? 0 : 1 ) || parsed > INT_MAX )
    throw std::invalid_argument( "Expected a representable nonnegative/positive integer" );
  return static_cast< int >( parsed );
}

Options options( int argc, char ** argv, int size )
{
  Options result;
  for( int i = 1; i < argc; ++i )
  {
    std::string const option = argv[i];
    auto next = [&]()
    {
      if( ++i == argc )
        throw std::invalid_argument( "Missing value for " + option );
      return argv[i];
    };
    if( option == "--cells" || option == "--process-grid" )
      for( int & value : option == "--cells" ? result.cells : result.grid )
        value = parseInteger( next() );
    else if( option == "--levels" )
      result.levels = parseInteger( next(), true );
    else if( option == "--point-components" )
      result.pointComponents = parseInteger( next(), true );
    else if( option == "--kind" )
      result.kind = next();
    else if( option == "--input" )
      result.input = next();
    else if( option == "--discovery" )
      result.discovery = next();
    else if( option == "--sharing" )
      result.sharing = next();
    else if( option == "--label" )
      result.label = next();
    else if( option == "--chunk-bytes" )
      result.chunkBytes = parseInteger( next() );
    else if( option == "--weak" )
      result.weak = true;
    else
      throw std::invalid_argument( "Unknown option " + option );
  }
  if( result.kind != "hex" && result.kind != "tet" && result.kind != "pyramid" )
    throw std::invalid_argument( "--kind must be hex, tet or pyramid" );
  if( result.discovery != "boundary" && result.discovery != "all" )
    throw std::invalid_argument( "--discovery must be boundary or all" );
  if( result.sharing != "interfaces" && result.sharing != "all" )
    throw std::invalid_argument( "--sharing must be interfaces or all" );
  if( !result.input.empty() && result.weak )
    throw std::invalid_argument( "--weak is only defined for generated meshes" );
#ifdef GEOS_USE_MPI
  if( result.grid[0] == 0 )
    MPI_Dims_create( size, 3, result.grid.data() );
#else
  if( result.grid[0] == 0 )
    result.grid = { 1, 1, 1 };
#endif
  std::uint64_t ranks = 1;
  for( int d = 0; d < 3; ++d )
  {
    if( ranks > static_cast< std::uint64_t >( size ) / result.grid[d] )
      throw std::invalid_argument( "Process-grid product exceeds MPI size" );
    ranks *= result.grid[d];
    if( result.weak )
    {
      if( result.cells[d] > INT_MAX / result.grid[d] )
        throw std::overflow_error( "Weak-scaling grid dimension overflow" );
      result.cells[d] *= result.grid[d];
    }
  }
  if( ranks != static_cast< std::uint64_t >( size ) )
    throw std::invalid_argument( "Process-grid product must equal MPI size" );
  if( !result.input.empty() )
  {
    if( result.pointComponents )
      throw std::invalid_argument( "--point-components applies only to generated meshes" );
    result.kind = "input";
  }
  return result;
}

std::uint64_t product( std::array< int, 3 > const & dimensions, int increment = 0 )
{
  std::uint64_t result = 1;
  for( int dimension : dimensions )
  {
    auto const factor = static_cast< std::uint64_t >( dimension ) + increment;
    if( result > static_cast< std::uint64_t >( std::numeric_limits< vtkIdType >::max() ) / factor )
      throw std::overflow_error( "Benchmark grid exceeds vtkIdType" );
    result *= factor;
  }
  return result;
}

struct Mesh
{
  std::vector< Coordinates > coordinates;
  Connectivity pointIds, cellIds;
  std::vector< Cell > cells;
  std::vector< EntitySupport > supports;
  bool sharingKnown = false;
  vtkSmartPointer< vtkPointData > pointData = vtkSmartPointer< vtkPointData >::New();
  vtkSmartPointer< vtkCellData > cellData = vtkSmartPointer< vtkCellData >::New();
};

Mesh generatedMesh( Options const & opt, int rank )
{
  Mesh mesh;
  auto const nodeCount = product( opt.cells, 1 ), cubeCount = product( opt.cells );
  int const perCube = opt.kind == "hex" ? 1 : 6;
  if( cubeCount > static_cast< std::uint64_t >( std::numeric_limits< vtkIdType >::max() ) / perCube ||
      ( opt.kind == "pyramid" && cubeCount > static_cast< std::uint64_t >( std::numeric_limits< vtkIdType >::max() ) - nodeCount ) )
    throw std::overflow_error( "Benchmark coarse IDs overflow" );
  std::array< int, 3 > begin{}, end{}, local{}, dims{};
  int remainder = rank;
  for( int d = 0; d < 3; ++d )
  {
    int const coordinate = remainder % opt.grid[d];
    remainder /= opt.grid[d];
    begin[d] = static_cast< int >( static_cast< std::int64_t >( opt.cells[d] ) * coordinate / opt.grid[d] );
    end[d] = static_cast< int >( static_cast< std::int64_t >( opt.cells[d] ) * ( coordinate + 1 ) / opt.grid[d] );
    dims[d] = end[d] - begin[d];
    local[d] = dims[d] + 1;
  }
  if( std::find( dims.begin(), dims.end(), 0 ) != dims.end() )
    return mesh;
  auto node = [&]( int x, int y, int z ) -> vtkIdType
  { return x + static_cast< vtkIdType >( local[0] ) * ( y + static_cast< vtkIdType >( local[1] ) * z ); };
  for( int z = begin[2]; z <= end[2]; ++z )
    for( int y = begin[1]; y <= end[1]; ++y )
      for( int x = begin[0]; x <= end[0]; ++x )
      {
        mesh.coordinates.push_back( { double( x ), double( y ), double( z ) } );
        mesh.pointIds.push_back( x + ( static_cast< vtkIdType >( opt.cells[0] ) + 1 ) *
                                         ( y + ( static_cast< vtkIdType >( opt.cells[1] ) + 1 ) * z ) );
      }
  for( int z = 0; z < dims[2]; ++z )
    for( int y = 0; y < dims[1]; ++y )
      for( int x = 0; x < dims[0]; ++x )
      {
        Cell const cube{ VTK_HEXAHEDRON,
                         { node( x, y, z ), node( x + 1, y, z ), node( x + 1, y + 1, z ), node( x, y + 1, z ), node( x, y, z + 1 ),
                           node( x + 1, y, z + 1 ), node( x + 1, y + 1, z + 1 ), node( x, y + 1, z + 1 ) },
                         0 };
        vtkIdType const cubeId =
            x + begin[0] +
            static_cast< vtkIdType >( opt.cells[0] ) * ( y + begin[1] + static_cast< vtkIdType >( opt.cells[1] ) * ( z + begin[2] ) );
        if( opt.kind == "hex" )
        {
          mesh.cells.push_back( cube );
          mesh.cellIds.push_back( cubeId );
        }
        else if( opt.kind == "tet" )
        {
          // Freudenthal subdivision: consistent face diagonals on the whole grid.
          constexpr int tetrahedra[6][4] = { { 0, 1, 2, 6 }, { 0, 2, 3, 6 }, { 0, 3, 7, 6 },
                                             { 0, 7, 4, 6 }, { 0, 4, 5, 6 }, { 0, 5, 1, 6 } };
          for( int i = 0; i < 6; ++i )
          {
            Connectivity corners;
            for( int vertex : tetrahedra[i] )
              corners.push_back( cube.points[vertex] );
            mesh.cells.push_back( { VTK_TETRA, std::move( corners ), 0 } );
            mesh.cellIds.push_back( 6 * cubeId + i );
          }
        }
        else
        {
          vtkIdType const center = mesh.coordinates.size();
          mesh.coordinates.push_back( { x + begin[0] + .5, y + begin[1] + .5, z + begin[2] + .5 } );
          mesh.pointIds.push_back( nodeCount + cubeId );
          auto faces = cellFaces( cube );
          for( std::size_t i = 0; i < faces.size(); ++i )
          {
            std::reverse( faces[i].begin(), faces[i].end() );
            faces[i].push_back( center );
            mesh.cells.push_back( { VTK_PYRAMID, std::move( faces[i] ), 0 } );
            mesh.cellIds.push_back( 6 * cubeId + i );
          }
        }
      }
  if( opt.pointComponents )
  {
    vtkNew< vtkDoubleArray > field;
    field->SetName( "benchmarkPointField" );
    field->SetNumberOfComponents( opt.pointComponents );
    field->SetNumberOfTuples( mesh.coordinates.size() );
    for( std::size_t i = 0; i < mesh.coordinates.size(); ++i )
      for( int component = 0; component < opt.pointComponents; ++component )
        field->SetTypedComponent( i, component,
                                  component + mesh.coordinates[i][0] + 2 * mesh.coordinates[i][1] + 3 * mesh.coordinates[i][2] );
    mesh.pointData->AddArray( field );
  }
  return mesh;
}

std::vector< EntitySupport > incidence( Mesh const & mesh, int rank )
{
  return meshEntitySupports( mesh.cells, mesh.cellIds, mesh.pointIds, rank );
}

std::uint64_t memoryKiB( bool peak )
{
  struct rusage usage{};
  if( peak )
  {
    if( getrusage( RUSAGE_SELF, &usage ) != 0 )
      throw std::runtime_error( "Cannot read peak resident memory" );
#ifdef __APPLE__
    return usage.ru_maxrss / 1024;
#else
    return usage.ru_maxrss;
#endif
  }
#ifdef __linux__
  std::ifstream input( "/proc/self/statm" );
  std::uint64_t virtualPages = 0, residentPages = 0;
  if( !( input >> virtualPages >> residentPages ) )
    throw std::runtime_error( "Cannot read current resident memory" );
  return residentPages * static_cast< std::uint64_t >( sysconf( _SC_PAGESIZE ) ) / 1024;
#else
  return 0; // Current RSS is only reported on Linux; peak RSS remains available.
#endif
}

std::string csv( std::string const & value )
{
  std::string result = "\"";
  for( char c : value )
  {
    if( c == '"' )
      result += '"';
    result += c;
  }
  return result + '"';
}

class Reporter
{
public:
  Reporter( Communication & comm, Options const & opt ) : m_comm( comm ), m_opt( opt )
  {
    if( comm.rank() == 0 )
      std::cout << "label,input,kind,weak,ranks,px,py,pz,nx,ny,nz,levels,point_components,chunk_bytes,level,phase,"
                   "seconds_min,seconds_mean,seconds_max,cells_sum,cells_min,cells_max,points_sum,owned_points_sum,shared_points_sum,"
                   "incidence_sum,neighbors_max,rss_kib_max,peak_rss_kib_max,payload_bytes_sum,count_bytes_sum,chunks_sum,"
                   "directory_exchanges,neighbor_exchanges,discovery,sharing\n"
                << std::setprecision( 17 );
  }

  void row( int level, std::string const & phase, double seconds, Mesh const & mesh, CommunicationStatistics const & before )
  {
    auto reduce = [&]( auto value, MpiWrapper::Reduction op ) { return MpiWrapper::allReduce( value, op, MPI_COMM_GEOS ); };
    using R = MpiWrapper::Reduction;
    double const minimum = reduce( seconds, R::Min ), mean = reduce( seconds, R::Sum ) / m_comm.size(), maximum = reduce( seconds, R::Max );
    auto const cells = static_cast< std::uint64_t >( mesh.cells.size() ), points = static_cast< std::uint64_t >( mesh.coordinates.size() );
    std::uint64_t owned = mesh.sharingKnown ? points : 0, shared = 0;
    for( auto const & support : mesh.supports )
      if( support.key.kind == EntityKind::vertex )
      {
        if( mesh.sharingKnown )
        {
          owned -= support.participants.front() != m_comm.rank();
          shared += support.participants.size() > 1;
        }
      }
    auto const cellSum = reduce( cells, R::Sum ), cellMin = reduce( cells, R::Min ), cellMax = reduce( cells, R::Max );
    auto const pointSum = reduce( points, R::Sum ), ownedSum = reduce( owned, R::Sum ), sharedSum = reduce( shared, R::Sum );
    auto const incidenceSum = reduce( static_cast< std::uint64_t >( mesh.supports.size() ), R::Sum );
    auto const neighbors = reduce( static_cast< std::uint64_t >( m_comm.neighbors().size() ), R::Max );
    std::uint64_t currentMemory = 0, peakMemory = 0;
    m_comm.checked( "benchmark memory counters",
                    [&]
                    {
                      currentMemory = memoryKiB( false );
                      peakMemory = memoryKiB( true );
                    } );
    m_peakResident = std::max( { m_peakResident, currentMemory, peakMemory } );
    auto const rss = reduce( currentMemory, R::Max ), peak = reduce( m_peakResident, R::Max );
    auto const & stats = m_comm.statistics();
    auto const payload = reduce( stats.payloadBytesSent - before.payloadBytesSent, R::Sum );
    auto const counts = reduce( stats.countBytesSent - before.countBytesSent, R::Sum );
    auto const chunks = reduce( stats.payloadChunksSent - before.payloadChunksSent, R::Sum );
    auto const directory = reduce( stats.directoryExchanges - before.directoryExchanges, R::Max );
    auto const neighbor = reduce( stats.neighborExchanges - before.neighborExchanges, R::Max );
    if( m_comm.rank() == 0 )
      std::cout << csv( m_opt.label ) << ',' << csv( m_opt.input ) << ',' << m_opt.kind << ',' << m_opt.weak << ',' << m_comm.size() << ','
                << m_opt.grid[0] << ',' << m_opt.grid[1] << ',' << m_opt.grid[2] << ',' << m_opt.cells[0] << ',' << m_opt.cells[1] << ','
                << m_opt.cells[2] << ',' << m_opt.levels << ',' << m_opt.pointComponents << ',' << m_opt.chunkBytes << ',' << level << ','
                << phase << ',' << minimum << ',' << mean << ',' << maximum << ',' << cellSum << ',' << cellMin << ',' << cellMax << ','
                << pointSum << ',' << ownedSum << ',' << sharedSum << ',' << incidenceSum << ',' << neighbors << ',' << rss << ',' << peak
                << ',' << payload << ',' << counts << ',' << chunks << ',' << directory << ',' << neighbor << ',' << m_opt.discovery << ','
                << m_opt.sharing << '\n';
  }

  void phase( int level, std::string const & name, Mesh const & mesh, std::function< void() > const & work )
  {
    MpiWrapper::barrier( MPI_COMM_GEOS );
    auto const stats = m_comm.statistics();
    auto const start = Clock::now();
    work();
    row( level, name, std::chrono::duration< double >( Clock::now() - start ).count(), mesh, stats );
  }

private:
  Communication & m_comm;
  Options const & m_opt;
  std::uint64_t m_peakResident{};
};

struct ReadIds
{
  Connectivity & ids;
  template < typename Array > void operator()( Array * array ) const
  {
    vtkDataArrayAccessor< Array > access( array );
    using Value = typename decltype( access )::APIType;
    for( vtkIdType i = 0; i < array->GetNumberOfTuples(); ++i )
    {
      Value const value = access.Get( i, 0 );
      if constexpr( std::is_signed< Value >::value )
        if( value < 0 )
          throw std::invalid_argument( "Negative benchmark active ID" );
      if( static_cast< std::uintmax_t >( value ) > static_cast< std::uintmax_t >( std::numeric_limits< vtkIdType >::max() ) )
        throw std::overflow_error( "Benchmark active ID exceeds vtkIdType" );
      ids.push_back( static_cast< vtkIdType >( value ) );
    }
  }
};

Connectivity exactIds( vtkDataArray & data, vtkIdType count )
{
  if( data.GetNumberOfComponents() != 1 || data.GetNumberOfTuples() != count )
    throw std::invalid_argument( "Benchmark active IDs must be scalar arrays with matching tuple counts" );
  Connectivity result;
  result.reserve( count );
  ReadIds const worker{ result };
  if( !vtkArrayDispatch::DispatchByValueType< vtkArrayDispatch::Integrals >::Execute( &data, worker ) )
    throw std::invalid_argument( "Benchmark active IDs must use supported integral storage" );
  return result;
}

void run( Options const & opt, Communication & comm )
{
  Reporter reporter( comm, opt );
  Mesh mesh;
  vtkSmartPointer< vtkUnstructuredGrid > input = vtkSmartPointer< vtkUnstructuredGrid >::New();
  reporter.phase( 0, "input", mesh,
                  [&]
                  {
                    comm.checked( "benchmark coarse input",
                                  [&]
                                  {
                                    if( opt.input.empty() )
                                      mesh = generatedMesh( opt, comm.rank() );
                                    else if( comm.rank() == 0 )
                                    {
                                      vtkNew< vtkXMLUnstructuredGridReader > reader;
                                      reader->SetFileName( opt.input.c_str() );
                                      reader->Update();
                                      if( reader->GetErrorCode() || reader->GetOutput()->GetNumberOfCells() == 0 )
                                        throw std::runtime_error( "Cannot read nonempty benchmark VTU" );
                                      input->ShallowCopy( reader->GetOutput() );
                                      // The benchmark can seed absent IDs before shipping the coarse input.
                                      if( !input->GetPointData()->GetGlobalIds() )
                                      {
                                        vtkNew< vtkIdTypeArray > ids;
                                        ids->SetName( "benchmarkPointIds" );
                                        ids->SetNumberOfTuples( input->GetNumberOfPoints() );
                                        for( vtkIdType i = 0; i < input->GetNumberOfPoints(); ++i )
                                          ids->SetValue( i, i );
                                        input->GetPointData()->SetGlobalIds( ids );
                                      }
                                      if( !input->GetCellData()->GetGlobalIds() )
                                      {
                                        vtkNew< vtkIdTypeArray > ids;
                                        ids->SetName( "benchmarkCellIds" );
                                        ids->SetNumberOfTuples( input->GetNumberOfCells() );
                                        for( vtkIdType i = 0; i < input->GetNumberOfCells(); ++i )
                                          ids->SetValue( i, i );
                                        input->GetCellData()->SetGlobalIds( ids );
                                      }
                                    }
                                  } );
                  } );
  if( !opt.input.empty() )
    reporter.phase( 0, "scatterAndNormalize", mesh,
                    [&]
                    {
                      array1d< integer > partitions( 3 );
                      for( int d = 0; d < 3; ++d )
                        partitions[d] = opt.grid[d];
                      auto scattered =
                          geos::vtk::scatterMesh( geos::vtk::ScatterMethod::cartesian, *input, partitions.toViewConst(), MPI_COMM_GEOS );
                      comm.checked( "benchmark local input",
                                    [&]
                                    {
                                      auto * pointIds = scattered->GetPointData()->GetGlobalIds();
                                      auto * cellIds = scattered->GetCellData()->GetGlobalIds();
                                      if( !pointIds || !cellIds )
                                        throw std::invalid_argument( "Scatter dropped benchmark global IDs" );
                                      mesh.pointIds = exactIds( *pointIds, scattered->GetNumberOfPoints() );
                                      mesh.cellIds = exactIds( *cellIds, scattered->GetNumberOfCells() );
                                      for( vtkIdType i = 0; i < scattered->GetNumberOfPoints(); ++i )
                                      {
                                        Coordinates coordinate;
                                        scattered->GetPoint( i, coordinate.data() );
                                        mesh.coordinates.push_back( coordinate );
                                      }
                                      PointRegistry points( mesh.coordinates, mesh.pointIds );
                                      for( vtkIdType i = 0; i < scattered->GetNumberOfCells(); ++i )
                                      {
                                        mesh.cells.push_back( normalizeCoarseCell( *scattered->GetCell( i ), points ) );
                                      }
                                      mesh.pointData->ShallowCopy( scattered->GetPointData() );
                                      mesh.cellData->ShallowCopy( scattered->GetCellData() );
                                    } );
                      input = nullptr;
                    } );
  std::vector< EntityKey > discoveryKeys;
  std::vector< EntityKey > volumeKeys;
  reporter.phase( 0, opt.discovery == "boundary" ? "coarseBoundary" : "coarseAllIncidence", mesh,
                  [&]
                  {
                    comm.checked( "benchmark coarse incidence",
                                  [&]
                                  {
                                    if( opt.discovery == "boundary" )
                                    {
                                      auto boundary = coarseBoundary( mesh.cells, mesh.cellIds, mesh.pointIds, comm.rank() );
                                      volumeKeys = std::move( boundary.volumeIds );
                                      discoveryKeys = volumeKeys;
                                      for( auto const & entity : boundary.entities )
                                        discoveryKeys.push_back( entity.key );
                                      mesh.supports =
                                          opt.sharing == "all" ? incidence( mesh, comm.rank() ) : std::move( boundary.entities );
                                    }
                                    else
                                    {
                                      mesh.supports = incidence( mesh, comm.rank() );
                                      for( auto const & entity : mesh.supports )
                                        discoveryKeys.push_back( entity.key );
                                    }
                                    CellCounts counts;
                                    for( auto const & cell : mesh.cells )
                                      if( cell.prismSides )
                                        ++counts.prisms.at( cell.prismSides );
                                      else
                                        switch( cell.vtkType )
                                        {
                                        case VTK_HEXAHEDRON:
                                          ++counts.hexahedra;
                                          break;
                                        case VTK_TETRA:
                                          ++counts.tetrahedra;
                                          break;
                                        case VTK_PYRAMID:
                                          ++counts.pyramids;
                                          break;
                                        case VTK_WEDGE:
                                          ++counts.wedges;
                                          break;
                                        default:
                                          throw std::invalid_argument( "Unsupported benchmark cell" );
                                        }
                                    for( int level = 0; level < opt.levels && counts.total(); ++level )
                                    {
                                      counts = counts.next();
                                      if( counts.total() > static_cast< std::uint64_t >( std::numeric_limits< vtkIdType >::max() ) )
                                        throw std::overflow_error( "Benchmark refined cell count exceeds vtkIdType" );
                                    }
                                  } );
                  } );
  reporter.phase( 0, opt.discovery == "boundary" ? "coarseDiscoveryBoundaryAndVolumeIds" : "coarseDiscoveryAllEntities", mesh,
                  [&]
                  {
                    auto const sharing = comm.discoverSharing( discoveryKeys );
                    comm.checked( "benchmark discovery install",
                                  [&]
                                  {
                                    for( auto & entity : mesh.supports )
                                    {
                                      auto const found = sharing.find( entity.key );
                                      if( found != sharing.end() )
                                        entity.participants = found->second;
                                      if( entity.key.kind == EntityKind::cell && entity.participants.size() != 1 )
                                        throw std::invalid_argument( "Duplicate distributed volume owner" );
                                    }
                                    for( auto const & key : volumeKeys )
                                      if( sharing.at( key ).size() != 1 )
                                        throw std::invalid_argument( "Duplicate distributed volume owner" );
                                    if( opt.sharing == "interfaces" )
                                    {
                                      std::vector< EntitySupport > shared;
                                      for( auto & entity : mesh.supports )
                                        if( entity.participants.size() > 1 )
                                          shared.push_back( std::move( entity ) );
                                      mesh.supports = std::move( shared );
                                    }
                                    mesh.sharingKnown = true;
                                    std::vector< EntityKey >{}.swap( discoveryKeys );
                                    std::vector< EntityKey >{}.swap( volumeKeys );
                                  } );
                  } );
  TransferPolicies policies;
  // VTK extraction provenance is identifier metadata, not a categorical field.
  // The production orchestrator must rebuild lineage; this component harness
  // deliberately excludes it along with the active IDs and collocation arrays.
  policies.excludedPointArrays.insert( "vtkOriginalPointIds" );
  policies.excludedCellArrays.insert( "vtkOriginalCellIds" );
  reporter.phase( 0, "coarsePointFields", mesh,
                  [&]
                  {
                    std::vector< PointCreation > requests;
                    std::unique_ptr< PointFieldLayout > layout;
                    comm.checked( "benchmark coarse point fields",
                                  [&]
                                  {
                                    PointRegistry points( mesh.coordinates, mesh.pointIds );
                                    mesh.pointData = transferPointData( *mesh.pointData, points, policies );
                                    layout = std::make_unique< PointFieldLayout >( *mesh.pointData );
                                    for( auto const & entity : mesh.supports )
                                      if( entity.key.kind == EntityKind::vertex && entity.participants.size() > 1 )
                                      {
                                        vtkIdType const local = entity.localCorners.front();
                                        requests.push_back(
                                            { entity.key, entity.participants, mesh.coordinates[local], layout->pack( local ) } );
                                      }
                                  } );
                    auto const records = comm.reconcileExistingPoints( requests );
                    comm.checked( "benchmark coarse point install",
                                  [&]
                                  {
                                    for( auto const & entity : mesh.supports )
                                      if( entity.key.kind == EntityKind::vertex && entity.participants.size() > 1 )
                                      {
                                        auto const & record = records.at( entity.key );
                                        mesh.coordinates[entity.localCorners.front()] = record.position;
                                        layout->install( entity.localCorners.front(), record.fields );
                                      }
                                  } );
                  } );
  for( int level = 1; level <= opt.levels; ++level )
  {
    Mesh next;
    std::unique_ptr< PointRegistry > points;
    std::unique_ptr< SharingInheritance > inherited;
    std::unique_ptr< InterfaceSharing > interfaces;
    std::unique_ptr< PointFieldLayout > layout;
    Connectivity parents;
    std::vector< double > fractions;
    std::vector< PointCreation > creations;
    std::map< EntityKey, PointRecord > records;
    IdRange cells{};
    CommunicationStatistics const before = comm.statistics();
    double levelSeconds = 0;
    auto phase = [&]( std::string const & name, std::function< void() > const & work )
    {
      MpiWrapper::barrier( MPI_COMM_GEOS );
      auto const stats = comm.statistics();
      auto const start = Clock::now();
      work();
      double const seconds = std::chrono::duration< double >( Clock::now() - start ).count();
      levelSeconds += seconds;
      reporter.row( level, name, seconds, mesh, stats );
    };
    phase( "levelPlanAndBuildChildren",
           [&]
           {
             GEOS_MARK_SCOPE( "uniformRefinement/levelPlan" );
             comm.checked( "benchmark level planning",
                           [&]
                           {
                             points = std::make_unique< PointRegistry >( mesh.coordinates, mesh.pointIds );
                             for( std::size_t i = 0; i < mesh.cells.size(); ++i )
                             {
                               auto split = subdivideCell( mesh.cells[i], mesh.cellIds[i], *points );
                               double const measure = signedMeasure( mesh.cells[i], *points );
                               for( auto const & child : split.children )
                                 fractions.push_back( signedMeasure( child, *points ) / measure );
                               parents.insert( parents.end(), split.children.size(), i );
                               next.cells.insert( next.cells.end(), std::make_move_iterator( split.children.begin() ),
                                                  std::make_move_iterator( split.children.end() ) );
                             }
                             if( opt.sharing == "interfaces" )
                               interfaces = std::make_unique< InterfaceSharing >( *points, mesh.supports, comm.rank() );
                             else
                               inherited = std::make_unique< SharingInheritance >( *points, mesh.supports, comm.rank() );
                           } );
           } );
    phase( "transferPointFieldsAndPack",
           [&]
           {
             GEOS_MARK_SCOPE( "uniformRefinement/transferFields" );
             comm.checked( "benchmark point field transfer",
                           [&]
                           {
                             next.pointData = transferPointData( *mesh.pointData, *points, policies );
                             layout = std::make_unique< PointFieldLayout >( *next.pointData );
                             for( vtkIdType i = points->originalSize(); i < static_cast< vtkIdType >( points->points().size() ); ++i )
                             {
                               auto const & point = points->points()[i];
                               auto const & participants =
                                   interfaces ? interfaces->participants( { i } ) : inherited->participants( { i } );
                               creations.push_back( { point.key, participants, point.position,
                                                      participants.size() > 1 ? layout->pack( i, FieldTupleFormat::valuesOnly ) : Bytes{} } );
                             }
                           } );
           } );
    phase( "pointAndCellIdsAndExchange",
           [&]
           {
             vtkIdType const maximum = mesh.pointIds.empty() ? -1 : *std::max_element( mesh.pointIds.begin(), mesh.pointIds.end() );
             records = comm.resolvePoints( level, creations, maximum );
             cells = comm.allocateRange( next.cells.size(), 0 );
           } );
    phase( "installPointsAndTransferCellFields",
           [&]
           {
             GEOS_MARK_SCOPE( "uniformRefinement/buildChildren" );
             comm.checked( "benchmark level construction",
                           [&]
                           {
                             next.pointIds = mesh.pointIds;
                             next.coordinates.reserve( points->points().size() );
                             for( vtkIdType i = 0; i < static_cast< vtkIdType >( points->points().size() ); ++i )
                             {
                               auto const & point = points->points()[i];
                               if( i < points->originalSize() )
                                 next.coordinates.push_back( point.position );
                               else
                               {
                                 auto const & record = records.at( point.key );
                                 next.coordinates.push_back( record.position );
                                 next.pointIds.push_back( record.globalId );
                                 auto const & participants =
                                     interfaces ? interfaces->participants( { i } ) : inherited->participants( { i } );
                                 if( participants.size() > 1 )
                                   layout->install( i, record.fields, FieldTupleFormat::valuesOnly );
                               }
                             }
                             next.cellIds.reserve( next.cells.size() );
                             for( std::size_t i = 0; i < next.cells.size(); ++i )
                               next.cellIds.push_back( cells.first + i );
                             next.cellData = transferCellData( *mesh.cellData, mesh.cells.size(), parents, fractions, policies );
                           } );
           } );
    phase( opt.sharing == "interfaces" ? "fineInterfacesAndSharing" : "fineIncidenceAndSharing",
           [&]
           {
             comm.checked( "benchmark fine incidence",
                           [&]
                           {
                             if( interfaces )
                               next.supports = interfaces->fineSupports( next.pointIds );
                             else
                             {
                               next.supports = incidence( next, comm.rank() );
                               for( auto & entity : next.supports )
                                 if( entity.key.kind != EntityKind::cell )
                                   entity.participants = inherited->participants( entity.localCorners );
                             }
                             next.sharingKnown = true;
                           } );
           } );
    phase( "commitAndReleaseParent",
           [&]
           {
             mesh = std::move( next );
             inherited.reset();
             interfaces.reset();
             points.reset();
             layout.reset();
             std::vector< PointCreation >{}.swap( creations );
             std::map< EntityKey, PointRecord >{}.swap( records );
             Connectivity{}.swap( parents );
             std::vector< double >{}.swap( fractions );
           } );
    reporter.row( level, "levelTotal", levelSeconds, mesh, before );
  }
}
} // namespace

int main( int argc, char ** argv )
{
  if( argc == 2 && std::string( argv[1] ) == "--help" )
  {
    std::cout << "Usage: benchmarkUniformRefinement [--cells NX NY NZ] [--process-grid PX PY PZ] [--weak]\n"
                 "  [--kind hex|tet|pyramid] [--levels L] [--point-components N] [--chunk-bytes N] [--label NAME]\n"
                 "  [--input FILE.vtu]   CSV on rank zero. Input/scatter measured separately.\n"
                 "  [--discovery boundary|all] [--sharing interfaces|all]\n"
                 "  Component benchmark. Defaults: boundary geometry + volume ID discovery; interface-only sharing.\n";
    return 0;
  }
  MpiWrapper::init( &argc, &argv );
  MPI_COMM_GEOS = MpiWrapper::commDup( MPI_COMM_WORLD );
#ifdef GEOS_USE_CHAI
  // Match normal GEOS setup and keep allocation diagnostics out of CSV output.
  chai::ArrayManager::getInstance()->disableCallbacks();
#endif
  int result = 0;
  try
  {
    Communication comm( MPI_COMM_GEOS );
    Options opt;
    comm.checked( "benchmark arguments", [&] { opt = options( argc, argv, comm.size() ); } );
    // Construct with the requested transport limit after arguments are agreed.
    Communication benchmark( MPI_COMM_GEOS, opt.chunkBytes );
    run( opt, benchmark );
  }
  catch( std::exception const & error )
  {
    if( MpiWrapper::commRank( MPI_COMM_GEOS ) == 0 )
      std::cerr << error.what() << '\n';
    result = 1;
  }
  MpiWrapper::finalize();
  return result;
}
