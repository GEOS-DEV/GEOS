#include "CellGraphPartitioner.hpp"
#include "mesh/generators/VTKUtilities.hpp"
#include "GraphPartitionEngine.hpp"

namespace geos
{
using namespace dataRepository;

CellGraphPartitioner::CellGraphPartitioner( string const & name, Group * const parent )
  : MeshPartitioner( name, parent )
{}

CellGraphPartitioner::~CellGraphPartitioner() = default;

array1d< int64_t > CellGraphPartitioner::computePartitioning( vtk::AllMeshes & mesh,
                                                              MPI_Comm const comm )
{
  GEOS_MARK_FUNCTION;

  pmet_idx_t const numParts = MpiWrapper::commSize( comm );
  return vtk::partitionByCellGraph( mesh, getEngine(), comm, numParts );
}

string CellGraphPartitioner::getInfoString() const
{
  string const engineName = getEngine().getCatalogName();

  return GEOS_FMT( "{} '{}': engine={}, refinements={}, numPartitions={}, numNeighbors={}, numColors={}, color={}",
                   catalogName(),
                   getName(),
                   engineName,
                   m_numRefinements,
                   m_numPartitions,
                   m_neighborsRank.size(),
                   m_numColors,
                   m_color );
}

// Register in DomainPartitioner catalog
REGISTER_CATALOG_ENTRY( DomainPartitioner, CellGraphPartitioner, string const &, Group * const )

} // namespace geos
