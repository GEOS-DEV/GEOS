#ifndef GEOS_MESH_MPICOMMUNICATIONS_MESHPARTITIONER_HPP
#define GEOS_MESH_MPICOMMUNICATIONS_MESHPARTITIONER_HPP

#include "DomainPartitioner.hpp"
#include "mesh/generators/VTKUtilities.hpp"

#include <memory>

namespace geos
{

class GraphPartitionEngine;

/**
 * @class MeshPartitioner
 * @brief Base class for mesh-based partitioners that use graph partitioning
 */
class MeshPartitioner : public DomainPartitioner
{
public:

  MeshPartitioner( string const & name, Group * const parent );

  virtual ~MeshPartitioner() = default;

  /**
   * @brief Partition and redistribute a mesh
   *
   * This is the main workflow method that:
   * 1. Calls computePartitioning() to get partition assignments (virtual, subclass-specific)
   * 2. Redistributes the mesh using VTK utilities
   * 3. Extracts neighbors from the redistributed mesh
   * 4. Computes graph coloring for communication scheduling
   *
   * @param mesh Input mesh (loaded from file or generated)
   * @param comm MPI communicator
   * @return Redistributed mesh
   *
   * @post m_neighborsRank is set
   * @post m_color and m_numColors are set
   */
  virtual vtk::AllMeshes partitionMeshes( vtk::AllMeshes & mesh, MPI_Comm const comm );

  /**
   * @brief Post-input initialization - creates default engine if needed
   */
  virtual void postInputInitialization() override;


  /**
   * @brief Set the graph partition engine
   *
   * @param engine Reference to the engine (ParMETIS or PTScotch)
   *
   * @note The engine must outlive this partitioner (typically owned by parent)
   */
  void setEngine( GraphPartitionEngine & engine );

  /**
   * @brief Get the graph partition engine
   *
   * @return Reference to the engine
   * @throws If engine has not been set
   */
  GraphPartitionEngine & getEngine() const;

  /**
   * @brief Process command-line overrides for partition counts
   *
   * @param xparCL X-direction partition count from command line (0 = not specified)
   * @param yparCL Y-direction partition count from command line (0 = not specified)
   * @param zparCL Z-direction partition count from command line (0 = not specified)
   */
  void processCommandLineOverrides( unsigned int const xparCL,
                                    unsigned int const yparCL,
                                    unsigned int const zparCL ) override;

  /**
   * @brief Compute graph coloring for communication scheduling
   */
  void color() override;

  /**
   * @brief Struct to serve as a container for variable strings and keys
   */
  struct viewKeyStruct : DomainPartitioner::viewKeyStruct
  {
    static constexpr char const * engineTypeString() { return "engine"; }
    static constexpr char const * numRefinementsString() { return "numRefinements"; }
  };

  /// @brief ViewKeys for MeshPartitioner
  struct viewKeysStruct
  {
    static constexpr char const * engineString() { return "engine"; }
    dataRepository::ViewKey engine = { engineString() };
  } viewKeys;

protected:

  /**
   * @brief Compute partition assignment for mesh elements
   *
   * This is the algorithm-specific part that subclasses must implement.
   * Different partitioners use different algorithms:
   * - CellGraphPartitioner: Pure dual graph
   * - StructuredMeshPartitioner: Area-layer decomposition
   *
   * @param mesh Input mesh
   * @param comm MPI communicator
   * @return Partition assignment for each element
   */
  virtual array1d< int64_t > computePartitioning( vtk::AllMeshes & mesh,
                                                  MPI_Comm const comm ) = 0;

  /**
   * @brief Extract neighbor ranks from redistributed mesh
   *
   * @param mesh Redistributed mesh
   * @param comm MPI communicator
   */
  void extractNeighborsFromMesh( vtk::AllMeshes const & mesh, MPI_Comm const comm );

  /**
   * @brief Creates a default graph partition engine appropriate for the current MPI configuration
   * @return Pointer to the created engine
   */
  GraphPartitionEngine * createDefaultEngine();

  /// Graph partitioning engine (NOT owned)
  GraphPartitionEngine * m_engine;

  /// Engine type ("parmetis", "ptscotch", or "noop")
  string m_engineType;

  /// Number of refinement iterations (ParMETIS only, ignored by others)
  int m_numRefinements;
};

} // namespace geos

#endif
