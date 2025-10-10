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
 * @class LayeredMeshPartitioner
 */
class LayeredMeshPartitioner : public MeshPartitioner
{
public:

  /**
   * @brief Constructor
   * @param name The name of this partitioner instance
   * @param parent The parent group
   */
  explicit LayeredMeshPartitioner( string const & name,
                                   dataRepository::Group * const parent );

  /**
   * @brief Destructor
   */
  virtual ~LayeredMeshPartitioner() override;

  /**
   * @brief Catalog name for factory registration
   * @return The catalog key
   */
  static string catalogName() { return "LayeredMeshPartitioner"; }

  /**
   * @brief Get the catalog name
   * @return The catalog key
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /**
   * @brief Process command-line partition count overrides
   *
   * For layered partitioner, Z-direction override sets numPartZ.
   *
   * @param xparCL X-direction partition count (0 = no override)
   * @param yparCL Y-direction partition count (0 = no override)
   * @param zparCL Z-direction partition count (0 = no override, sets numPartZ)
   */
  virtual void processCommandLineOverrides( unsigned int const xparCL,
                                            unsigned int const yparCL,
                                            unsigned int const zparCL ) override;

  /**
   * @brief Struct to serve as a container for variable strings and keys
   */
  struct viewKeyStruct : MeshPartitioner::viewKeyStruct
  {
    /// String key for the structured index array name
    static constexpr char const * indexArrayNameString() { return "indexArrayName"; }

    /// String key for the number of partitions in Z direction
    static constexpr char const * numPartZString() { return "numPartZ"; }
  };

protected:

  /**
   * @brief Post-processing after input has been read
   */
  virtual void postInputInitialization() override;

  /**
   * @brief Compute partition assignment using area-layer algorithm
   *
   * @param mesh Input VTK mesh with structured index array
   * @param comm MPI communicator
   * @return Partition assignment for each element (main mesh only)
   */
  virtual array1d< int64_t > computePartitioning( vtk::AllMeshes & mesh,
                                                  MPI_Comm const comm ) override;

private:

  /// Name of VTK cell data array containing [area, layer] indices
  string m_indexArrayName;

  /// Number of partitions in layer Z-direction
  int m_numPartZ;
};
