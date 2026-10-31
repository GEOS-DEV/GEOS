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
 * @file ConsistencyAdaptation.hpp
 */

#ifndef GEOS_MIXEDMIMETIC_CONSISTENCY_CONSISTENCYADAPTATION_HPP_
#define GEOS_MIXEDMIMETIC_CONSISTENCY_CONSISTENCYADAPTATION_HPP_

#include "common/DataTypes.hpp"
#include "mesh/ElementRegionManager.hpp"

#include <utility>

namespace geos
{

class MeshLevel;
class NeighborCommunicator;

/**
 * @class ConsistencyAdaptation
 * @brief Selects the inner product of every cell of the mixed mimetic discretization, eta = 1 for the
 *        consistent (MFD) product and eta = 0 for the diagonal (TPFA) product, through three layers
 *        applied in order: the residual-based consistency layer, the prescription read from the mesh
 *        and the degeneracy layer (admissibility of the consistent product). A prescribed cell is
 *        final: the two other layers act on the free cells only. The faces are then labelled: 0 when
 *        both cells use the diagonal product (the flux dof is condensed), 1 otherwise (saddle point).
 */
class ConsistencyAdaptation
{
public:

  /// Parameters of the three layers
  struct Parameters
  {
    bool adaptiveConsistency = true;            ///< run the consistency layer
    real64 consistencyTolerance = 1e-3;         ///< eta = 1 where the consistency indicator exceeds it
    real64 degeneracyTolerance = 0.1;           ///< eta = 0 below this percentage of the node-star volume
    real64 nominalGradient[3] = { 1.0, 1.0, 1.0 }; ///< gradient of the probe pressure field
    real64 lengthTolerance = 0.0;               ///< geometric tolerance of the face computations
    bool effectiveTpfa = false;                 ///< the selected product is itself TPFA: every face is condensed
  };

  /// Locally-owned cell counts of one classification (reduce over the ranks to report)
  struct Report
  {
    localIndex numCells = 0;        ///< cells of the target regions
    localIndex numConsistent = 0;   ///< cells with eta = 1 after the consistency layer
    localIndex numPrescribed0 = 0;  ///< cells prescribed the diagonal product
    localIndex numPrescribed1 = 0;  ///< cells prescribed the consistent product
    localIndex numDegenerate = 0;   ///< free cells switched to the diagonal product by the degeneracy layer
    localIndex numPrescribedDegenerate = 0; ///< cells prescribed the consistent product below the degeneracy tolerance (kept)
    localIndex numConsistentFinal = 0; ///< cells with eta = 1 after the three layers

    /**
     * @brief Accumulate the counts of another mesh level.
     * @param other the counts to add
     */
    void add( Report const & other )
    {
      numCells += other.numCells;
      numConsistent += other.numConsistent;
      numPrescribed0 += other.numPrescribed0;
      numPrescribed1 += other.numPrescribed1;
      numDegenerate += other.numDegenerate;
      numPrescribedDegenerate += other.numPrescribedDegenerate;
      numConsistentFinal += other.numConsistentFinal;
    }
  };

  /// Cell permeability by region and subregion
  using PermeabilityAccessor = ElementRegionManager::ElementViewConst< arrayView3d< real64 const > >;

  /**
   * @brief Run the three layers on the cells of the target regions and label the faces.
   *        Collective over the ranks: the fields are synchronized on the ghost cells.
   * @param mesh the mesh level
   * @param regionNames the target regions
   * @param regionFilter the target region indices
   * @param permeability the cell permeability
   * @param params the parameters of the layers
   * @param neighbors the neighbor communicators of the domain
   * @return the locally-owned cell counts
   */
  static Report classify( MeshLevel & mesh,
                          string_array const & regionNames,
                          SortedArrayView< localIndex const > const & regionFilter,
                          PermeabilityAccessor const & permeability,
                          Parameters const & params,
                          stdVector< NeighborCommunicator > & neighbors );

private:

  /**
   * @brief Consistency layer: eta = 0 where the two-point product reproduces the probe field within tolerance.
   * @return the locally-owned cells left with eta = 1
   */
  static localIndex applyConsistencyLayer( MeshLevel & mesh,
                                           string_array const & regionNames,
                                           SortedArrayView< localIndex const > const & regionFilter,
                                           PermeabilityAccessor const & permeability,
                                           Parameters const & params,
                                           stdVector< NeighborCommunicator > & neighbors );

  /**
   * @brief Prescription layer: eta follows prescribedMfdFlag where it is not negative.
   * @return the locally-owned cells prescribed 0 and prescribed 1
   */
  static std::pair< localIndex, localIndex > applyPrescription( MeshLevel & mesh,
                                                                string_array const & regionNames,
                                                                stdVector< NeighborCommunicator > & neighbors );

  /**
   * @brief Degeneracy layer on the free cells: eta = 0 where the cell volume is below the tolerance
   *        (percent of its node star); a prescribed cell is left unchanged.
   * @return the locally-owned free cells switched, and the prescribed eta = 1 cells below the tolerance
   */
  static std::pair< localIndex, localIndex > applyDegeneracyLayer( MeshLevel & mesh,
                                                                   string_array const & regionNames,
                                                                   real64 const tolerance,
                                                                   stdVector< NeighborCommunicator > & neighbors );

  /**
   * @brief Label the faces from eta: 0 when both cells use the diagonal product, 1 otherwise.
   */
  static void labelFaces( MeshLevel & mesh,
                          SortedArrayView< localIndex const > const & regionFilter,
                          bool const effectiveTpfa );
};

} // namespace geos

#endif /* GEOS_MIXEDMIMETIC_CONSISTENCY_CONSISTENCYADAPTATION_HPP_ */
