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

/*
 * @file InternalWellGenerator.hpp
 *
 */

#ifndef GEOS_MESH_GENERATORS_INTERNALWELLGENERATOR_HPP_
#define GEOS_MESH_GENERATORS_INTERNALWELLGENERATOR_HPP_

#include "WellGeneratorBase.hpp"

#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"


namespace geos
{

/**
 * @class InternalWellGenerator
 *
 * This class processes the data of a single well from the XML and generates the well geometry
 */
class InternalWellGenerator : public WellGeneratorBase
{
public:


  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param name name of the object in the data hierarchy.
   * @param parent pointer to the parent group in the data hierarchy.
   */
  InternalWellGenerator( const string & name,
                         Group * const parent );

  /**
   * @brief Get the catalog name.
   * @return the name of this type in the catalog
   */
  static string catalogName() { return "InternalWell"; }


  /**
   * @brief Main function of the class that generates the well geometry
   */
  void generateWellGeometry( ) override;


  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  // getters for element data

  /**
   * @brief Get the global number of well elements.
   * @return the global number of elements
   */
  globalIndex numElements() const override { return m_numElems; }

  /**
   * @brief Getter to the Segment to PolyNode mapping
   * @return The Segment to PolyNode mapping as a 2D array
   */
  const array2d< globalIndex > & getSegmentToPolyNodeMap() const override { return m_segmentToPolyNodeMap; };

  /**
   * @brief Get the number of nodes per well element
   * @return the number of nodes per well element
   */
  globalIndex numNodesPerElement() const override { return m_numNodesPerElem; }

  /**
   * @brief Get the Coordinates of the polyline nodes
   * @return the Coordinates of the polyline nodes
   */
  const array2d< real64 > & getPolyNodeCoord() const override { return m_polyNodeCoords; }

  /**
   * @return The minimum segment length
   */
  real64 getMinSegmentLength() const override { return m_minSegmentLength; }

  /**
   * @return The minimum element length
   */
  real64 getMinElemLength() const override { return m_minElemLength; }

  /**
   * @brief Get the physical location of the centers of well elements.
   * @return list of center locations of the well elements
   */
  arrayView2d< real64 const > getElemCoords() const override { return m_elemCenterCoords; }

  /**
   * @brief Get the global indices mapping an element to the next.
   * @return list providing the global index of the next element for each element
   */
  arrayView1d< globalIndex const > getNextElemIndex() const override { return m_nextElemId; }

  /**
   * @brief Get the global indices mapping an element to the previous ones.
   * @return list providing the global indices of the previous elements for each element
   */
  arrayView1d< arrayView1d< globalIndex const > const > getPrevElemIndices() const override { return m_prevElemId.toNestedViewConst(); }

  /**
   * @brief Get the global indices of the well nodes nodes connected to each element.
   * @return list providing the global index of the well nodes for each well element
   */
  arrayView2d< globalIndex const > getElemToNodesMap() const override { return m_elemToNodesMap; }

  /**
   * @brief Get the volume of the well elements.
   * @return list of volumes of the well elements
   */
  arrayView1d< real64 const > getElemVolume() const override { return m_elemVolume; }

  /**
   * @brief Get the radius in the well.
   * @return the radius in the well
   */
  real64 getElementRadius() const override { return m_radius; }

  // getters for node data

  /**
   * @brief Get the global number of well nodes.
   * @return the global number of nodes
   */
  globalIndex numNodes() const override { return m_numNodes; }

  /**
   * @brief Get the physical location of the centers of well elements.
   * @return list of center locations of the well elements
   */
  arrayView2d< real64 const > getNodeCoords() const override { return m_nodeCoords; }

  // getters for perforation data

  /**
   * @brief Get the locations of the perforations.
   * @return list of locations of all the perforations on the well
   */
  arrayView2d< real64 const > getPerfCoords() const override { return m_perfCoords; }

  /**
   * @brief Get the well transmissibility at the perforations.
   * @return list of well transmissibility at all the perforations on the well
   */
  arrayView1d< real64 const > getPerfTransmissibility() const override { return m_perfTransmissibility; }

  /**
   * @brief Get the global indices of the well elements connected to each perforation.
   * @return list providing the global index of the connected well element for each perforation
   */
  arrayView1d< globalIndex const > getPerfElemIndex() const override { return m_perfElemId; }

  /**
   * @returns The number of physical dimensions
   */
  int getPhysicalDimensionsNumber() const override { return m_nDims; }

  ///@}

  const string getWellRegionName() const override { return m_wellRegionName; }
  const string getWellControlsName() const override { return m_wellControlsName; }



protected:

  /**
   * @brief This function provides capability to post process input values prior to
   * any other initialization operations.
   */
  void postProcessInput() override final;

private:

  /**
   * @name Helper functions to construct the geometry of the well
   */
  ///@{

  /**
   * @brief Map each polyline node to the polyline segment(s) it is connected to.
   */
  void constructPolylineNodeToSegmentMap();

  /**
   * @brief Find the head node of the well (i.e., top node of the polyline).
   */
  void findPolylineHeadNodeIndex();

  /**
   * @brief Discretize the polyline by placing well elements.
   */
  void discretizePolyline();

  /**
   * @brief Map each perforation to a well element.
   */
  void connectPerforationsToWellElements();

  /**
   * @brief Make sure that the perforation locations are valid:
   *   - for partitioning purposes
   *   - to have a well-posed problem
   */
  void checkPerforationLocationsValidity();

  /**
   * @brief Merge perforations on the elements with multiple perforations.
   */
  void mergePerforations( array1d< array1d< localIndex > > const & elemToPerfMap );

  /**
   * @brief At a given node, find the next segment going in the direction of the bottom of the well.
   * @param[in] topSegId index of the top segment
   * @param[in] currentNodeId index of the current node
   */
  globalIndex getNextSegmentIndex( globalIndex topSegId,
                                   globalIndex currentNodeId ) const;

  ///@}

  /// @cond DO_NOT_DOCUMENT
  void debugWellGeometry() const;
  /// @endcond

  // XML Input

  /// Connectivity between the polyline nodes
  array2d< globalIndex > m_segmentToPolyNodeMap;

  /// Number of well elements per polyline interval
  int m_numElemsPerSegment;

  /// Min segment length
  real64 m_minSegmentLength;

  /// Min well element length
  real64 m_minElemLength;

  /// Radius area of the well (assumed to be valid for the entire well)
  real64 m_radius;

  /// Name of the corresponding well region
  string m_wellRegionName;

  /// Name of the constraints associated with this well
  string m_wellControlsName;



  // Geometry of the well (later passed to the WellElementSubRegion)

  // well element data

  /// Global number of well elements
  globalIndex m_numElems;

  /// Physical location of the center of the well element
  array2d< real64 > m_elemCenterCoords;

  /// Global index of the next well element
  array1d< globalIndex > m_nextElemId;

  /// Global indices of the prev well elements (maybe need multiple prevs for branching)
  array1d< array1d< globalIndex > > m_prevElemId;

  /// Connectivity between elements and nodes
  array2d< globalIndex > m_elemToNodesMap;

  /// Volume of well elements
  array1d< real64 > m_elemVolume;


  // well node data

  /// Number of nodes per well element
  globalIndex const m_numNodesPerElem;

  /// Global number of well nodes
  globalIndex m_numNodes;

  /// Physical location of the nodes
  array2d< real64 > m_nodeCoords;

  // perforation data

  /// Absolute physical location of the perforation
  array2d< real64 > m_perfCoords;

  /// Well Peaceman index at the perforation
  array1d< real64 > m_perfTransmissibility;

  /// Global index of the well element
  array1d< globalIndex > m_perfElemId;



  // Auxiliary data

  // Number of physical dimensions
  const int m_nDims;

  /// Coordinates of the polyline nodes
  array2d< real64 > m_polyNodeCoords;

  /// Map from the polyline nodes to the polyline nodes
  array1d< SortedArray< globalIndex > > m_polyNodeToSegmentMap;

  /// Index of the node at the well head
  globalIndex m_polylineHeadNodeId;

  /// Physical location of the polyline node wrt to well head
  array1d< real64 > m_nodeDistFromHead;

  // Perforation data

  /// Physical location of the perforation wrt to well head
  array1d< real64 > m_perfDistFromHead;

};
}
#endif /* GEOS_MESH_GENERATORS_INTERNALWELLGENERATOR_HPP_ */
