/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2020-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_WELLBLOCK_HPP
#define GEOSX_WELLBLOCK_HPP

#include "mesh/generators/WellBlockABC.hpp"
#include "mesh/generators/InternalWellGenerator.hpp"


namespace geosx
{

/**
 * Implementation of the WellBlock responsible for modification/creation capabilities.
 */
class WellBlock : public WellBlockABC
{
public:
  /**
   * @brief Constructor
   * @param internallWellGenerator
   */
  WellBlock( const InternalWellGenerator & internalWellGenrator );


  /**
   * @name Getters / Setters
   */
  ///@{

  // getters for element data

  /**
   * @brief Get the global number of well elements.
   * @return the global number of elements
   */
  globalIndex getNumElements() const { return m_numElems; }

  /**
   * @brief Get the physical location of the centers of well elements.
   * @return list of center locations of the well elements
   */
  arrayView2d< real64 const > getElemCoords() const { return m_elemCenterCoords; }

  /**
   * @brief Get the global indices mapping an element to the next.
   * @return list providing the global index of the next element for each element
   */
  arrayView1d< globalIndex const > getNextElemIndex() const { return m_nextElemId; }

  /**
   * @brief Get the global indices mapping an element to the previous ones.
   * @return list providing the global indices of the previous elements for each element
   */
  arrayView1d< arrayView1d< globalIndex const > const > getPrevElemIndices() const { return m_prevElemId.toNestedViewConst(); }

  /**
   * @brief Get the global indices of the well nodes nodes connected to each element.
   * @return list providing the global index of the well nodes for each well element
   */
  arrayView2d< globalIndex const > getElemToNodesMap() const { return m_elemToNodesMap; }

  /**
   * @brief Get the volume of the well elements.
   * @return list of volumes of the well elements
   */
  arrayView1d< real64 const > getElemVolume() const { return m_elemVolume; }

  /**
   * @brief Get the radius in the well.
   * @return the radius in the well
   */
  real64 getElementRadius() const { return m_radius; }

  // getters for node data

  /**
   * @brief Get the global number of well nodes.
   * @return the global number of nodes
   */
  globalIndex getNumNodes() const { return m_numNodes; }

  /**
   * @brief Get the physical location of the centers of well elements.
   * @return list of center locations of the well elements
   */
  arrayView2d< real64 const > getNodeCoords() const { return m_nodeCoords; }



  // getters for perforation data

  /**
   * @brief Get the global number of perforations on this well.
   * @return the global number of elements
   */
  globalIndex getNumPerforations() const { return m_numPerforations; }

  /**
   * @brief Get the locations of the perforations.
   * @return list of locations of all the perforations on the well
   */
  arrayView2d< real64 const > getPerfCoords() const { return m_perfCoords; }

  /**
   * @brief Get the well transmissibility at the perforations.
   * @return list of well transmissibility at all the perforations on the well
   */
  arrayView1d< real64 const > getPerfTransmissibility() const { return m_perfTransmissibility; }

  /**
   * @brief Get the global indices of the well elements connected to each perforation.
   * @return list providing the global index of the connected well element for each perforation
   */
  arrayView1d< globalIndex const > getPerfElemIndex() const { return m_perfElemId; }

  ///@}

  const string getWellRegionName() const { return m_wellRegionName; }
  const string getWellControlsName() const { return m_wellControlsName; }
  /// @endcond
  
  
protected:

  // /**
  //  * @brief This function provides capability to post process input values prior to
  //  * any other initialization operations.
  //  */
  // void postProcessInput() override final;

private:

  /**
   * @name Helper functions to construct the geometry of the well
   */
  ///@{

  /**
   * @brief Map each polyline node to the polyline segment(s) it is connected to.
   */
  void constructPolylineNodeToSegmentMap( const InternalWellGenerator & internalWellGenerator );

  /**
   * @brief Find the head node of the well (i.e., top node of the polyline).
   */
  void findPolylineHeadNodeIndex( const InternalWellGenerator & internalWellGenerator );

  /**
   * @brief Discretize the polyline by placing well elements.
   */
  void discretizePolyline( const InternalWellGenerator & internalWellGenerator );

  /**
   * @brief Map each perforation to a well element.
   */
  void connectPerforationsToWellElements( const InternalWellGenerator & internalWellGenerator );

  /**
   * @brief Make sure that the perforation locations are valid:
   *   - for partitioning purposes
   *   - to have a well-posed problem
   */
  void checkPerforationLocationsValidity( const InternalWellGenerator & internalWellGenerator );

  /**
   * @brief Merge perforations on the elements with multiple perforations.
   */
  void mergePerforations( array1d< array1d< localIndex > > const & elemToPerfMap, const InternalWellGenerator & internalWellGenerator );

  /**
   * @brief At a given node, find the next segment going in the direction of the bottom of the well.
   * @param[in] topSegId index of the top segment
   * @param[in] currentNodeId index of the current node
   */
  globalIndex getNextSegmentIndex( globalIndex topSegId,
                                   globalIndex currentNodeId ) const;

  ///@}

  /// @cond DO_NOT_DOCUMENT
  void debugWellGeometry( const InternalWellGenerator & internalWellGenerator ) const;
  /// @endcond

  // XML Input

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

  /// Name of the mesh body associated with this well
  string m_meshBodyName;



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

  /// Global number of perforations
  globalIndex m_numPerforations;

  /// Absolute physical location of the perforation
  array2d< real64 > m_perfCoords;

  /// Well Peaceman index at the perforation
  array1d< real64 > m_perfTransmissibility;

  /// Global index of the well element
  array1d< globalIndex > m_perfElemId;



  // Auxiliary data

  // Number of physical dimensions
  const int m_nDims;


  /// Map from the polyline nodes to the polyline nodes
  array1d< SortedArray< globalIndex > > m_polyNodeToSegmentMap;

  /// Index of the node at the well head
  globalIndex m_polylineHeadNodeId;

  /// Physical location of the polyline node wrt to well head
  array1d< real64 > m_nodeDistFromHead;

  // Perforation data

  /// List of perforation names
  string_array m_perforationList;

  /// Physical location of the perforation wrt to well head
  array1d< real64 > m_perfDistFromHead;

};
}
#endif
