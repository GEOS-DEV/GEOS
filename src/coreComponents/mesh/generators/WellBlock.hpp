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
  globalIndex numElements() const { return m_numElems; }

  /**
   * @brief Get the physical location of the centers of well elements.
   * @return list of center locations of the well elements
   */
  arrayView2d< real64 const > elemCoords() const { return m_elemCenterCoords; }

  /**
   * @brief Get the global indices mapping an element to the next.
   * @return list providing the global index of the next element for each element
   */
  arrayView1d< globalIndex const > nextElemIndex() const { return m_nextElemId; }

  /**
   * @brief Get the global indices mapping an element to the previous ones.
   * @return list providing the global indices of the previous elements for each element
   */
  arrayView1d< arrayView1d< globalIndex const > const > getPrevElemIndices() const { return m_prevElemId.toNestedViewConst(); }

  /**
   * @brief Get the global indices of the well nodes nodes connected to each element.
   * @return list providing the global index of the well nodes for each well element
   */
  arrayView2d< globalIndex const > elemToNodesMap() const { return m_elemToNodesMap; }

  /**
   * @brief Get the volume of the well elements.
   * @return list of volumes of the well elements
   */
  arrayView1d< real64 const > elemVolume() const { return m_elemVolume; }

  /**
   * @brief Get the radius in the well.
   * @return the radius in the well
   */
  real64 elementRadius() const { return m_radius; }

  // getters for node data

  /**
   * @brief Get the global number of well nodes.
   * @return the global number of nodes
   */
  globalIndex numNodes() const { return m_numNodes; }

  /**
   * @brief Get the physical location of the centers of well elements.
   * @return list of center locations of the well elements
   */
  arrayView2d< real64 const > nodeCoords() const { return m_nodeCoords; }



  // getters for perforation data

  /**
   * @brief Get the global number of perforations on this well.
   * @return the global number of elements
   */
  globalIndex numPerforations() const { return m_numPerforations; }

  /**
   * @brief Get the locations of the perforations.
   * @return list of locations of all the perforations on the well
   */
  arrayView2d< real64 const > perfCoords() const { return m_perfCoords; }

  /**
   * @brief Get the well transmissibility at the perforations.
   * @return list of well transmissibility at all the perforations on the well
   */
  arrayView1d< real64 const > perfTransmissibility() const { return m_perfTransmissibility; }

  /**
   * @brief Get the global indices of the well elements connected to each perforation.
   * @return list providing the global index of the connected well element for each perforation
   */
  arrayView1d< globalIndex const > perfElemIndex() const { return m_perfElemId; }

  ///@}

  const string getWellRegionName() const { return m_wellRegionName; }
  const string getWellControlsName() const { return m_wellControlsName; }
  /// @endcond
  
  

private:


  // XML Input

  /// Number of well elements per polyline interval
  int const m_numElemsPerSegment;

  /// Min segment length
  real64 const m_minSegmentLength;

  /// Min well element length
  real64 const m_minElemLength;

  /// Radius area of the well (assumed to be valid for the entire well)
  real64 const m_radius;

  /// Name of the corresponding well region
  string const m_wellRegionName;

  /// Name of the constraints associated with this well
  string const m_wellControlsName;



  // Geometry of the well (later passed to the WellElementSubRegion)

  // well element data

  /// Global number of well elements
  const globalIndex m_numElems;

  /// Physical location of the center of the well element
  arrayView2d< real64 const > m_elemCenterCoords;

  /// Global index of the next well element
  arrayView1d< globalIndex const > m_nextElemId;

  /// Global indices of the prev well elements (maybe need multiple prevs for branching)
  arrayView1d< arrayView1d< globalIndex const > const > m_prevElemId;

  /// Connectivity between elements and nodes
  arrayView2d< globalIndex const > m_elemToNodesMap;

  /// Volume of well elements
  arrayView1d< real64 const > m_elemVolume;


  // well node data

  /// Number of nodes per well element
  globalIndex const m_numNodesPerElem;

  /// Global number of well nodes
  globalIndex const m_numNodes;

  /// Physical location of the nodes
  arrayView2d< real64 const > m_nodeCoords;

  // perforation data

  /// Global number of perforations
  globalIndex const m_numPerforations;

  /// Absolute physical location of the perforation
  arrayView2d< real64 const > m_perfCoords;

  /// Well Peaceman index at the perforation
  arrayView1d< real64 const > m_perfTransmissibility;

  /// Global index of the well element
  arrayView1d< globalIndex const > m_perfElemId;



  // Auxiliary data

  // Number of physical dimensions
  const int m_nDims;

  // Perforation data

  /// List of perforation names
  const string_array m_perforationList;
};
}
#endif
