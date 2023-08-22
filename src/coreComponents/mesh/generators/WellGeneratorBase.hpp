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
 * @file WellGeneratorBase.hpp
 *
 */

#ifndef GEOS_MESH_GENERATORS_WELLGENERATORBASE_HPP_
#define GEOS_MESH_GENERATORS_WELLGENERATORBASE_HPP_

#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"


namespace geos
{

/**
 * @class WellGeneratorBase
 *
 * This class processes the data of a single well from the XML and generates the well geometry
 */
class WellGeneratorBase : public dataRepository::Group
{
public:

  /**
   * @brief Constructor.
   * @param name name of the object in the data hierarchy.
   * @param parent pointer to the parent group in the data hierarchy.
   */
  WellGeneratorBase( const string & name,
                     Group * const parent );

  /**
   * @brief Get the catalog name.
   * @return the name of this type in the catalog
   */
  static string catalogName() { return "WellGeneratorBase"; }

  /// This function is used to expand any catalogs in the data structure
  virtual void expandObjectCatalogs() override;

  /// using alias for templated Catalog meshGenerator type
  using CatalogInterface = dataRepository::CatalogInterface< WellGeneratorBase, string const &, Group * const >;


  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Accessor for the singleton Catalog object
   * @return a static reference to the Catalog object
   */
  static CatalogInterface::CatalogType & getCatalog();

  /**
   * @brief Main function of the class that generates the well geometry
   */
  virtual void generateWellGeometry( ) = 0;

  /**
   * @name Getters / Setters
   */
  ///@{

  // getters for element data

  /**
   * @brief Get the global number of well elements.
   * @return the global number of elements
   */
  virtual globalIndex numElements() const = 0;

  /**
   * @brief Getter to the Segment to PolyNode mapping
   * @return The Segment to PolyNode mapping as a 2D array
   */
  virtual const array2d< globalIndex > & getSegmentToPolyNodeMap() const = 0;

  /**
   * @brief Get the number of nodes per well element
   * @return the number of nodes per well element
   */
  virtual globalIndex numNodesPerElement() const = 0;

  /**
   * @brief Get the Coordinates of the polyline nodes
   * @return the Coordinates of the polyline nodes
   */
  virtual const array2d< real64 > & getPolyNodeCoord() const = 0;

  /**
   * @return The minimum segment length
   */
  virtual real64 getMinSegmentLength() const = 0;

  /**
   * @return The minimum element length
   */
  virtual real64 getMinElemLength() const = 0;

  /**
   * @return The list of perforation names
   */
  const string_array & getPerforationList() const { return m_perforationList; }

  /**
   * @brief Get the physical location of the centers of well elements.
   * @return list of center locations of the well elements
   */
  virtual arrayView2d< real64 const > getElemCoords() const = 0;

  /**
   * @brief Get the global indices mapping an element to the next.
   * @return list providing the global index of the next element for each element
   */
  virtual arrayView1d< globalIndex const > getNextElemIndex() const = 0;

  /**
   * @brief Get the global indices mapping an element to the previous ones.
   * @return list providing the global indices of the previous elements for each element
   */
  virtual arrayView1d< arrayView1d< globalIndex const > const > getPrevElemIndices() const = 0;

  /**
   * @brief Get the global indices of the well nodes nodes connected to each element.
   * @return list providing the global index of the well nodes for each well element
   */
  virtual arrayView2d< globalIndex const > getElemToNodesMap() const = 0;

  /**
   * @brief Get the volume of the well elements.
   * @return list of volumes of the well elements
   */
  virtual arrayView1d< real64 const > getElemVolume() const = 0;

  /**
   * @brief Get the radius in the well.
   * @return the radius in the well
   */
  virtual real64 getElementRadius() const = 0;

  // getters for node data

  /**
   * @brief Get the global number of well nodes.
   * @return the global number of nodes
   */
  virtual globalIndex numNodes() const = 0;

  /**
   * @brief Get the physical location of the centers of well elements.
   * @return list of center locations of the well elements
   */
  virtual arrayView2d< real64 const > getNodeCoords() const = 0;



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
  virtual arrayView2d< real64 const > getPerfCoords() const = 0;

  /**
   * @brief Get the well transmissibility at the perforations.
   * @return list of well transmissibility at all the perforations on the well
   */
  virtual arrayView1d< real64 const > getPerfTransmissibility() const = 0;

  /**
   * @brief Get the global indices of the well elements connected to each perforation.
   * @return list providing the global index of the connected well element for each perforation
   */
  virtual arrayView1d< globalIndex const > getPerfElemIndex() const = 0;

  /**
   * @returns The number of physical dimensions
   */
  virtual int getPhysicalDimensionsNumber() const = 0;

  /**
   * Getter for the associated well region name
   * @return  the associated well region name
   */
  virtual const string getWellRegionName() const = 0;

  /**
   * Getter for the associated well control name
   * @return  the associated well control name
   */
  virtual const string getWellControlsName() const = 0;

  ///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    constexpr static char const * polylineNodeCoordsString() { return "polylineNodeCoords"; }
    constexpr static char const * polylineSegmentConnString() { return "polylineSegmentConn"; }
    constexpr static char const * numElementsPerSegmentString() { return "numElementsPerSegment"; }
    constexpr static char const * minSegmentLengthString() { return "minSegmentLength"; }
    constexpr static char const * minElementLengthString() { return "minElementLength"; }
    constexpr static char const * radiusString() { return "radius"; }
    constexpr static char const * wellRegionNameString() { return "wellRegionName"; }
    constexpr static char const * wellControlsNameString() { return "wellControlsName"; }
    constexpr static char const * meshNameString() { return "meshName"; }
    constexpr static char const * perforationString() { return "Perforation"; }
  };
  /// @endcond

protected:
  /// Global number of perforations
  globalIndex m_numPerforations;

  /// List of perforation names
  string_array m_perforationList;
};
}
#endif /* GEOS_MESH_GENERATORS_WELLGENERATORBASE_HPP_ */
