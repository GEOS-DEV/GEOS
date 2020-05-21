/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_CORECOMPONENTS_WELLS_WELLGENERATORBASE_HPP_
#define GEOSX_CORECOMPONENTS_WELLS_WELLGENERATORBASE_HPP_


#include "dataRepository/Group.hpp"
#include "meshUtilities/MeshGeneratorBase.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const nElems           = "numElementsPerSegment";
string const radius           = "radius";
string const wellRegionName   = "wellRegionName";
string const wellControlsName = "wellControlsName";
string const meshBodyName     = "meshName";
}
}

class WellGeneratorBase : public MeshGeneratorBase
{
public:

  /**
   * @brief Struct to define the top and bottom node of a segment.
   */
  struct NodeLocation
  {
    static constexpr integer TOP = 0;
    static constexpr integer BOTTOM = 1;
  };

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param name name of the object in the data hierarchy.
   * @param parent pointer to the parent group in the data hierarchy.
   */
  WellGeneratorBase( const std::string & name,
                     Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~WellGeneratorBase() override;
  //
  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this type in the catalog
   */
  static string CatalogName() { return "WellGeneratorBase"; }
  ///@}

  /**
   * @name Overriding functions defined in MeshGeneratorBase and above
   */
  ///@{

  /**
   * @brief Expand any catalogs in the data structure.
   */
  virtual void ExpandObjectCatalogs() override;
  
  /// not implemented
  virtual void GenerateElementRegions( DomainPartition & GEOSX_UNUSED_PARAM( domain ) ) override {}

  /**
   * @brief Creates a new sub-Group using the ObjectCatalog functionality.
   * @param[in] childKey The name of the new object type's key in the
   *                     ObjectCatalog.
   * @param[in] childName The name of the new object in the collection of
   *                      sub-Groups.
   * @return A pointer to the new Group created by this function.
   */
  virtual Group * CreateChild( string const & childKey,
                               string const & childName ) override;

  /**
   * @brief Main function of the class that generates the well geometry
   * @param[in] domain the domain object
   */
  virtual void GenerateMesh( DomainPartition * const domain ) override final;

  /// not implemented
  virtual void GetElemToNodesRelationInBox ( std::string const & GEOSX_UNUSED_PARAM( elementType ),
                                             int const * GEOSX_UNUSED_PARAM( index ),
                                             int const & GEOSX_UNUSED_PARAM( iEle ),
                                             int * GEOSX_UNUSED_PARAM( nodeIDInBox ),
                                             int const GEOSX_UNUSED_PARAM( size )) override {}


  /// not implemented
  virtual void RemapMesh ( dataRepository::Group * const GEOSX_UNUSED_PARAM( domain ) ) override {}

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Get the global number of well elements.
   * @return the global number of elements
   */
  globalIndex GetNumElements() const { return m_numElems; }

  /**
   * @brief Get the physical location of the centers of well elements.
   * @return list of center locations of the well elements
   */
  arrayView1d< R1Tensor const > const & GetElemCoords() const { return m_elemCenterCoords; }

  /**
   * @brief Get the global indices mapping an element to the previous ones.
   * @return list providing the global indices of the previous elements for each element
   */
  arrayView1d< arrayView1d< globalIndex const > const > const & GetPrevElemIndices() const { return m_prevElemId.toViewConst(); }

  /**
   * @brief Get the global indices mapping an element to the next.
   * @return list providing the global index of the next element for each element
   */
  arrayView1d< globalIndex const > const & GetNextElemIndex() const { return m_nextElemId; }

  /**
   * @brief Get the global indices of the well nodes nodes connected to each element.
   * @return list providing the global index of the well nodes for each well element
   */
  arrayView2d< globalIndex const > const & GetElemToNodesMap() const { return m_elemToNodesMap; }

  /**
   * @brief Get the physical location of the centers of well elements.
   * @return list of center locations of the well elements
   */
  arrayView1d< R1Tensor const > const & GetNodeCoords() const { return m_nodeCoords; }

  /**
   * @brief Get the volume of the well elements.
   * @return list of volumes of the well elements
   */
  arrayView1d< real64 const > const & GetElemVolume() const { return m_elemVolume; }

  // getters for perforation data

  /**
   * @brief Get the well transmissibility at the perforations.
   * @return list of well transmissibility at all the perforations on the well
   */
  arrayView1d< real64 const > const & GetPerfTransmissibility() const { return m_perfTransmissibility; }

  /**
   * @brief Get the locations of the perforations.
   * @return list of locations of all the perforations on the well
   */
  arrayView1d< R1Tensor const > const & GetPerfCoords() const { return m_perfCoords; }

  /**
   * @brief Get the global indices of the well elements connected to each perforation.
   * @return list providing the global index of the connected well element for each perforation
   */
  arrayView1d< globalIndex const > const & GetPerfElemIndex() const { return m_perfElemId; }

  /**
   * @brief Get the radius in the well.
   * @return the radius in the well
   */
  real64 GetElementRadius() const { return m_radius; }

  /**
   * @brief Get the global number of perforations on this well.
   * @return the global number of elements
   */
  globalIndex GetNumPerforations() const { return m_numPerforations; }

  /**
   * @brief Get the global number of well nodes.
   * @return the global number of nodes
   */
  globalIndex GetNumNodes() const { return m_numNodes; }
  
  ///@}

  /// @cond DO_NOT_DOCUMENT
  void DebugWellGeometry() const;
  /// @endcond

protected:
  /**
   * @brief This function provides capability to post process input values prior to
   * any other initialization operations.
   */
  void PostProcessInput() override;

  virtual void GeneratePolyLine() = 0;

  /**
   * @name Helper functions to construct the geometry of the well
   */
  ///@{
  
  /**
   * @brief Map each polyline node to the polyline segment(s) it is connected to.
   */
  void ConstructPolylineNodeToSegmentMap();

  /**
   * @brief Find the head node of the well (i.e., top node of the polyline).
   */
  void FindPolylineHeadNodeIndex();

  /**
   * @brief Discretize the polyline by placing well elements.
   */
  void DiscretizePolyline();

  /**
   * @brief Map each perforation to a well element.
   */
  void ConnectPerforationsToWellElements();

  /**
   * @brief Merge perforations on the elements with multiple perforations.
   */
  void MergePerforations();

  /**
   * @brief At a given node, find the next segment going in the direction of the bottom of the well
   */
  globalIndex GetNextSegmentIndex( globalIndex topSegId,
                                   globalIndex currentNodeId ) const;
  
  ///@}

protected:

  // Basic information on the well, given by the input

  /// Coordinates of the polyline
  array2d< real64 >      m_inputPolyNodeCoords;

  /// Connectivity between the polyline nodes
  array2d< globalIndex > m_segmentToPolyNodeMap;

  /// Number of well elements per polyline interval
  int m_numElemsPerSegment;

  /// Radius of the well (assumed to be valid for the entire well)
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
  array1d< R1Tensor >    m_elemCenterCoords;

  /// Global indices of the prev well elements (maybe need multiple prevs for branching)
  array1d< array1d< globalIndex > > m_prevElemId;

  /// Global index of the next well element
  array1d< globalIndex > m_nextElemId;

  /// Connectivity between elements and nodes
  array2d< globalIndex > m_elemToNodesMap;

  /// Volume of well elements
  array1d< real64 > m_elemVolume;



  /// Number of nodes per well element
  globalIndex const m_numNodesPerElem;

  /// Global number of well nodes
  globalIndex m_numNodes;

  /// Global number of perforations
  globalIndex m_numPerforations;

  /// Physical location of the nodes
  array1d< R1Tensor >    m_nodeCoords;

  /// Absolute physical location of the perforation
  array1d< R1Tensor >    m_perfCoords;

  /// Well Peaceman index at the perforation
  array1d< real64 >      m_perfTransmissibility;

  /// Global index of the well element
  array1d< globalIndex > m_perfElemId;


  // Other data

  // Number of physical dimensions
  const int m_nDims;

  /// Coordinates of the polyline nodes in R1Tensor format
  array1d< R1Tensor >           m_polyNodeCoords;

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
} // namespace

#endif
