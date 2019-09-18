/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
string const nodeCoords       = "polylineNodeCoords";
string const segmentConn      = "polylineSegmentConn";
string const nElems           = "numElementsPerSegment";
string const crossSectionArea = "crossSectionArea";
string const wellRegionName   = "wellRegionName";
string const wellControlsName = "wellControlsName";
string const meshBodyName     = "meshName";
}
}

class WellGeneratorBase : public MeshGeneratorBase
{
  public:

  // define the top and bottom node of a segment
  struct NodeLocation
  {
    static constexpr integer TOP = 0;
    static constexpr integer BOTTOM = 1;
  };
  WellGeneratorBase( const std::string& name,
                     Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~WellGeneratorBase() override;

  /**
   * @return the name of this type in the catalog
   */  
  static string CatalogName() { return "WellGeneratorBase"; }

  /// not implemented
  virtual void GenerateElementRegions( DomainPartition& domain ) override {}

  virtual Group * CreateChild( string const & childKey, 
                               string const & childName ) override;


  /**
   * @brief main function of this class: processes the well input and creates the globla well topology
   * @param domain the physical domain object
   */  
  virtual void GenerateMesh( DomainPartition * const domain ) override = 0;

  /// not implemented 
  virtual void GetElemToNodesRelationInBox ( const std::string& elementType,
                                             const int index[],
                                             const int& iEle,
                                             int nodeIDInBox[],
                                             const int size) override {}

  /// not implemented
  virtual void RemapMesh ( dataRepository::Group * const domain ) override {}
     
  // getters for element data

  /**
   * @brief Getter for the global number of well elements
   * @return the global number of elements
   */
  globalIndex GetNumElements() const { return m_numElems; }

  /**
   * @brief Getter for the physical location of the centers of well elements 
   * @return list of center locations of the well elements
   */
  array1d<R1Tensor const> const & GetElemCoords() const { return m_elemCenterCoords; }

  /**
   * @brief Getter for the global indices mapping an element to the previous ones
   * @return list providing the global indices of the previous elements for each element
   */
  array1d< array1d<globalIndex> const> const & GetPrevElemIndices() const { return m_prevElemId; }

  /**
   * @brief Getter for the global indices mapping an element to the next
   * @return list providing the global index of the next element for each element
   */
  array1d<globalIndex const> const & GetNextElemIndex() const { return m_nextElemId; }

  /**
   * @brief Getter for the global indices of the well nodes nodes connected to each element
   * @return list providing the global index of the well nodes for each well element
   */
  array2d<globalIndex const> const & GetElemToNodesMap() const { return m_elemToNodesMap; }

  /**
   * @brief Getter for the physical location of the centers of well elements 
   * @return list of center locations of the well elements
   */
  array1d<R1Tensor const> const & GetNodeCoords() const { return m_nodeCoords; }

  /**
   * @brief Getter for the volume of the well elements 
   * @return list of volumes of the well elements
   */
  array1d<real64 const> const & GetElemVolume() const { return m_elemVolume; }

  // getters for perforation data

  /**
   * @brief Getter for the transmissibility at the perforations
   * @return list of transmissibilities at all the perforations on the well
   */
  array1d<real64 const> const & GetPerfTransmissibility() const { return m_perfTrans; }
  
  /**
   * @brief Getter for the locations of the perforations
   * @return list of locations of all the perforations on the well
   */
  array1d<R1Tensor const> const & GetPerfCoords() const { return m_perfCoords; }

  /**
   * @brief Getter for the global indices of the well elements connected to each perforation
   * @return list providing the global index of the connected well element for each perforation
   */
  array1d<globalIndex const> const & GetPerfElemIndex() const { return m_perfElemId; }

  /**
   * @brief Getter for the global number of perforations on this well
   * @return the global number of elements
   */
  globalIndex GetNumPerforations() const { return m_numPerforations; }

  /**
   * @brief Getter for the global number of well nodes
   * @return the global number of nodes
   */
  globalIndex GetNumNodes() const { return m_numNodes; }

  void DebugWellGeometry() const;

  protected:
  void PostProcessInput() override;

  protected:

  // Basic information on the well, given by the input
  
  /// Coordinates of the polyline
  array2d<real64>      m_inputPolyNodeCoords;

  /// Connectivity between the polyline nodes
  array2d<globalIndex> m_segmentToPolyNodeMap;

  /// Number of well elements per polyline interval
  int m_numElemsPerSegment;

  /// Cross section area of the well (assumed to be valid for the entire well)
  real64 m_crossSectionArea;

  /// Name of the corresponding well region
  string m_wellRegionName;

  /// Name of the constraints associated with this well
  string m_wellControlsName;

  /// Name of the mesh body associated with this well
  string m_meshBodyName;
  
  // Geometry of the well (later passed to the WellElementSubRegion)

  // well element data

  /// Global number of well elements
  globalIndex          m_numElems;
  /// Physical location of the center of the well element
  array1d<R1Tensor>    m_elemCenterCoords;

  /// Global indices of the prev well elements (maybe need multiple prevs for branching)
  array1d< array1d<globalIndex> > m_prevElemId;

  /// Global index of the next well element
  array1d<globalIndex> m_nextElemId;

  /// Connectivity between elements and nodes
  array2d<globalIndex> m_elemToNodesMap;

  /// Volume of well elements
  array1d<real64> m_elemVolume;



  /// Number of nodes per well element
  globalIndex const    m_numNodesPerElem;

  /// Global number of well nodes
  globalIndex          m_numNodes;
 
  /// Global number of perforations
  globalIndex          m_numPerforations;
 
  /// Physical location of the nodes
  array1d<R1Tensor>    m_nodeCoords;

  /// Absolute physical location of the perforation 
  array1d<R1Tensor>    m_perfCoords;

  /// Transmissibility at the perforation
  array1d<real64>      m_perfTrans;

  /// Global index of the well element
  array1d<globalIndex> m_perfElemId;


  // Other data

  // Number of physical dimensions
  const int m_nDims;
  
  /// Coordinates of the polyline nodes in R1Tensor format
  array1d<R1Tensor>           m_polyNodeCoords;
  
  /// Map from the polyline nodes to the polyline nodes
  array1d< set<globalIndex> > m_polyNodeToSegmentMap;

  /// Index of the node at the well head
  globalIndex     m_polylineHeadNodeId;

  /// Physical location of the polyline node wrt to well head
  array1d<real64> m_nodeDistFromHead;

  // Perforation data

  /// List of perforation names
  string_array    m_perforationList;

  /// Physical location of the perforation wrt to well head
  array1d<real64> m_perfDistFromHead;

};
} // namespace

#endif
