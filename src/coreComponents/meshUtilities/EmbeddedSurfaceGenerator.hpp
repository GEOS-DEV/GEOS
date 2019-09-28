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

/**
 * @file EmbeddedSurfaceGenerator.hpp
 */
#ifndef SRC_COMPONENTS_MESHUTILITIES_EMBEDDEDSURFACEGENERATOR_HPP_
#define SRC_COMPONENTS_MESHUTILITIES_EMBEDDEDSURFACEGENERATOR_HPP_



#include "dataRepository/Group.hpp"
#include "meshUtilities/MeshGeneratorBase.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const nVector          = "unitNormalVector";
string const planeCenter      = "planeCenter";
string const meshBodyName     = "meshName";
}
}

/**
 * @class InternalWellGenerator
 *
 * This class processes the data of a single well from the XML and generates the well geometry
 */
class EmbeddedSurfaceGenerator : public MeshGeneratorBase
{
public:

  EmbeddedSurfaceGenerator( const std::string& name,
                         Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~EmbeddedSurfaceGenerator() override;

  /**
   * @return the name of this type in the catalog
   */
  static string CatalogName() { return "EmbeddedSurface"; }

  /// not implemented
  virtual void GenerateElementRegions( DomainPartition & GEOSX_UNUSED_ARG( domain ) ) override {}

  virtual Group * CreateChild( string const & childKey,
                                      string const & childName ) override;

  /**
   * @brief main function of this class: processes the well input and creates the globla well topology
   * @param domain the physical domain object
   */
  virtual void GenerateMesh( DomainPartition * const domain ) override;

  /// not implemented
  virtual void GetElemToNodesRelationInBox ( std::string const & GEOSX_UNUSED_ARG( elementType ),
                                             int const * GEOSX_UNUSED_ARG( index ),
                                             int const & GEOSX_UNUSED_ARG( iEle ),
                                             int * GEOSX_UNUSED_ARG( nodeIDInBox ),
                                             int const GEOSX_UNUSED_ARG( size )) override {}

  /// not implemented
  virtual void RemapMesh ( dataRepository::Group * const GEOSX_UNUSED_ARG( domain ) ) override {}

  // getters for element data

  /**
   * @brief Getter for the global number of well elements
   * @return the global number of elements
   */
  globalIndex GetNumElements() const { return m_numElems; }

  /**
   * @brief Getter for the global indices of the well nodes nodes connected to each element
   * @return list providing the global index of the well nodes for each well element
   */
  array2d<globalIndex const> const & GetElemToNodesMap() const { return m_elemToNodesMap; }

  /**
   * @brief Getter for the volume of the well elements
   * @return list of volumes of the well elements
   */
  array1d<real64 const> const & GetElemVolume() const { return m_elemVolume; }


  // getters for node data

  /**
   * @brief Getter for the global number of well nodes
   * @return the global number of nodes
   */
  globalIndex GetNumNodes() const { return m_numNodes; }

  /**
   * @brief Getter for the physical location of the centers of well elements
   * @return list of center locations of the well elements
   */
  array1d<R1Tensor const> const & GetNodeCoords() const { return m_nodeCoords; }

protected:

  void PostProcessInput() override final;

private:

  // XML Input

  /// normal vector
  R1Tensor      m_nVector;

  // location (point determining where the plane is)
  R1Tensor      m_planeCenter;

  string        m_meshBodyName;

  // Geometry of the embedded surface (later passed to the EmbeddedSurfaceSubregion)

  // well element data

  /// Global number of well elements
  globalIndex          m_numElems;

  /// Connectivity between elements and nodes
  array2d<globalIndex> m_elemToNodesMap;

  /// Volume of well elements
  array1d<real64> m_elemVolume;

  // Surface node data

  /// Number of nodes per well element
  globalIndex          m_numNodesPerElem;

  /// Global number of  nodes
  globalIndex          m_numNodes;

  /// Physical location of the nodes
  array1d<R1Tensor>    m_nodeCoords;

  // Auxiliary data

  // Number of physical dimensions
  const int m_nDims;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_MESHUTILITIES_EMBEDDEDSURFACEGENERATOR_HPP_ */
