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

/*
 * @file InternalWellGenerator.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_WELLS_INTERNALWELLGENERATOR_HPP_
#define GEOSX_CORECOMPONENTS_WELLS_INTERNALWELLGENERATOR_HPP_

#include "dataRepository/Group.hpp"
#include "WellGeneratorBase.hpp"

namespace geosx
{


/**
 * @class InternalWellGenerator
 *
 * This class processes the data of a single well from the XML and generates the well geometry
 */  
class InternalWellGenerator : public WellGeneratorBase
{
public:

  InternalWellGenerator( const std::string& name,
                         Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~InternalWellGenerator() override;

  /**
   * @return the name of this type in the catalog
   */  
  static string CatalogName() { return "InternalWell"; }

  /// not implemented
  virtual void GenerateElementRegions( DomainPartition& domain ) override {}

  /**
   * @brief main function of this class: processes the well input and creates the globla well topology
   * @param domain the physical domain object
   */  
  virtual void GenerateMesh( DomainPartition * const domain ) override;

  /// not implemented 
  virtual void GetElemToNodesRelationInBox ( const std::string& elementType,
                                             const int index[],
                                             const int& iEle,
                                             int nodeIDInBox[],
                                             const int size) override {}

  /// not implemented
  virtual void RemapMesh ( dataRepository::Group * const domain ) override {}

protected:

  void PostProcessInput() override final;

private:

  /**
   * @brief Map each polyline node to the polyline segment(s) it is connected to
   */
  void ConstructPolylineNodeToSegmentMap();

  /**
   * @brief Find the head node of the well (i.e., top node of the polyline)
   */
  void FindPolylineHeadNodeIndex();

  /**
   * @brief Discretize the polyline by placing well elements
   */  
  void DiscretizePolyline();

  /**
   * @brief Map each perforation to a well element
   */  
  void ConnectPerforationsToWellElements();
 
  /**
   * @brief Merge perforations on the elements with multiple perforations
   */
  void MergePerforations();

  /**
   * @brief At a given node, find the next segment going in the direction of the bottom of the well
   */  
  globalIndex GetNextSegmentIndex( globalIndex topSegId,
                                   globalIndex currentNodeId ) const;

  void DebugWellGeometry() const;


  private:
  // Auxiliary data


};
}

#endif /* GEOSX_CORECOMPONENTS_WELLS_INTERNALWELLGENERATOR_HPP_ */
