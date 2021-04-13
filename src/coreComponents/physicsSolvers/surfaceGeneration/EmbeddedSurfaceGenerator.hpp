/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file EmbeddedSurfaceGenerator.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SURFACEGENERATION_EMBEDDEDSURFACEGENERATOR_HPP_
#define GEOSX_PHYSICSSOLVERS_SURFACEGENERATION_EMBEDDEDSURFACEGENERATOR_HPP_

#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "mesh/DomainPartition.hpp"


namespace geosx
{

//struct ModifiedObjectLists
//{
//  std::set<localIndex> newNodes;
//  std::set<localIndex> newEdges;
//  std::set<localIndex> newFaces;
//  map< std::pair<localIndex,localIndex>, std::set<localIndex> > newElements;
//  map< std::pair<localIndex,localIndex>, std::set<localIndex> > modifiedElements;
//};


class SpatialPartition;

class NodeManager;
class EdgeManager;
class FaceManager;
class ExternalFaceManager;
class ElementRegionManager;
class ElementRegionBase;

/**
 * @class SurfaceGenerator
 *
 * This solver manages the mesh topology splitting methods.
 *
 */
class EmbeddedSurfaceGenerator : public SolverBase
{
public:
  EmbeddedSurfaceGenerator( const string & name,
                            Group * const parent );
  ~EmbeddedSurfaceGenerator() override;


  static string catalogName() { return "EmbeddedSurfaceGenerator"; }

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const GEOSX_UNUSED_PARAM( eventCounter ),
                        real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                        DomainPartition & domain ) override
  {
    solverStep( time_n, dt, cycleNumber, domain );
    return false;
  }

  /**
   * @brief xxx
   * @param[in] ...
   * @param[in] ...
   * @param[in] ...
   * @return ...
   *
   */
  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override;

  /**@}*/

protected:

  virtual void initializePostSubGroups() override final;

  virtual void initializePostInitialConditionsPreSubGroups() override final;

private:

  void addToFractureStencil( DomainPartition & domain );

  void setGlobalIndices( ElementRegionManager & elemManager,
                         EmbeddedSurfaceSubRegion & embeddedSurfaceSubregion );

  void addEmbeddedElementsToSets( ElementRegionManager const & elemManager,
                                  EmbeddedSurfaceSubRegion & embeddedSurfaceSubregion );

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * solidMaterialNameString() {return "solidMaterialNames"; }
    constexpr static char const * fractureRegionNameString() {return "fractureRegion"; }
    //TODO: rock toughness should be a material parameter, and we need to make rock toughness to KIC a constitutive
    // relation.
    constexpr static char const * rockToughnessString() {return "rockToughness"; }
  };

  // solid solver name
  array1d< string > m_solidMaterialNames;
  // fracture region name
  string m_fractureRegionName;
};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_SURFACEGENERATION_EMBEDDEDSURFACEGENERATOR_HPP_ */
