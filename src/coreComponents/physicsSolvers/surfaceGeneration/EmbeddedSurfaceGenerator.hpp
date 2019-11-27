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

#ifndef SRC_COMPONENTS_SURFACEGENERATION_EMBEDDEDSURFACEGENERATOR_HPP_
#define SRC_COMPONENTS_SURFACEGENERATION_EMBEDDEDSURFACEGENERATOR_HPP_

#include "mpiCommunications/NeighborCommunicator.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "managers/DomainPartition.hpp"


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
  EmbeddedSurfaceGenerator( const std::string& name,
                    Group * const parent );
  ~EmbeddedSurfaceGenerator() override;


  static string CatalogName() { return "EmbeddedSurfaceGenerator"; }

  virtual void RegisterDataOnMesh( Group * const MeshBody ) override final;

  virtual void Execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const GEOSX_UNUSED_ARG( eventCounter ),
                        real64 const  GEOSX_UNUSED_ARG( eventProgress ),
                        dataRepository::Group * domain ) override
  {
    SolverStep( time_n, dt, cycleNumber, domain->group_cast<DomainPartition*>());
  }

  /**
       * @brief xxx
       * @param[in] ...
       * @param[in] ...
       * @param[in] ...
       * @return ...
       *
       */
  virtual real64 SolverStep( real64 const& time_n,
                             real64 const& dt,
                             integer const cycleNumber,
                             DomainPartition * domain ) override;

  /**@}*/

protected:

  /**
       * @brief xxx
       * @param[in] ...
       * @param[in] ...
       * @param[in] ...
       * @return ...
       *
       */
  virtual void InitializePostSubGroups( Group * const problemManager ) override final;
  /**
         * @brief xxx
         * @param[in] ...
         * @param[in] ...
         * @param[in] ...
         * @return ...
         *
         */
  virtual void InitializePostInitialConditions_PreSubGroups( Group * const problemManager ) override final;
  virtual void postRestartInitialization( Group * const domain ) override final;

private:

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto solidMaterialNameString = "solidMaterialName";
    constexpr static auto fractureRegionNameString = "fractureRegion";
    //TODO: rock toughness should be a material parameter, and we need to make rock toughness to KIC a constitutive relation.
    constexpr static auto rockToughnessString = "rockToughness";
  }; //SurfaceGenViewKeys;

private:
  // solid solver name
  string m_solidMaterialName;
  // fracture region name
  string m_fractureRegionName;
};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_SURFACEGENERATION_EMBEDDEDSURFACEGENERATOR_HPP_ */
