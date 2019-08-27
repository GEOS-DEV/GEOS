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

#include "MPI_Communications/NeighborCommunicator.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "managers/DomainPartition.hpp"

namespace geosx
{

struct ModifiedObjectLists
{
  std::set<localIndex> newNodes;
  std::set<localIndex> newEdges;
  std::set<localIndex> newFaces;
  std::set<localIndex> modifiedNodes;
  std::set<localIndex> modifiedEdges;
  std::set<localIndex> modifiedFaces;
  map< std::pair<localIndex,localIndex>, std::set<localIndex> > newElements;
  map< std::pair<localIndex,localIndex>, std::set<localIndex> > modifiedElements;

  void clearNewFromModified();

  void insert( ModifiedObjectLists const & lists );
};


class SpatialPartition;

class NodeManager;
class EdgeManager;
class FaceManager;
class ExternalFaceManager;
class ElementRegionManager;
class ElementRegion;

/**
 * @class EmbeddedSurfaceGenerator
 *
 * This solver manages the methods to create embedded surface elements.
 *
 */
class EmbeddedSurfaceGenerator : public SolverBase
{
public:
  EmbeddedSurfaceGenerator( const std::string& name,
                    ManagedGroup * const parent );
  ~EmbeddedSurfaceGenerator() override;


  static string CatalogName() { return "EmbeddedSurfaceGenerator"; }

  virtual void RegisterDataOnMesh( ManagedGroup * const MeshBody ) override final;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void Execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::ManagedGroup * domain ) override
  {
    SolverStep( time_n, dt, cycleNumber, domain->group_cast<DomainPartition*>());
  }

  virtual real64 SolverStep( real64 const& time_n,
                             real64 const& dt,
                             integer const cycleNumber,
                             DomainPartition * domain ) override;

  /**@}*/

protected:
  virtual void InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager ) override final;
  virtual void postRestartInitialization( ManagedGroup * const domain ) override final;

private:

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto ruptureStateString = "ruptureState";
    constexpr static auto failCriterionString = "failCriterion";
    constexpr static auto degreeFromCrackString = "degreeFromCrack";
    constexpr static auto fractureRegionNameString = "fractureRegion";
  }; //SurfaceGenViewKeys;

private:
  /// choice of failure criterion
  integer m_failCriterion=1;

  /// set of separable faces
  localIndex_set m_separableFaceSet;

  /// name of the element region to place all new fractures
  string m_fractureRegionName;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_SURFACEGENERATION_EMBEDDEDSURFACEGENERATOR_HPP_ */
