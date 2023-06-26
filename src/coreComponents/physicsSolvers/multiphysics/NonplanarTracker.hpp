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

/**
 * @file NonplanarTracker.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_NONPLANARTRACKER_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_NONPLANARTRACKER_HPP_

#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "PhaseFieldFractureSolver.hpp"

namespace geos
{

class NonplanarTracker : public PhaseFieldFractureSolver
{
public:

  using Base = PhaseFieldFractureSolver;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  NonplanarTracker( const string & name,
                    Group * const parent );

  ~NonplanarTracker() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "NonplanarTracker";
  }

  void buildBaseToPatchMaps( const MeshLevel & base, 
                             const MeshLevel & patch );

  virtual void registerDataOnMesh( dataRepository::Group & meshBodies ) override;                             

  void initializeFracturedElements( MeshLevel & base, MeshLevel & patch );      

  void initializeCrackFront( MeshLevel & base );      

  void cutDamagedElements( MeshLevel & base,
                           MeshLevel const & patch );

  virtual real64 sequentiallyCoupledSolverStep( real64 const & time_n,
                                                real64 const & dt,
                                                int const cycleNumber,
                                                DomainPartition & domain ) override;                                                                 

private:

  //set of base elements in the crack front
  SortedArray<localIndex> m_crackFront;

  //indicator of insertion of new fracture elements
  bool m_addedFractureElements;

  //map from base to patch elements
  map<globalIndex, set<globalIndex>> m_baseToPatchElementRelation;

  //map from patch to base elements
  map<globalIndex, globalIndex> m_patchToBaseElementRelation;

  //map from base edges to patch nodes
  map<globalIndex, vector<globalIndex>> m_baseEdgeToPatchNodeRelation;

  //set of fully fractured elements of the base mesh  
  SortedArrayView< localIndex const > const m_fracturedElements

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_NONPLANARTRACKER_HPP_ */